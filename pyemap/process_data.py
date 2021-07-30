# PyeMap: A python package for automatic identification of electron and hole transfer pathways in proteins.
# Copyright(C) 2017-2020 Ruslan Tazhigulov, James Gayvert, Ksenia Bravaya (Boston University, USA)
"""Processes parsed .pdb/.mmcif file, and generates a graph based on user selected options.

Collects requested residues and constructs a distance matrix, from which the graph is generated. Edges are filtered out based
on user specifications. Finally, the selected criteria for surface exposure is used to classify residues as buried or exposed.
Results are stored in the emap object which was passed in.

"""
import math
import sys
import Bio.PDB
import networkx as nx
import numpy as np
from Bio.PDB.DSSP import DSSP
from Bio.PDB.ResidueDepth import get_surface, residue_depth
from scipy.spatial import distance_matrix
import warnings
from .data import res_name_to_char, TRP_sc, TYR_sc, PHE_sc, HIS_sc

# Monkey patches detach self to save original ID upon re-assignment to custom residue


def detach_parent(self):
    if self.parent:
        self.original_id = self.parent.full_id
    self.parent = None


Bio.PDB.Atom.Atom.detach_parent = detach_parent

# Monkey patches Atom __init__ to include a variable which tracks its original ID and assigns a unique ID


def atm_init(self, name, coord, bfactor, occupancy, altloc, fullname, serial_number, element=None):
    self.level = "A"
    # Reference to the residue
    self.parent = None
    # the atomic data
    self.name = name  # eg. CA, spaces are removed from atom name
    self.fullname = fullname  # e.g. " CA ", spaces included
    self.coord = coord
    self.bfactor = bfactor
    self.occupancy = occupancy
    self.altloc = altloc
    # (structure id, model id, chain id, residue id, atom id)
    self.full_id = None
    # id of atom is the atom name (e.g. "CA")
    self.id = name + str(serial_number)
    self.disordered_flag = 0
    self.anisou_array = None
    self.original_id = None
    self.siguij_array = None
    self.sigatm_array = None
    self.serial_number = serial_number
    # Dictionary that keeps additional properties
    self.xtra = {}
    assert not element or element == element.upper(), element
    self.element = self._assign_element(element)
    self.mass = self._assign_atom_mass()
    # For atom sorting (protein backbone atoms first)
    self._sorting_keys = {'N': 0, 'CA': 1, 'C': 2, 'O': 3}


Bio.PDB.Atom.Atom.__init__ = atm_init

# monkey patches mangling of disordered atoms


def get_unpacked_list(self):
    """
     Returns all atoms from the residue,
     in case of disordered, keep only first alt loc and remove the alt-loc tag
    """
    atom_list = self.get_list()
    undisordered_atom_list = []
    for atom in atom_list:
        if atom.is_disordered():
            atom.altloc = " "
            undisordered_atom_list.append(atom)
        else:
            undisordered_atom_list.append(atom)
    return undisordered_atom_list


Bio.PDB.Residue.Residue.get_unpacked_list = get_unpacked_list


def pathways_model(dist, coef_alpha, exp_beta, r_offset):
    """Applies penalty function parameters and returns score.

        :math:`\\epsilon =\\alpha \\exp(-\\beta(R-R_{offset}))`

    Parameters
    ----------
    dist: float
        Actual distance in angstroms
    coef_alpha,exp_beta,r_offset:float
        Penalty function parameters

    Returns
    -------
    mod_penalty: float
    
    """
    penalty = coef_alpha * np.exp(-exp_beta * (dist - r_offset))
    mod_penalty = -np.log10(penalty)
    return mod_penalty


def calculate_residue_depth(aromatic_residues, model):
    """Returns a list of surface exposed residues as determined by residue depth.

    Parameters
    ----------
    aromatic_residues: list of :class:`Bio.PDB.Residue.Residue`
        residues included in the analysis
    model: :class:`Bio.PDB.Model.Model`
        BioPython object containing a model of a PDB or MMCIF file

    Returns
    -------
    surface_exposed_res: list of str
        List of residue names corresponding to the surface exposed residues

    """
    try:
        surface = get_surface(model)
        cutoff = 3.03
        surface_exposed_res = []
        for residue in aromatic_residues:
            depth = residue_depth(residue, surface)
            if (depth <= cutoff):
                surface_exposed_res.append(residue.node_label)
        return surface_exposed_res
    except Exception as e:
        warnings.warn("Unable to calculate residue depth. Check that MSMS is installed. Please note that MSMS is not compatible with MacOS Catalina.", RuntimeWarning,stacklevel=2)
        return []


def calculate_asa(model, filename, node_list):
    """Returns a list of surface exposed residues as determined by relative solvent accessibility.

    Only standard protein residues are currently supported. Non-protein and user specified custom residues cannot be
    classified as surface exposed using this criteria.

    Parameters
    ---------
    model: :class:`Bio.PDB.Model.Model`
        Model which contains chains and residues of protein strucutre
    filename: str
        Name of pdb file to be analyzed
    AROM_LIST : list of str
        List containing which standard residues are included in analysis
    chain_list: list of str
        Chains are included in analysis

    Notes
    -----
    The relative accessible surface area (RSA) of each residue is calculated using the Bio.PDB.DSSP module. A residue
    with an RSA value of 0.05 or higher is classified as surface exposed.

    References
    ---------
    Tien, M. Z.; Meyer, A. G.; Sydykova, D. K.; Spielman, S. J.; Wilke, C. O. PLoS ONE 2013, 8 (11).
        Reference for relative solvent accessibility cutoff of 0.05, and for MaxASA values

    """
    cutoff = .05
    surface_exposed_res = []
    try:
        dssp = DSSP(model, filename, acc_array="Wilke")
        keys = list(dssp.keys())
        for key in keys:
            goal_str = dssp[key][1] + str(key[1][1]) + "(" + str(key[0]) + ")"
            if goal_str in node_list and dssp[key][3] >= cutoff:
                surface_exposed_res.append(goal_str)
    except Exception as e:
        warnings.warn("Unable to calculate solvent accessibility. Check that DSSP is installed.", RuntimeWarning,stacklevel=2)
    return surface_exposed_res


def dist(x, y):
    """Returns Euclidean distance between two 3 dimensional points.

    Parameters
    ----------
    x: 3D numpy array of float
        Point 1
    y: 3D numpy array of float
        Point 2

    Returns
    -------
    float:
        Euclidean distance between two points

    """
    return np.sqrt(np.sum((x - y)**2))

def get_atom_list(res):
    """ Recovers side chain atoms only of standard aromatic residues, returns all atoms for other residues.
    
    Parameters
    -----------
    res: :class:`Bio.PDB.Residue.Residue`
        A BioPython residue object
    
    Returns
    -------
    atom_list: list of :class:`Bio.PDB.Atom.Atom`
        List of atoms to be used in distance matrix calculation
    """
    atom_list = []
    for atm in res.get_atoms():
        if res.resname == "TRP" and atm.name in TRP_sc:
            atom_list.append(atm)
        elif res.resname == "TYR" and atm.name in TYR_sc:
            atom_list.append(atm)
        elif res.resname == "PHE" and atm.name in PHE_sc:
            atom_list.append(atm)
        elif res.resname == "HIS" and atm.name in HIS_sc:
            atom_list.append(atm)
        elif res.resname not in ["TRP","TYR","PHE","HIS"]:
            atom_list.append(atm)
    return atom_list

def get_full_atom_distance_matrix(residues):
    com_d = []
    atoms_per_res = []
    for i in range(len(residues)):
        res = residues[i]
        res.get_full_id()
        atm_list = get_atom_list(res)
        atoms_per_res.append(len(atm_list))
        for j in range(len(atm_list)):
            cur_atom = atm_list[j]
            com_d.append(np.array([cur_atom.coord[0], cur_atom.coord[1], cur_atom.coord[2]]))
    return distance_matrix(com_d, com_d), atoms_per_res


def get_com_distance_matrix(residues):
    com_d = []
    # Compute COMs (x0, y0, z0) of side-chains
    for i in range(len(residues)):
        res = residues[i]
        res.get_full_id()
        atm_list = get_atom_list(res)
        x_wsum, y_wsum, z_wsum, mass_sum = 0.0, 0.0, 0.0, 0.0
        for j in range(len(atm_list)):
            cur_atom = atm_list[j]
            mass = cur_atom.mass
            x = cur_atom.coord[0]
            y = cur_atom.coord[1]
            z = cur_atom.coord[2]
            x_wsum += x * mass
            y_wsum += y * mass
            z_wsum += z * mass
            mass_sum += mass
        com_x = x_wsum / mass_sum
        com_y = y_wsum / mass_sum
        com_z = z_wsum / mass_sum
        com_d.append(np.array([com_x, com_y, com_z]))
    return distance_matrix(com_d, com_d)


def closest_atom_dmatrix(residues, coef_alpha, exp_beta, r_offset):
    """Constructs distance matrix based on closest atom distance.
    
    Parameters
    ----------
    residues: list of :class:`Bio.PDB.Residue.Residue`
        List of BioPython residues
    coef_alpha,exp_beta,r_offset:float
        Penalty funciton parameters
    Returns
    -------
    distance_matrix: numpy.array of float
        Distance matrix of residues
    pathways_matrix: numpy.array of float
        Modified distance matrix of residues using scores determined by penalty function parameters
    """
    dmat, atoms_per_res = get_full_atom_distance_matrix(residues)
    distance_matrix = np.zeros((len(residues), len(residues)))
    pathways_matrix = np.zeros((len(residues), len(residues)))
    slice1_idx = 0
    for i in range(0, len(residues)):
        slice2_idx = slice1_idx + atoms_per_res[i]
        slice3_idx = slice2_idx
        for j in range(i + 1, len(residues)):
            slice4_idx = slice3_idx + atoms_per_res[j]
            my_slice = dmat[slice1_idx:slice2_idx, slice3_idx:slice4_idx]
            min_val = np.min(my_slice)
            distance_matrix[i][j] = min_val
            distance_matrix[j][i] = min_val
            pathways_matrix[i][j] = pathways_model(min_val, coef_alpha, exp_beta, r_offset)
            pathways_matrix[j][i] = pathways_matrix[i][j]
            slice3_idx = slice4_idx
        slice1_idx = slice2_idx
    return distance_matrix, pathways_matrix


def com_dmatrix(residues, coef_alpha, exp_beta, r_offset):
    """Constructs distance matrix based on distances between centers of mass.

    Parameters
    ----------
    residues: list of :class:`Bio.PDB.Residue.Residue`
        List of residues to be included in analysis
    coef_alpha,exp_beta,r_offset:float
        Penalty funciton parameters

    Returns
    -------
    node_label: dict of int:str
        List of node labels for graph
    distance_matrix: numpy array of :class:`Bio.PDB.Residue.Residue`
        Distance matrix of residues

    """
    com_d = get_com_distance_matrix(residues)
    # calculate matrix with penalty functions
    pathways_matrix = np.zeros((len(com_d), len(com_d)))
    for i in range(len(com_d)):
        for j in range(i + 1, len(com_d)):
            dist_i_j = com_d[i][j]
            pathways_matrix[i][j] = pathways_model(dist_i_j, coef_alpha, exp_beta, r_offset)
            pathways_matrix[j][i] = pathways_matrix[i][j]
    return com_d, pathways_matrix


def process_standard_residues(standard_residue_list):
    """Generates customized Bio.PDB.Residue.Residue objects. Residues containing no side chain atoms will be removed.

    Parameters
    ----------
    standard_residue_list: list of :class:`Bio.PDB.Residue.Residue`
        List of Bio.PDB.Residue.Residue objects corresponding to protein residues

    Returns
    -------
    res_list: list of :class:`Bio.PDB.Residue.Residue`

    """
    res_list = []
    for i in range(0, len(standard_residue_list)):
        res = standard_residue_list[i]
        res_letter = res_name_to_char.get(res.resname)
        chain = res.full_id[2]
        resnum = res.full_id[3][1]
        res.node_label = res_letter + str(resnum) + "(" + chain + ")"
        atm_ids = []
        atm_names = []
        valid = False
        for atm in res.get_atoms():
            atm_ids.append(atm.id)
            atm_names.append(atm.name)
        for k in range(0, len(atm_ids)):
            if res.resname == "TRP":
                if atm_names[k] in TRP_sc:
                    valid = True
            elif res.resname == "TYR":
                if atm_names[k] in TYR_sc:
                    valid = True
                    break
            elif res.resname == "PHE":
                if atm_names[k] in PHE_sc:
                    valid = True
                    break
            elif res.resname == "HIS":
                if atm_names[k] in HIS_sc:
                    valid = True
                    break
        if valid:
            res_list.append(res)
        else:
            if res.resname == "TRP":
                sca = str(TRP_sc)
            elif res.resname == "TYR":
                sca = str(TYR_sc)
            elif res.resname == "PHE":
                sca = str(PHE_sc)
            elif res.resname == "HIS":
                sca = str(HIS_sc)
            side_chain_atms = sca
            def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
                        return '%s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
            warnings.formatwarning = warning_on_one_line
            warnings.warn("The record for residue " + res.node_label + " did not contain any of the following side chain atoms: " + side_chain_atms + " and therefore is not included in the graph.",RuntimeWarning,stacklevel=2)
    return res_list


def create_user_res(serial_list, all_atoms, chain_selected, used_atoms, user_res_names):
    """Creates a customized BioPython Residue object corresponding to a user specified residue.

    Parameters
    ----------
    serial_list: list of int
        List of atom serial numbers included in residue
    all_atoms: list of :class:`Bio.PDB.Atom.Atom`
        All atoms in protein on selected chains
    chain_selected: list of str
        Chains included in analysis
    used_atoms: list of :class:`Bio.PDB.Atom.Atom`
        Atoms already included in analysis
    user_res_names: list of str
        User residue names already included in the analysis.

    Returns
    -------
    user_res: :class:`Bio.PDB.Residue.Residue`
        Customized residue object corresponding to the atoms in serial_list

    """
    source_res = []
    atm_list = []
    for atm in all_atoms:
        if atm.serial_number in serial_list:
            if atm.serial_number in used_atoms:
                message = "Invalid atom serial number range. Atom " + str(
                    atm.serial_number) + " is already included in another residue."
                raise ValueError(message)
            if atm.parent.parent.id not in chain_selected:
                message = "Invalid atom serial number range. Atom " + str(
                    atm.serial_number) + " is in chain " + str(atm.parent.parent.id) + "."
                raise ValueError(message)
            atm.get_full_id()
            if not source_res:
                source_res = atm.parent
            atm.get_full_id()
            atm_copy = atm.copy()
            atm_list.append(atm_copy)
    source_res.get_full_id()
    user_res = source_res.copy()
    for atm in list(source_res.get_atoms()):
        user_res.detach_child(atm.id)
    for atm in atm_list:
        user_res.add(atm)
    k = 1
    name = "CUST"
    name += "-"
    while name + str(k) in user_res_names:
        k += 1
    user_res_names.append(name + str(k))
    user_res.resname = name + str(k)
    user_res.node_label = user_res.resname
    return user_res


def get_standard_residues(all_residues, chain_list, include_Trp, include_Tyr, include_Phe, include_His):
    """Returns list of standard aromatic residues.

    Runs through every residue in the structure, keeping only the TRP, TYR
    (optionally PHE and HIS if specified by user) residues on the chains
    that were chosen by the user.

    Parameters
    ----------
    all_residues: list of :class:`Bio.PDB.Residue.Residue`
        List of every residue in structure
    chain_list: list of str
        Chains to be included in analysis
    include_Trp,include_Tyr,include_Phe,include_His: bool
        True if residue type is included

    Returns
    -------
    residue_list: list of :class:`Bio.PDB.Residue.Residue`
        Residues to be included in analysis
    AROM_LIST: list of str
        Types of aromatic residues included in graph

    """
    AROM_LIST = ['TRP', 'HIS', 'PHE', 'TYR']
    if not include_Trp:
        AROM_LIST.remove('TRP')
    if not include_Tyr:
        AROM_LIST.remove('TYR')
    if not include_Phe:
        AROM_LIST.remove('PHE')
    if not include_His:
        AROM_LIST.remove('HIS')
    residue_list = []
    for res in all_residues:
        if res.resname in AROM_LIST and res.parent.id in chain_list:
            res.get_full_id()
            arom_res = res.copy()
            arom_res.get_full_id()
            residue_list.append(arom_res)
    residue_list = process_standard_residues(residue_list)
    return residue_list, AROM_LIST


def get_user_residues(custom, all_atoms, chain_selected, used_atoms):
    """Generates customized Bio.PDB.Residue.Residue objects from a custom atom string.

    Users may not choose atoms that are already part of standard residues or eta moieties, nor those that are not
    part of the selected chains.

    Parameters
    ----------
    custom: str
        Specified by user to select atoms for custom residues
    all_atoms: list of :class:`Bio.PDB.Atom.Atom`
        All atoms in protein on selected chains
    chain_selected: list of str
        Chains included in analysis
    used_atoms: list of :class:`Bio.PDB.Atom.Atom`
        Atoms already included in analysis

    Returns
    -------
    res_list: list of :class:`Bio.PDB.Residue.Residue`
        Custom residues specified by user

    Notes
    -----
    The custom atom string syntax is as follows:
    Custom residues are specified based on PDB atom serial number. Ranges are specified with a dash (ex. 10-20).
    Individual atoms may also be selected, separated by commas (ex. 1-10,15,17). To select multiple custom residues,
    the user encloses the atom ranges with parentheses, and then separate them with commas. In the graph and
    visualization, these residues are named CUS-1, CUS-2 etc. based on the order they were specified.

    """
    user_res_names = []
    custom = custom.strip()
    if custom:
        res_list = []
        l1 = custom.split("),(")
        for a in l1:
            serial_number_list = []
            # get rid of the (
            atm_str = a[:]
            atm_str = atm_str.replace("(", '')
            atm_str = atm_str.replace(")", '')
            l2 = atm_str.split(",")
            for atm in l2:
                if "-" in atm:
                    index = atm.index("-")
                    start_atm = int(atm[:index])
                    end_atm = atm[index + 1:]
                    end_atm = end_atm.replace(")", "")
                    end_atm = int(end_atm)
                    for i in range(start_atm, end_atm + 1):
                        serial_number_list.append(i)
                else:
                    serial_number_list.append(int(atm))
            if serial_number_list:
                new_res = create_user_res(
                    serial_number_list, all_atoms, chain_selected, used_atoms, user_res_names)
            else:
                raise SyntaxError(
                    "Invalid atom serial number range. See the manual for proper syntax.")
            res_list.append(new_res)
            for atm in new_res.get_atoms():
                used_atoms.append(atm.serial_number)
        return res_list
    return []


def finish_graph(G, surface_exposed_res, chain_list):
    """Sets surface exposed residues as boxes in graph and removes disconnected vertices.

    Parameters
    ----------
    G: :class:`networkx.graph`
        Graph object constructed from distance matrix of residues
    surface_exposed_res: list of str
        Names of surface exposed residues
    chain_list: list of str
        Names of chains included in analysis

    """
    for goal in surface_exposed_res:
        G.nodes[goal]['margin'] = '0.11'
        G.nodes[goal]['shape'] = 'box'
    # get rid of all disconnected nodes
    all_nodes = list(G.nodes())
    for node in all_nodes:
        if G[node] == {}:
            G.remove_node(node)


def filter_edges(G, G_pathways, distance_cutoff, percent_edges, num_st_dev_edges):
    '''Applies distance based filters to edges, removing those edges which do not fit criteria.

    Parameters
    ----------
    G: :class:`networkx.graph`
        graph where edge weights are pure distances
    G_pathways: :class:`networkx.graph`
        graph where edge weights are distance dependent penalty functions
    distance_cutoff,percent_edges,num_st_dev_edges: float
        Parameters that determine which edges are kept.
    '''
    # keep 99th percentile and ~20A filter on all edges
    included_edges = []
    for node in G.nodes():
        edge_length_per_node = []
        weights = []
        for neighbor in G[node]:
            weights.append(G.get_edge_data(node, neighbor)['weight'])
        weights = sorted(weights)
        thresh_index = math.ceil(len(weights) * percent_edges / 100)
        for neighbor in G[node]:
            if weights.index(G.get_edge_data(node, neighbor)['weight']) <= thresh_index and \
                    G.get_edge_data(node, neighbor)['weight'] <= distance_cutoff:
                edge_length_per_node.append(
                    G.get_edge_data(node, neighbor)['weight'])

        len_average, len_st_dev = 0.0, 0.0
        if edge_length_per_node != []:
            edge_length_per_node = np.array(
                edge_length_per_node, dtype='float64')
            len_average = np.average(edge_length_per_node)
            len_st_dev = np.std(edge_length_per_node)

        for neighbor in G[node]:
            if weights.index(G.get_edge_data(node, neighbor)['weight']) <= thresh_index and \
                    G.get_edge_data(node, neighbor)['weight'] <= distance_cutoff and \
                    np.round(G.get_edge_data(node, neighbor)['weight'], 8) <= np.round((len_average + num_st_dev_edges * len_st_dev), 8):
                included_edges.append([node, neighbor])

    excluded_edges = []
    for edge in G.edges():
        node1 = edge[0]
        node2 = edge[1]
        if ([node1, node2] not in included_edges) and ([node2, node1] not in included_edges):
            excluded_edges.append(edge)

    for edge in excluded_edges:
        node1 = edge[0]
        node2 = edge[1]
        if G.has_edge(node1, node2):
            G.remove_edge(node1, node2)
        if G_pathways.has_edge(node1, node2):
            G_pathways.remove_edge(node1, node2)


def create_graph(dmatrix, pathways_matrix, node_labels, distance_cutoff, percent_edges, num_st_dev_edges, eta_moieties):
    """Constructs the graph from the distance matrix and node labels.

    Parameters
    ----------
    dmatrix: numpy.array of float
        Distance matrix of aromatic residues.
    node_label: list of str
        Labels for residues in the graph.
    distance_cutoff,percent_edges,num_st_dev_edges: float
        Parameters that determine which edges are kept.
    eta_moieties: list of str
        Non standard residues that were automatically identified

    Returns
    -------
    G: :class:`networkx.Graph`
        Graph of aromatic residues in protein

    References
    ----------
    Gray, H. B.; Winkler, J. R. Long-Range Electron Transfer. Proc. Natl. Acad. Sci. U. S. A. 2005, 102 (10),
        3534 LP-3539.
        Reference for 20A filter on edges
    """
    np.set_printoptions(threshold=sys.maxsize)
    G = nx.from_numpy_array(dmatrix)
    minval_pathways = np.min(pathways_matrix[pathways_matrix.nonzero()])
    G_pathways = nx.from_numpy_array(pathways_matrix)
    filter_edges(G, G_pathways, distance_cutoff,
                 percent_edges, num_st_dev_edges)
    for u, v, d in G_pathways.edges(data=True):
        d['len'] = d['weight'] / minval_pathways
        d['distance'] = G[u][v]['weight']
    G = G_pathways
    G = nx.relabel_nodes(G, node_labels)
    for name_node in G.nodes():
        G.nodes[name_node]['style'] = 'filled'
        G.nodes[name_node]['fontname'] = 'Helvetica-Bold'
        G.nodes[name_node]['fontsize'] = 32
        G.nodes[name_node]['shape'] = "oval"
        G.nodes[name_node]['margin'] = '0.04'
        G.nodes[name_node]['fontcolor'] = "#000000"
        G.nodes[name_node]['color'] = '#708090'
        G.nodes[name_node]['penwidth'] = 2.0
        if(name_node[1].isdigit()):
            if name_node not in eta_moieties:
                if 'Y' == name_node[0]:
                    G.nodes[name_node]['fillcolor'] = '#96c8f0'
                elif 'W' == name_node[0]:
                    G.nodes[name_node]['fillcolor'] = '#f07878'
                elif 'F' == name_node[0]:
                    G.nodes[name_node]['fillcolor'] = '#f09664'
                elif 'H' == name_node[0]:
                    G.nodes[name_node]['fillcolor'] = '#c8f0c8'
            else:
                G.nodes[name_node]['fillcolor'] = '#FFC0CB'
        else:
            G.nodes[name_node]['fillcolor'] = '#FFC0CB'
    for edge in G.edges():
        name_node1, name_node2 = edge[0], edge[1]
        G[name_node1][name_node2]['color'] = '#778899'
        G[name_node1][name_node2]['penwidth'] = 1.5
        G[name_node1][name_node2]['style'] = 'dashed'
    return G


def process(emap,
            chains="All",
            eta_moieties="All",
            dist_def=0,
            sdef=1,
            trp=True,
            tyr=True,
            phe=False,
            his=False,
            custom="",
            distance_cutoff=20,
            percent_edges=1.0,
            num_st_dev_edges=1.0,
            coef_alpha=1.0,
            exp_beta=2.3,
            r_offset=0.0):
    """Constructs emap graph theory model based on user specs, and saves it to the emap object.

    Parameters
    ---------
    emap: :class:`~pyemap.emap`
        Object for storing state of emap analysis.
    chains: list of str
        List of strings corresponding to chains included in analysis
    eta_moieties: list of str
        List of strings corresponding to residue names of eta moieties
    dist_def: int, optional
        0 for center of mass, 1 for closest atom
    sdef: int, optional
        0 for residue depth, 1 for solvent accessibility
    trp,tyr,phe,his: bool, optional
        True if residue type is included
    custom: str, optional
        Custom atom string specified by user
    distance_cutoff: float
         Defines a pure distance threshold. PyeMap will only keep edges with distances less than or equal distance_cutoff.
    percent_edges: float
        Specifies a percentage of the shortest edges per vertex to keep.
    num_st_dev_edges: float
        For a particular vertex, only edges with length < average length + num_st_dev_edges * SD included
    coef_alpha,exp_beta,r_offset: float, optional
        Penalty function parameters.
    Raises
    ------
    RuntimeError:
        Not enough residues to construct a graph

    """
    emap._reset_process()
    if eta_moieties == "All":
        eta_moieties = emap.eta_moieties.keys()
    if chains == "All":
        chains = emap.chains
    model = emap.structure[0]
    all_residues = list(model.get_residues())
    all_atoms = list(model.get_atoms())
    aromatic_residues, AROM_LIST = get_standard_residues(
        all_residues, chains, trp, tyr, phe, his)
    for resname in eta_moieties:
        AROM_LIST.append(resname)
        aromatic_residues.append(emap.eta_moieties[resname])
    used_atoms = []
    for res in aromatic_residues:
        for atm in res.get_atoms():
            used_atoms.append(atm.serial_number)
    user_residues = []
    if custom:
        user_residues = get_user_residues(
            custom, all_atoms, chains, used_atoms)
    for res in user_residues:
        AROM_LIST.append(res.resname)
        emap.user_residues[res.resname] = res
    aromatic_residues += user_residues
    node_labels = {}
    residue_numbers = {}
    aligned_residue_numbers = {}
    for i in range(0, len(aromatic_residues)):
        node_labels[i] = aromatic_residues[i].node_label
        residue_numbers[aromatic_residues[i].node_label] = aromatic_residues[i].full_id[3][1]
        if hasattr(aromatic_residues[i],'aligned_residue_number'):
            aligned_residue_numbers[aromatic_residues[i].node_label] = aromatic_residues[i].aligned_residue_number
    if int(dist_def) == 0:
        dmatrix, pathways_matrix = com_dmatrix(aromatic_residues, coef_alpha, exp_beta, r_offset)
    else:
        dmatrix, pathways_matrix = closest_atom_dmatrix(aromatic_residues, coef_alpha, exp_beta, r_offset)
    G = create_graph(dmatrix, pathways_matrix, node_labels,
                     distance_cutoff, percent_edges, num_st_dev_edges,emap.eta_moieties.keys())
    if len(G.edges()) == 0:
        raise RuntimeError(
            "Not enough edges to construct a graph.")
    # define surface exposed residues
    if int(sdef) == 0:
        surface_exposed_res = calculate_residue_depth(aromatic_residues, model)
    else:
        pdb_file = emap.filename
        surface_exposed_res = calculate_asa(model, pdb_file, node_labels.values())
    finish_graph(G, surface_exposed_res, chains)
    for res in aromatic_residues:
        emap._add_residue(res)
    emap._store_initial_graph(G)
