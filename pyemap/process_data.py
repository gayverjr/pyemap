# PyeMap: A python package for automatic identification of electron and hole transfer pathways in proteins.
# Copyright(C) 2017-2020 Ruslan Tazhigulov, James Gayvert, Ksenia Bravaya (Boston University, USA)
"""Processes parsed .pdb/.mmcif file, and generates a graph based on user selected options.

Collects requested residues and constructs a distance matrix, from which the graph is generated. Edges are filtered out based
on user specifications. Finally, the selected criteria for surface exposure is used to classify residues as buried or exposed.
Results are stored in the emap object which was passed in.

"""

import sys
import Bio.PDB
import networkx as nx
import numpy as np
from Bio.PDB.DSSP import DSSP
from Bio.PDB.ResidueDepth import get_surface, residue_depth
from scipy.spatial import distance_matrix
import warnings
from .data import res_name_to_char, side_chain_atoms


# Monkey patches detach self to save original ID upon re-assignment to custom residue
def detach_parent(self):
    if self.parent:
        self.original_id = self.parent.full_id
    self.parent = None

Bio.PDB.Atom.Atom.detach_parent = detach_parent
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


def calculate_residue_depth(filename, aromatic_residues, rd_cutoff):
    """Returns a list of surface exposed residues as determined by residue depth.

    Parameters
    ----------
    filename: str
        Name of pdb file to be analyzed
    aromatic_residues: list of :class:`Bio.PDB.Residue.Residue`
        residues included in the analysis
    Returns
    -------
    surface_exposed_res: list of str
        List of residue names corresponding to the surface exposed residues

    """
    try:
        surface = get_surface(model)
        cutoff = rd_cutoff
        surface_exposed_res = []
        for residue in aromatic_residues:
            depth = residue_depth(residue, surface)
            if (depth <= cutoff):
                surface_exposed_res.append(residue.node_label)
        return surface_exposed_res
    except Exception as e:
        warnings.warn(
            "Unable to calculate residue depth. Check that MSMS is installed.",
            RuntimeWarning,
            stacklevel=2)
        return []


def calculate_asa(filename, model, node_list, asa_cutoff):
    """Returns a list of surface exposed residues as determined by relative solvent accessibility.

    Only standard protein residues are currently supported. Non-protein and user specified custom residues cannot be
    classified as surface exposed using this criteria.

    Parameters
    ---------
    filename: str
        Name of pdb file to be analyzed
    node_list : list of str
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
    cutoff = asa_cutoff
    surface_exposed_res = []
    try:
        dssp = DSSP(model, filename, acc_array="Wilke")
        for key in dssp.keys():
            goal_str = dssp[key][1] + str(key[1][1]) + "(" + str(key[0]) + ")"
            if goal_str in node_list and dssp[key][3] >= cutoff:
                surface_exposed_res.append(goal_str)
    except Exception as e:
        warnings.warn("Unable to calculate solvent accessibility. Check that DSSP is installed.",
                      RuntimeWarning,
                      stacklevel=2)
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
    if res.resname in side_chain_atoms:
        atom_list = []
        sca = side_chain_atoms[res.resname]
        for atm in res:
            if atm.name in sca:
                atom_list.append(atm)
        return atom_list
    else:
        return res.get_atoms()


def get_full_atom_distance_matrix(residues):
    atm_d = []
    atoms_per_res = []
    for res in residues:
        res.get_full_id()
        atoms = get_atom_list(res)
        num_atoms = 0
        for cur_atom in atoms:
            num_atoms+=1
            atm_d.append(np.array([cur_atom.coord[0], cur_atom.coord[1], cur_atom.coord[2]]))
        atoms_per_res.append(num_atoms)
    return distance_matrix(atm_d, atm_d), atoms_per_res


def closest_atom_dmatrix(residues):
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
    """
    dmat, atoms_per_res = get_full_atom_distance_matrix(residues)
    distance_matrix = np.zeros((len(residues), len(residues)))
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
            slice3_idx = slice4_idx
        slice1_idx = slice2_idx
    return distance_matrix

def com_dmatrix(residues):
    """Constructs distance matrix based on distances between centers of mass.

    Parameters
    ----------
    residues: list of :class:`Bio.PDB.Residue.Residue`
        List of residues to be included in analysis
    Returns
    -------
    distance_matrix: numpy array of :class:`Bio.PDB.Residue.Residue`
        Distance matrix of residues

    """
    com_d = []
    # Compute COMs (x0, y0, z0) of side-chains
    for res in residues:
        res.get_full_id()
        atoms = get_atom_list(res)
        x_wsum, y_wsum, z_wsum, mass_sum = 0.0, 0.0, 0.0, 0.0
        for cur_atom in atoms:
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
    for res in standard_residue_list:
        res_letter = res_name_to_char.get(res.resname)
        chain = res.full_id[2]
        resnum = res.full_id[3][1]
        res.node_label = res_letter + str(resnum) + "(" + chain + ")"
        atm_ids = []
        atm_names = []
        valid = False
        for atm in res:
            sca = side_chain_atoms[res.resname]
            if atm.name in sca:
                valid = True
                break
        if valid:
            res_list.append(res)
        else:
            def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
                return '%s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
            warnings.formatwarning = warning_on_one_line
            warnings.warn("The record for residue " + res.node_label +
                          " did not contain any of the following side chain atoms: " + str(sca) +
                          " and therefore is not included in the graph.",
                          RuntimeWarning,
                          stacklevel=2)
    return res_list


def create_user_res(serial_list, all_atoms, chain_selected, used_atoms, user_res_names):
    """Creates a customized BioPython Residue object corresponding to a user specified residue.

    Parameters
    ----------
    serial_list: list of int
        List of atom serial numbers included in residue
    all_atoms: iterator of :class:`Bio.PDB.Atom.Atom`
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
                message = "Invalid atom serial number range. Atom " + str(atm.serial_number) + " is in chain " + str(
                    atm.parent.parent.id) + "."
                raise ValueError(message)
            atm.get_full_id()
            if not source_res:
                source_res = atm.parent
            atm.get_full_id()
            atm_copy = atm.copy()
            atm_list.append(atm_copy)
    source_res.get_full_id()
    user_res = source_res.copy()
    for atm in source_res:
        if atm not in atm_list:
            user_res.detach_child(atm.id)
    k = 1
    name = "CUST"
    name += "-"
    while name + str(k) in user_res_names:
        k += 1
    user_res_names.append(name + str(k))
    user_res.resname = name + str(k)
    user_res.node_label = user_res.resname
    return user_res


def get_standard_residues(all_residues, chain_list, include_residues):
    """Returns list of standard aromatic residues.

    Runs through every residue in the structure, keeping only the TRP, TYR
    (optionally PHE and HIS if specified by user) residues on the chains
    that were chosen by the user.

    Parameters
    ----------
    all_residues: iterator of :class:`Bio.PDB.Residue.Residue`
        List of every residue in structure
    chain_list: list of str
        Chains to be included in analysis
    include_residues: list of str
        Included amino acids specified by 3 letter code

    Returns
    -------
    residue_list: list of :class:`Bio.PDB.Residue.Residue`
        Residues to be included in analysis

    """
    residue_list = []
    for res in all_residues:
        if res.resname in include_residues and res.parent.id in chain_list:
            res.get_full_id()
            arom_res = res.copy()
            arom_res.get_full_id()
            residue_list.append(arom_res)
    residue_list = process_standard_residues(residue_list)
    return residue_list


def get_user_residues(custom, all_atoms, chain_selected, used_atoms):
    """Generates customized Bio.PDB.Residue.Residue objects from a custom atom string.

    Users may not choose atoms that are already part of standard residues or eta moieties, nor those that are not
    part of the selected chains.

    Parameters
    ----------
    custom: str
        Specified by user to select atoms for custom residues
    all_atoms: iterator of :class:`Bio.PDB.Atom.Atom`
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
                new_res = create_user_res(serial_number_list, all_atoms, chain_selected, used_atoms, user_res_names)
            else:
                raise SyntaxError("Invalid atom serial number range. See the manual for proper syntax.")
            res_list.append(new_res)
            for atm in new_res:
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

def filter_edges(G,coef_alpha, exp_beta, r_offset, distance_cutoff, percentile, max_degree, include_residues,eta_moieties):
    '''Applies distance based filters to edges, removing those edges which do not fit criteria.

    Parameters
    ----------
    G: :class:`networkx.graph`
        protein graph
    distance_cutoff,percentile,max_degree: float
        Parameters that determine which edges are kept.
    '''
    G2 = G.copy()
    edge_types = [[[] for i in range(len(include_residues))] for j in range(len(include_residues))]
    edge_dist_dict = {}
    letter_codes = [res_name_to_char[res.upper()] for res in include_residues]
    minval = min(dict(G.edges).items(), key=lambda x: x[1]['weight'])[1]['weight']
    # should never happen, but just in case
    if minval == 0:
        minval = 1
    remove_edges = []
    # impose hard cutoff, collect other edge weights
    for u, v, d in G.edges(data=True):
        if d['weight'] > distance_cutoff:
            remove_edges.append((u,v))
        else:
            d['distance'] = d['weight']
            d['weight'] = pathways_model(d['weight'],coef_alpha, exp_beta, r_offset)
            d['len'] = d['weight'] / minval  # scaling factor for prettier graphs
            if u not in eta_moieties and v not in eta_moieties:
                idx1 = letter_codes.index(u[0])
                idx2 = letter_codes.index(v[0])
                edge_types[idx1][idx2].append(d['distance'])
    G.remove_edges_from(remove_edges)
    '''
    remove_edges = []
    # collect percentiles for each edge type
    for i in range(0,len(edge_types)):
        for j in range(i,len(edge_types)):
            key = letter_codes[i]+letter_codes[j]
            if(len(edge_types[i][j])) > 1:
                value = np.percentile(edge_types[i][j],percentile)
            else:
                value = 0.0
            edge_dist_dict[key] = value
    default_cutoff = np.ma.masked_equal(list(edge_dist_dict.values()), 0.0, copy=False).max()
    # remove edges which don't fit criteria
    for u, v, d in G.edges(data=True):
        if u in eta_moieties or v in eta_moieties:
            cutoff = default_cutoff
        else:
            idx1 = letter_codes.index(u[0])
            idx2 = letter_codes.index(v[0])
            key = letter_codes[idx1]+letter_codes[idx2]
            if key not in edge_dist_dict:
                key = letter_codes[idx2]+letter_codes[idx1]
            cutoff = edge_dist_dict[key]
        if d['distance'] > cutoff:
            remove_edges.append((u,v))
    G.remove_edges_from(remove_edges)
    '''
    remove_edges = []
    for node in G.nodes:
        if G.degree(node) > max_degree:
            for edge in sorted(list(G.edges(node)), key=lambda x: G.edges[x]['distance'])[max_degree:]:
                if edge not in remove_edges and edge[::-1] not in remove_edges:
                    remove_edges.append(edge)
    remove_edges = sorted(remove_edges, key=lambda x: G.edges[x]['distance'])
    for u,v in remove_edges:
        if G.degree(u) > max_degree or G.degree(v) > max_degree:
            G.remove_edge(u,v)
    

def create_graph(dmatrix,node_labels, coef_alpha, exp_beta, r_offset, distance_cutoff,percentile,max_degree,eta_moieties,include_residues):
    """Constructs the graph from the distance matrix and node labels.

    Parameters
    ----------
    dmatrix: numpy.array of float
        Distance matrix of aromatic residues.
    node_label: list of str
        Labels for residues in the graph.
    residue_numbers:
        res numbers
    distance_cutoff,percentile, max_degree: float
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
    G = nx.from_numpy_matrix(dmatrix)
    G = nx.relabel_nodes(G, node_labels)
    filter_edges(G,coef_alpha, exp_beta, r_offset,distance_cutoff,percentile,max_degree,include_residues,eta_moieties)
    for name_node in G.nodes():
        G.nodes[name_node]['style'] = 'filled'
        G.nodes[name_node]['fontname'] = 'Helvetica-Bold'
        G.nodes[name_node]['fontsize'] = 32
        G.nodes[name_node]['shape'] = "oval"
        G.nodes[name_node]['margin'] = '0.04'
        G.nodes[name_node]['fontcolor'] = "#000000"
        G.nodes[name_node]['color'] = '#708090'
        G.nodes[name_node]['penwidth'] = 2.0
        if (name_node[1].isdigit()):
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
                    G.nodes[name_node]['fillcolor'] = '#708090'
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
            chains=None,
            eta_moieties=None,
            dist_def=1,
            sdef=1,
            include_residues=["TYR", "TRP"],
            custom="",
            distance_cutoff=20,
            percent_edges=100,
            coef_alpha=1.0,
            max_degree = 4,
            exp_beta=2.3,
            r_offset=0.0,
            rd_thresh=3.03,
            asa_thresh=0.2):
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
    include_residues: list of str
        Included amino acids specified by 3 letter code
    custom: str, optional
        Custom atom string specified by user
    distance_cutoff: float
         Defines a pure distance threshold. PyeMap will only keep edges with distances less than or equal distance_cutoff.
    percent_edges: float
        Specifies a percentage of the shortest edges per edge type to keep.
    max_degree: int
        Maximum degree of any vertex.
    coef_alpha,exp_beta,r_offset: float, optional
        Penalty function parameters.
    Raises
    ------
    RuntimeError:
        Not enough residues to construct a graph

    """
    emap._reset_process()
    pdb_file = emap.file_path
    if chains == None:
        chains = [emap.chains[0]]
    if eta_moieties == None:
        eta_moieties = []
        for resname,moiety in emap.eta_moieties.items():
            if moiety.get_full_id()[2] in chains:
                eta_moieties.append(resname)
    else:
        for resname in eta_moieties:
            if resname[resname.index("(")+1:resname.index(")")] not in chains:
                warnings.warn(str(resname) + " is not in a selected chain, so it will not be included in the analysis.")
                eta_moieties.remove(resname)
    for i, res in enumerate(include_residues):
        if res.upper() in res_name_to_char:
            include_residues[i] = res.upper()
        else:
            raise RuntimeError("Error: " + str(res) + " is not a valid 3 letter amino acid code.")
    model = emap._structure[0]
    aromatic_residues = get_standard_residues(model.get_residues(), chains, include_residues)
    for resname in eta_moieties:
        aromatic_residues.append(emap.eta_moieties[resname])
    used_atoms = []
    for res in aromatic_residues:
        for atm in res:
            used_atoms.append(atm.serial_number)
    user_residues = []
    if custom:
        user_residues = get_user_residues(custom, model.get_atoms(), chains, used_atoms)
    for res in user_residues:
        emap.user_residues[res.resname] = res
    aromatic_residues += user_residues
    node_labels = {}
    for i in range(0, len(aromatic_residues)):
        node_labels[i] = aromatic_residues[i].node_label
    if int(dist_def) == 0:
        dmatrix = com_dmatrix(aromatic_residues)
    else:
        dmatrix = closest_atom_dmatrix(aromatic_residues)
    G = create_graph(dmatrix, node_labels, coef_alpha, exp_beta, r_offset,distance_cutoff,percent_edges, int(max_degree),
                     emap.eta_moieties.keys(),include_residues)
    G.graph['pdb_id'] = emap.pdb_id
    if len(G.edges()) == 0:
        raise RuntimeError("Not enough edges to construct a graph.")
    # define surface exposed residues
    if sdef and int(sdef) == 0:
        surface_exposed_res = calculate_residue_depth(pdb_file,aromatic_residues,rd_thresh)
    elif sdef and int(sdef) == 1:
        surface_exposed_res = calculate_asa(pdb_file, model, node_labels.values(), asa_thresh)
    else:
        surface_exposed_res = []
    finish_graph(G, surface_exposed_res, chains)
    for res in aromatic_residues:
        emap._add_residue(res)
    emap._store_initial_graph(G)
    return emap
