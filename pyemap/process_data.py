# This Python script is a part of
# eMAP: online mapping of electron transfer channels in biomolecules
# Copyright(C) 2017-2018 Ruslan Tazhigulov, James Gayvert, Melissa Wei, Ksenia Bravaya (Boston University, USA)
"""Processes parsed .pdb/.mmcif file, and generates a graph based on user selected options.

Main module used for constructing the graph. Takes in all user specifications to gather and prepare the appropriate
residues. A distance matrix is constructed for these residues, and then the graph is generated. Finally, the selected
criteria for surface exposure is used to classify residues as buried or exposed. 

Usage
-----
Called by the views module to process the parsed file according to the user specifications, and generate the graph.

"""
import math
import os
import sys
import time
import traceback
import ast
import Bio.PDB
import networkx as nx
import numpy as np
from Bio.PDB.DSSP import DSSP
from Bio.PDB.ResidueDepth import get_surface, residue_depth
from networkx.drawing.nx_agraph import to_agraph
from scipy.spatial import distance_matrix
from .data import *
"""str: module level directory path for writing to file"""

TRP_sc = ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
TYR_sc = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH']
PHE_sc = ['CG', 'CD1', 'CD2', 'CE1', 'CZ', 'CE2']
HIS_sc = ['CG', 'ND1', 'CD2', 'CD1', 'NE1']
"""module level lists of side chain atoms for respective residues"""

def pathways_model(dist, coef_alpha, exp_beta, r_offset):
    penalty = coef_alpha * np.exp(-exp_beta * (dist - r_offset))
    mod_penalty = -np.log10(penalty)
    return mod_penalty


def pathways_model_thrspace_penalty_COM(com_d, coef_alpha, exp_beta, r_offset):
    mod_penalty_matrix = np.zeros((len(com_d), len(com_d)))
    for i in range(len(com_d)):
        for j in range(i + 1, len(com_d)):
            dist_i_j = dist(com_d[i], com_d[j])
            mod_penalty_matrix[i][j] = pathways_model(dist_i_j, coef_alpha, exp_beta, r_offset)
            mod_penalty_matrix[j][i] = mod_penalty_matrix[i][j]
    return mod_penalty_matrix

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


def calculateResidueDepth(aromatic_residues, model):
    """Returns a list of surface exposed residues as determined by residue depth.

    Parameters
    ----------
    aromatic_residues: array-like
        List of BioPython Residue objects included in the analysis
    model: object
        BioPython object containing a model of a PDB or MMCIF file

    Returns
    -------
    surface_exposed_res: array-like
        List of residue names corresponding to the surface exposed residues

    """
    surface = get_surface(model)
    surface_exposed_res = []
    cutoff = 3.03
    for residue in aromatic_residues:
        depth = residue_depth(residue, surface)
        if (depth <= cutoff):
            # custom residue
            if res_name_to_char.get(residue.get_resname()) == None:
                name = residue.get_resname()
                if "(" not in name:
                    name += "(" + residue.full_id[2] + ")"
            # standard residue
            else:
                name = res_name_to_char.get(residue.get_resname())
                name += str(residue.id[1])
                name += "(" + residue.full_id[2] + ")"
            surface_exposed_res.append(name)
    return surface_exposed_res


def calculate_asa(model, filename, AROM_LIST, chain_list):
    """Returns a list of surface exposed residues as determined by relative solvent accessibility.

    Only standard protein residues are currently supported. Non-protein and user specified custom residues will never be
    classified as surface exposed using this criteria.

    Parameters
    ---------
    aromatic_residues: array-like
        List of BioPython Residue objects included in the analysis
    model: object
        BioPython object containing a model of a PDB or MMCIF file
    AROM_LIST: array-like
        List of aromatic residue names included in the analysis

    See Also
    -------
    Bio.PDB.DSSP: module used to calculate relative solvent accessibility

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
    letter_list = []
    for res_name in AROM_LIST:
        if res_name_to_char.get(res_name):
            letter_list.append(res_name_to_char.get(res_name))
    try:
        dssp = DSSP(model, filename, acc_array="Wilke")
        keys = list(dssp.keys())
        for key in keys:
            if key[0] in chain_list and dssp[key][3] >= cutoff and dssp[key][1] in letter_list:
                goal_str = dssp[key][1] + str(key[1][1]) + "(" + str(key[0]) + ")"
                surface_exposed_res.append(goal_str)
    except Exception as e:
        pass
    return surface_exposed_res


def dist(x, y):
    """Returns Euclidean distance between two 3 dimensional points.

    Parameters
    ----------
    x: 3D numpy array
        Point 1
    y: 3D numpy array
        Point 2

    Returns
    -------
    float:
        Euclidean distance between two points

    """
    return np.sqrt(np.sum((x - y)**2))


def closestAtomDMatrix(residues, coef_alpha, exp_beta, r_offset):
    """Constructs distance matrix based on closest atom distance.

    Parameters
    ----------
    residues: array-like
        List of BioPython residues

    Returns
    -------
    node_label: dictionary
        List of node labels for graph
    distanceMatrix: array-like
        Distance matrix of residues

    """
    node_label = {}
    distanceMatrix = np.zeros((len(residues), len(residues)))
    pathwaysMatrix = np.zeros((len(residues), len(residues)))
    # iterate over all residues
    for i in range(0, len(residues)):
        res = residues[i]
        atm_list = list(res.get_atoms())
        resname = res.resname
        resnum = res.full_id[3][1]
        chain = res.full_id[2]
        # iterate over all residues that follow
        for j in range(i + 1, len(residues)):
            res2 = residues[j]
            atm_list_2 = list(res2.get_atoms())
            bestDist = sys.maxsize
            # iterate over all atoms of my residue
            for k in range(0, len(atm_list)):
                # iterate over all atoms of neighbor
                for l in range(0, len(atm_list_2)):
                    cur_atom = atm_list[k]
                    next_atom = atm_list_2[l]
                    x1 = cur_atom.coord[0]
                    y1 = cur_atom.coord[1]
                    z1 = cur_atom.coord[2]
                    x2 = next_atom.coord[0]
                    y2 = next_atom.coord[1]
                    z2 = next_atom.coord[2]
                    a = np.array((x1, y1, z1))
                    b = np.array((x2, y2, z2))
                    dist_a_b = dist(a, b)
                    if dist_a_b < bestDist:
                        bestDist = dist_a_b
            distanceMatrix[i][j] = bestDist
            distanceMatrix[j][i] = bestDist
            pathwaysMatrix[i][j] = pathways_model(bestDist, coef_alpha, exp_beta, r_offset)
            pathwaysMatrix[j][i] = pathwaysMatrix[i][j]
        # Node names in graph
        res_letter = res_name_to_char.get(resname)
        if res_letter:
            node_label[i] = res_letter + str(resnum) + "(" + chain + ")"
        else:
            node_label[i] = resname + "(" + chain + ")"
    return node_label, distanceMatrix, pathwaysMatrix


def comDMatrix(residues, coef_alpha, exp_beta, r_offset):
    """Constructs distance matrix based on distances between centers of mass.

    Parameters
    ----------
    residues: array-like
        List of BioPython residues

    Returns
    -------
    node_label: dictionary
        List of node labels for graph
    distanceMatrix: array-like
        Distance matrix of residues

    """
    try:
        node_label = {}
        com_d = []
        # Compute COMs (x0, y0, z0) of side-chains
        for i in range(len(residues)):
            res = residues[i]
            res.get_full_id()
            atm_list = list(res.get_atoms())
            resname = res.resname
            resnum = res.full_id[3][1]
            chain = res.full_id[2]
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
            # Node names in graph
            res_letter = res_name_to_char.get(resname)
            if res_letter:
                node_label[i] = res_letter + str(resnum) + "(" + chain + ")"
            else:
                node_label[i] = resname + "(" + chain + ")"
            com_d.append(np.array([com_x, com_y, com_z]))

        return node_label, distance_matrix(com_d, com_d), pathways_model_thrspace_penalty_COM(
            com_d, coef_alpha, exp_beta, r_offset)

    except Exception as e:
        raise Exception(e)


def process_standard_residues(standard_residue_list):
    """Generates customized BioPython residue objects with only side chain atoms included.

    Parameters
    ----------
    standard_residue_list: array-like
        List of BioPython residue objects corresponding to protein residues

    Returns
    -------
    res_list: array-like
        List of BioPython residue objects with only side chain atoms included.

    """
    res_list = []
    for i in range(0, len(standard_residue_list)):
        # make copy of residue
        standard_residue_list[i].get_full_id()
        res = standard_residue_list[i].copy()
        atm_ids = []
        atm_names = []
        for atm in res.get_atoms():
            atm_ids.append(atm.id)
            atm_names.append(atm.name)
        # if not side chain atom, remove from residue
        for k in range(0, len(atm_ids)):
            if res.resname == "TRP":
                if atm_names[k] not in TRP_sc:
                    res.detach_child(atm_ids[k])
            elif res.resname == "TYR":
                if atm_names[k] not in TYR_sc:
                    res.detach_child(atm_ids[k])
            elif res.resname == "PHE":
                if atm_names[k] not in PHE_sc:
                    res.detach_child(atm_ids[k])
            elif res.resname == "HIS":
                if atm_names[k] not in HIS_sc:
                    res.detach_child(atm_ids[k])
        res_list.append(res)
    return res_list


def get_user_res(serial_list, all_atoms, chain_selected, used_atoms, user_res_names):
    """Creates a customized BioPython Residue object corresponding to a user specified residue.

    Parameters
    ----------
    serial_list: array-like
        List of atom serial numbers included in residue
    all_atoms: array-like
        List of BioPython Atom objects
    chain_selected: array-like
        List of chains included in analysis
    used_atoms: array-like
        List of BioPython Atom objects already included in the analysis
    user_res_names: array-like
        List of user residue names already included in the analysis.

    Returns
    -------
    user_res: BioPython Residue object
        Customized Residue object corresponding to the atoms in serial_list

    """
    try:
        source_res = []
        atm_list = []
        for atm in all_atoms:
            if atm.serial_number in serial_list:
                if atm.serial_number in used_atoms:
                    message = "Error. Invalid atom serial number range. Atom " + str(
                        atm.serial_number) + " is in a specified standard or custom residue."
                    raise Exception(message)
                if atm.parent.parent.id not in chain_selected:
                    message = "Error. Invalid atom serial number range. Atom " + str(
                        atm.serial_number) + " is in chain " + str(atm.parent.parent.id) + "."
                    raise Exception(message)
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
        done = False
        k = 1
        name = "CUST"
        name += "-"
        while name + str(k) in user_res_names:
            k += 1
        user_res_names.append(name + str(k))
        user_res.resname = name + str(k)
        return user_res
    except Exception as e:
        message = []
        if "has no attribute" in str(e):
            message = "Error. One or more of the atoms specified do not exist. Please check your pdb file."
        if message:
            raise Exception(message)
        else:
            raise Exception(e)


def get_residues(all_residues, AROM_LIST, chain_list, eta_moieties,emap):
    """Returns list of standard aromatic residues.

    Runs through every residue in the structure, keeping only the TRP, TYR
    (optionally PHE and HIS if specified by user) residues on the chains
    that were chosen by the user.

    Parameters
    ----------
    all_residues: array-like
        List of every Residue in Structure
    AROM_LIST: array-like
        List of standard aromatic residues names to be included in graph
    chain_list: array-like
        List of chains chosen by user for analysis

    Returns
    -------
    residue_list: array-like
        List of standard aromatic Bio.PDB Residue objects to
        be included in graph.
    used_atoms: array-like
        List of atoms already included in selected standard or custom residues

    """
    residue_list = []
    for res in all_residues:
        if res.resname in AROM_LIST and res.parent.id in chain_list:
            res.get_full_id()
            arom_res = res.copy()
            residue_list.append(arom_res)
    residue_list = process_standard_residues(residue_list)
    return residue_list

def get_user_residues(custom_atm_string, all_atoms, chain_selected, used_atoms):
    """Generates customized BioPython Residue objects from a custom atom string.

    Users may not choose atoms that are already part of standard/included custom residues, nor those that are not
    part of the selected chains.

    Parameters
    ----------
    custom_atm_string: str
        Specified by user to select atoms for custom residues
    all_atoms: array-like
        List of all Atoms in structure
    chain_selected: array-like
        List of chains included in the analysis
    used_atoms: array-like
        List of atoms already included in the analyis

    Returns
    -------
    res_list: array-like
        List of customized Residue objects corresponding to the atoms in custom_atm_string

    Raises
    ------
    Exception:
    Invalid atom serial number range.

    Notes
    -----
    The custom atom string syntax is as follows:
    Custom residues are specified based on PDB atom serial number. Ranges are specified with a dash (ex. 10-20).
    Individual atoms may also be selected, separated by commas (ex. 1-10,15,17). To select multiple custom residues,
    the user encloses the atom ranges with parentheses, and then separate them with commas. In the graph and
    visualization, these residues are named CUS-1, CUS-2 etc. based on the order they were specified.

    """
    user_res_names = []
    try:
        custom_atm_string = custom_atm_string.strip()
        if not custom_atm_string == "-1":
            res_list = []
            l1 = custom_atm_string.split("),(")
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
                new_res = get_user_res(serial_number_list, all_atoms, chain_selected, used_atoms, user_res_names)
                res_list.append(new_res)
                for atm in new_res.get_atoms():
                    used_atoms.append(atm.serial_number)
            return res_list
        return []
    except Exception as e:
        if "Error" not in str(e):
            raise Exception("Error. Invalid atom serial number range. See the manual for proper syntax.")
        else:
            raise Exception(e)


def finish_graph(G, surface_exposed_res, chain_list, filename):
    """Draws and writes out the graph to file in downloadable form and for later use by the application, and
    returns the node labels.

    Parameters
    ----------
    G: NetworkX Graph
        Graph object constructed from distance matrix of residues
    surface_exposed_res: array-like
        List of surface exposed residues
    chain_list: array-like
        list of chains included in analysis
    filename: string
        File hash used for writing to file

    Returns
    -------
    A: pygraphviz agraph
        Graph object representing emap model

    See Also
    -----
    NetworkX: module used for graph theory representation
    PyGraphViz: program used to draw the graph

    Notes
    -----
    Disconnected vertices are not included in the final graph. For graphs with a single chain, the chain ID
    is omitted from the node label. The graph is drawn using PyGraphViz. The graph saved in the following formats:

    """
    for goal in surface_exposed_res:
        G.node[goal]['margin'] = '0.11'
        G.node[goal]['shape'] = 'box'
    # get rid of all disconnected nodes
    all_nodes = list(G.nodes())
    for node in all_nodes:
        if G[node] == {}:
            G.remove_node(node)


def process_user_options(include_Trp, include_Tyr, include_Phe, include_His, custom_residues, chain_selected,emap):
    """Processes the graph drawing options specified by the user.

    Removes PHE and HIS if not specified from AROM_LIST, and gathers the chains specified by the user.

    Parameters
    ----------
    include_Trp: str
        User specification. True if tryptophan included
    include_Tyr: str
        User specification. True if tyrosine included
    include_Phe: str
        User specification. True if phenylalanine included
    include_His: str
        User specification. True if histidine included
    chain_selected: array-like
        User specification. List of chains included in the analysis
    custom_residues: array-like
        List of strings corresponding to ETA moieties selected by user for analysis

    Returns
    -------
    AROM_LIST: array-like
        List of aromatic residue names included in the analysis
    chain_list: array-like
        List of chains included in the analysis

    """
    AROM_LIST = ['TRP', 'HIS', 'PHE', 'TYR']
    if include_Trp == "False":
        AROM_LIST.remove('TRP')
    if include_Tyr == "False":
        AROM_LIST.remove('TYR')
    if include_Phe == "False":
        AROM_LIST.remove('PHE')
    if include_His == "False":
        AROM_LIST.remove('HIS')
    chain_list = []
    for chain in chain_selected:
        chain = chain.strip()
        chain_list.append(chain)
    return AROM_LIST, chain_list


def create_graph(dmatrix, pathways_matrix, node_label, distanceCutoff, percentEdges, numStDevEdges):
    """Constructs the graph from the distance matrix and node labels.

    Parameters
    ----------
    dmatrix: array-like
        Distance matrix of aromatic residues
    node_label: array-like
        List of labels for residues in the graph
    distanceCutoff: float
        Distance cutoff for edges in the graph
    percentEdges: float
        Percentage of top (percentEdges)% edges per vertex included in the graph
    numStDevEdges: float
        For particular vertex, only edges with length < average length + numStDevEdges * SD included

    Returns
    -------
    G: NetworkX graph
        Graph of aromatic residues in protein

    Notes
    -----
    99th percentile of edges are kept with a ~20A filter on all edges.

    References
    ----------
    Gray, H. B.; Winkler, J. R. Long-Range Electron Transfer. Proc. Natl. Acad. Sci. U. S. A. 2005, 102 (10),
        3534 LP-3539.
        Reference for 20A filter on edges

    """
    np.set_printoptions(threshold=sys.maxsize)
    minval = np.min(dmatrix[dmatrix.nonzero()])
    G = nx.from_numpy_matrix(dmatrix)

    minval_pathways = np.min(pathways_matrix[pathways_matrix.nonzero()])
    G_pathways = nx.from_numpy_matrix(pathways_matrix)

    # keep 99th percentile and ~20A filter on all edges
    included_edges = []
    for node in G.nodes():
        edge_length_per_node = []
        weights = []
        for neighbor in G[node]:
            weights.append(G.get_edge_data(node, neighbor)['weight'])
        weights = sorted(weights)
        thresh_index = math.ceil(len(weights) * percentEdges / 100)
        for neighbor in G[node]:
            if weights.index(G.get_edge_data(node, neighbor)['weight']) <= thresh_index and \
                    G.get_edge_data(node, neighbor)['weight'] <= distanceCutoff:
                edge_length_per_node.append(G.get_edge_data(node, neighbor)['weight'])

        len_average, len_st_dev = 0.0, 0.0
        if edge_length_per_node != []:
            edge_length_per_node = np.array(edge_length_per_node, dtype='float64')
            len_average = np.average(edge_length_per_node)
            len_st_dev = np.std(edge_length_per_node)

        for neighbor in G[node]:
            if weights.index(G.get_edge_data(node, neighbor)['weight']) <= thresh_index and \
                    G.get_edge_data(node, neighbor)['weight'] <= distanceCutoff and \
                    np.round(G.get_edge_data(node, neighbor)['weight'], 8) <= np.round((len_average + numStDevEdges * len_st_dev), 8):
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

    G = G_pathways

    # scaling lengths
    for u, v, d in G.edges(data=True):
        d['len'] = d['weight'] / minval_pathways

    keys, values = list(node_label.keys()), list(node_label.values())
    for i in range(0, len(values)):
        val = values[i]
        if val.count("(") > 1:
            values[i] = val[:val.rindex("(")]
    node_label = dict(zip(keys, values))
    G = nx.relabel_nodes(G, node_label)

    for name_node in G.nodes():
        G.node[name_node]['style'] = 'filled'
        G.node[name_node]['fontname'] = 'Helvetica-Bold'
        G.node[name_node]['fontsize'] = 32
        G.node[name_node]['shape'] = "oval"
        G.node[name_node]['margin'] = '0.04'
        G.node[name_node]['fontcolor'] = "#000000"
        G.node[name_node]['color'] = '#708090'
        G.node[name_node]['penwidth'] = 2.0
        try:
            val = int(name_node[1])
            if 'Y' == name_node[0]:
                G.node[name_node]['fillcolor'] = '#96c8f0'
            elif 'W' == name_node[0]:
                G.node[name_node]['fillcolor'] = '#f07878'
            elif 'F' == name_node[0]:
                G.node[name_node]['fillcolor'] = '#f09664'
            elif 'H' == name_node[0]:
                G.node[name_node]['fillcolor'] = '#c8f0c8'
            else:
                G.node[name_node]['fillcolor'] = '#FFC0CB'
        except ValueError:
            G.node[name_node]['fillcolor'] = '#FFC0CB'
    for edge in G.edges():
        name_node1, name_node2 = edge[0], edge[1]
        G[name_node1][name_node2]['color'] = '#778899'
        G[name_node1][name_node2]['penwidth'] = 1.5
        G[name_node1][name_node2]['style'] = 'dashed'

    return G


def process(emap,
            chains,
            eta_moieties,
            distance_criteria=0,
            surface_exposed_bool=0,
            trp="True",
            tyr="True",
            phe="False",
            his="False",
            custom_atm_string="",
            distanceCutoff=20,
            percentEdges=1.0,
            numStDevEdges=1.0,
            coef_alpha=1.0,
            exp_beta=2.3,
            r_offset=0.0,
            graph_dest=""):
    """Constructs emap graph theory model based on user specs, and saves it to the emap object.

    Parameters
    ---------
    emap: emap object
        Object for storing state of emap analysis.
    chains: array-like
        List of strings corresponding to chains included in alaysis
    eta_moieties: array-like
        List of strings corresponding to residue names of BioPython custom residue objects
    distance_criteria: int, optional
        User specification. 0 for center of mass, 1 for closest atom
    surface_exposed_bool: int, optional
        User specification. 0 for residue depth, 1 for solvent accessibility
    trp: str, optional
        User specification. True if tryptophan included
    tyr: str, optional
        User specification. True if tyrosine included
    phe: str, optional
        User specification. True if phenylalanine included
    his: str, optional 
        User specification. True if histidine included
    custom_atm_string: str, optional
        User specification. Custom atom string specified by user
    distanceCutoff: float, optional
        User specification. Default is 20.0
    percentEdges: float, optional
        User specification. Default is 1.0
    numStDevEdges: float, optional
        User specification. Default is 1.0
    coef_alpha: float, optional
        User specification. Default is 1.0
    exp_beta: float, optional
        User specification. Default is 2.3
    r_offset: float, optional
        User specification. Default is 0.0
    store_graph: bool, optional
        Whether to store the graph in the emap object or not
    graph_dest: str, optional
        Destination to write the graph to disk instead of saving to emap object
    Raises
    ------
    Exception:
    Not enough residues to construct a graph

    """
    try:
        emap.reset_process()
        AROM_LIST, chain_list = process_user_options(trp, tyr, phe, his, eta_moieties, chains,emap)
        model = emap.structure[0]
        all_residues = list(model.get_residues())
        all_atoms = list(model.get_atoms())
        aromatic_residues = get_residues(all_residues, AROM_LIST, chain_list, eta_moieties,emap)
        for res in aromatic_residues:
            res_letter = res_name_to_char.get(res.resname)
            res_num = res.full_id[3][1]
            res_chain = res.full_id[2]
            emap.add_standard_residue(res, res_letter + str(res_num) + "(" + res_chain + ")")
        for resname in eta_moieties:
            AROM_LIST.append(resname)
            aromatic_residues.append(emap.get_residue(resname))
        used_atoms=[]
        for res in aromatic_residues:
            for atm in res.get_atoms():
                used_atoms.append(atm.serial_number)
        user_residues=[]
        if custom_atm_string:
            user_residues = get_user_residues(custom_atm_string, all_atoms, chains, used_atoms)    
        for res in user_residues:
            AROM_LIST.append(res.resname)
            emap.add_user_residue(res)
        if len(aromatic_residues) < 2:
            raise Exception(
                "Not enough residues to construct a graph. Please try different options or a different protein.")
        if int(distance_criteria) == 0:
            node_labels, dmatrix, pathways_matrix = comDMatrix(aromatic_residues, coef_alpha, exp_beta, r_offset)
        else:
            node_labels, dmatrix, pathways_matrix = closestAtomDMatrix(aromatic_residues, coef_alpha, exp_beta,
                                                                       r_offset)
        G = create_graph(dmatrix, pathways_matrix, node_labels, distanceCutoff, percentEdges, numStDevEdges)
        # define surface exposed residues
        if int(surface_exposed_bool) == 0:
            surface_exposed_res = calculateResidueDepth(aromatic_residues, model)
        else:
            pdb_file = emap.filename
            surface_exposed_res = calculate_asa(model, pdb_file, AROM_LIST, chain_list)
        finish_graph(G, surface_exposed_res, chain_list, emap.filename)
        emap.store_initial_graph(G)
        if graph_dest:
            emap.save_init_graph(dest=graph_dest+".svg",)
            emap.save_init_graph(dest=graph_dest+".png")
    except Exception as e:
        raise (Exception(e))
