## Copyright (c) 2017-2022, James Gayvert, Ruslan Tazhigulov, Ksenia Bravaya
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
import math
from collections import OrderedDict
import warnings
from .data import res_name_to_char, side_chain_atoms, char_to_res_name
from .pyemap_exceptions import *
from .utils import validate_binary_params


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


def calculate_residue_depth(model, aromatic_residues, rd_cutoff):
    """Returns a list of surface exposed residues as determined by residue depth.

    Parameters
    ----------
    filename: str
        Name of pdb file to be analyzed
    aromatic_residues: list of :class:`Bio.PDB.Residue.Residue`
        residues included in the analysis
    rd_cutoff: float
        Cutoff for buried/surface exposed
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
    except Exception:
        warnings.warn("Unable to calculate residue depth. Check that MSMS is installed.", RuntimeWarning, stacklevel=2)
        return []


def calculate_rsa(filename, model, node_list, rsa_cutoff):
    """Returns a list of surface exposed residues as determined by relative solvent accessibility.

    Only standard protein residues are currently supported. Non-protein and user specified custom residues cannot be
    classified as surface exposed using this criteria.

    Parameters
    ---------
    filename: str
        Name of pdb file to be analyzed
    model: :class:`Bio.PDB.Model.Model`
        Model under analysis
    node_list : list of str
        List containing which standard residues are included in analysis
    rsa_cutoff: float
        Cutoff for buried/surface exposed

    References
    ---------
    Tien, M. Z.; Meyer, A. G.; Sydykova, D. K.; Spielman, S. J.; Wilke, C. O. PLoS ONE 2013, 8 (11).
        Reference for relative solvent accessibility cutoff of 0.05, and for MaxASA values
    """
    surface_exposed_res = []
    try:
        dssp = DSSP(model, filename, acc_array="Wilke")
        for key in dssp.keys():
            goal_str = dssp[key][1] + str(key[1][1]) + "(" + str(key[0]) + ")"
            if goal_str in node_list and dssp[key][3] >= rsa_cutoff:
                surface_exposed_res.append(goal_str)
    except Exception:
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
            num_atoms += 1
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
        valid = False
        for atm in res:
            sca = side_chain_atoms[res.resname]
            if atm.name in sca:
                valid = True
                break
        if valid:
            res_list.append(res)
        else:
            warnings.warn("The record for residue " + res.node_label +
                          " did not contain any of the following side chain atoms: " + str(sca) +
                          " and therefore is not included in the graph.",
                          RuntimeWarning,
                          stacklevel=2)
    return res_list


def create_user_res(serial_list, used_atoms, serial_dict, user_res_names):
    """Creates a customized BioPython Residue object corresponding to a user specified residue.

    Parameters
    ----------
    serial_list: list of int
        List of atom serial numbers included in residue
    used_atoms: list of :class:`Bio.PDB.Atom.Atom`
        Atoms already included in analysis
    serial_dict: dict of int, :class:`Bio.PDB.Atom.Atom`
        Dictionary of serial numbers and :class:`Bio.PDB.Atom.Atom` objects
    user_res_names: list of str
        User residue names already included in the analysis.

    Returns
    -------
    user_res: :class:`Bio.PDB.Residue.Residue`
        Customized residue object corresponding to the atoms in serial_list

    """
    source_res = serial_dict[serial_list[0]].parent
    source_res.get_full_id()
    user_res = source_res.copy()
    for atm in list(source_res.get_atoms()):
        user_res.detach_child(atm.id)
    for serial_number in serial_list:
        if serial_number in used_atoms:
            message = "Invalid atom serial number range. Atom " + str(
                serial_number) + " is already included in another residue."
            raise PyeMapUserResidueException(message)
        if serial_number not in serial_dict:
            message = str(serial_number) + " is not a valid serial number."
            raise PyeMapUserResidueException(message)
        user_res.add(serial_dict[serial_number])
    k = 1
    name = "CUST"
    name += "-"
    while name + str(k) in user_res_names:
        k += 1
    user_res_names.append(name + str(k))
    user_res.resname = name + str(k)
    user_res.node_label = user_res.resname
    return user_res


def get_standard_residues(all_residues, chain_list, res_names):
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
    res_names: list of str
        Included amino acids specified by 3 letter code

    Returns
    -------
    residue_list: list of :class:`Bio.PDB.Residue.Residue`
        Residues to be included in analysis

    """
    residue_list = []
    for res in all_residues:
        if res.resname in res_names and res.parent.id in chain_list:
            res.get_full_id()
            arom_res = res.copy()
            arom_res.get_full_id()
            residue_list.append(arom_res)
    residue_list = process_standard_residues(residue_list)
    return residue_list


def get_user_residues(custom, used_atoms, serial_dict):
    """Generates customized Bio.PDB.Residue.Residue objects from a custom atom string.

    Users may not choose atoms that are already part of standard residues or eta moieties, nor those that are not
    part of the selected chains.

    Parameters
    ----------
    custom: str
        Specified by user to select atoms for custom residues
    used_atoms: list of :class:`Bio.PDB.Atom.Atom`
        Atoms already included in analysis
    serial_dict: dict of int, :class:`Bio.PDB.Atom.Atom`
        Dictionary of serial numbers and :class:`Bio.PDB.Atom.Atom` objects

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
            if serial_number_list != []:
                serial_number_list = sorted(list(OrderedDict.fromkeys(serial_number_list)))
                new_res = create_user_res(serial_number_list, used_atoms, serial_dict, user_res_names)
            else:
                raise UserResidueException("Invalid atom serial number range. See the manual for proper syntax.")
            res_list.append(new_res)
            for atm in new_res:
                used_atoms.append(atm.serial_number)
        return res_list
    return []


def finish_graph(G, surface_exposed_res):
    """Sets surface exposed residues as boxes in graph and removes disconnected vertices.

    Parameters
    ----------
    G: :class:`networkx.graph`
        Graph object constructed from distance matrix of residues
    surface_exposed_res: list of str
        Names of surface exposed residues
    """
    for goal in surface_exposed_res:
        G.nodes[goal]['margin'] = '0.11'
        G.nodes[goal]['shape'] = 'box'
    # get rid of all disconnected nodes
    all_nodes = list(G.nodes())
    for node in all_nodes:
        if G[node] == {}:
            G.remove_node(node)


def filter_by_percent(G, percent_edges, num_st_dev_edges, distance_cutoff, coef_alpha, exp_beta, r_offset):
    included_edges = []
    minval = min(dict(G.edges).items(), key=lambda x: x[1]['weight'])[1]['weight']
    # should never happen, but just in case
    if minval == 0:
        minval = 1
    for _, _, d in G.edges(data=True):
        d['distance'] = d['weight']
        d['weight'] = pathways_model(d['weight'], coef_alpha, exp_beta, r_offset)
        d['len'] = d['weight'] / minval  # scaling factor for prettier graphs
    for node in G.nodes():
        edge_length_per_node = []
        weights = []
        for neighbor in G[node]:
            weights.append(G.edges[(node, neighbor)]['weight'])
        weights = sorted(weights)
        thresh_index = math.ceil(len(weights) * percent_edges / 100)
        for neighbor in G[node]:
            if weights.index(G.edges[(node, neighbor)]['weight']) <= thresh_index and \
                    G.edges[(node, neighbor)]['distance'] <= distance_cutoff:
                edge_length_per_node.append(G.edges[(node, neighbor)]['weight'])
        len_average, len_st_dev = 0.0, 0.0
        if edge_length_per_node != []:
            edge_length_per_node = np.array(edge_length_per_node, dtype='float64')
            len_average = np.average(edge_length_per_node)
            len_st_dev = np.std(edge_length_per_node)
        for neighbor in G[node]:
            if weights.index(G.edges[(node, neighbor)]['weight']) <= thresh_index and \
                    G.edges[(node, neighbor)]['distance'] <= distance_cutoff and \
                    np.round(G.edges[(node, neighbor)]['weight'], 8) <= np.round((len_average + num_st_dev_edges * len_st_dev), 8):
                included_edges.append([node, neighbor])
    excluded_edges = []
    for edge in G.edges():
        node1 = edge[0]
        node2 = edge[1]
        if ([node1, node2] not in included_edges) and ([node2, node1] not in included_edges):
            excluded_edges.append(edge)
    G.remove_edges_from(excluded_edges)


def filter_by_degree(G, max_degree, distance_cutoff, coef_alpha, exp_beta, r_offset):
    minval = min(dict(G.edges).items(), key=lambda x: x[1]['weight'])[1]['weight']
    # should never happen, but just in case
    if minval == 0:
        minval = 1
    remove_edges = []
    # impose hard cutoff, collect other edge weights
    for u, v, d in G.edges(data=True):
        if d['weight'] > distance_cutoff:
            remove_edges.append((u, v))
        else:
            d['distance'] = d['weight']
            d['weight'] = pathways_model(d['weight'], coef_alpha, exp_beta, r_offset)
            d['len'] = d['weight'] / minval  # scaling factor for prettier graphs
    G.remove_edges_from(remove_edges)
    remove_edges = []
    for node in G.nodes:
        if G.degree(node) > max_degree:
            for edge in sorted(list(G.edges(node)), key=lambda x: G.edges[x]['distance'])[max_degree:]:
                if edge not in remove_edges and edge[::-1] not in remove_edges:
                    remove_edges.append(edge)
    remove_edges = sorted(remove_edges, key=lambda x: G.edges[x]['distance'])
    for u, v in remove_edges:
        if G.degree(u) > max_degree or G.degree(v) > max_degree:
            G.remove_edge(u, v)


def create_graph(dmatrix, node_labels, edge_prune, coef_alpha, exp_beta, r_offset, distance_cutoff, percent_edges,
                 num_st_dev_edges, max_degree, eta_moieties):
    """Constructs the graph from the distance matrix and node labels.

    Parameters
    ----------
    dmatrix: numpy.array of float
        Distance matrix of aromatic residues.
    node_label: list of str
        Labels for residues in the graph.
    edge_prune: int
        0 for degree, 1 for percent
    residue_numbers:
        res numbers
    distance_cutoff,max_degree: float
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
    G = nx.relabel_nodes(G, node_labels)
    if edge_prune == 'DEGREE':
        filter_by_degree(G, max_degree, distance_cutoff, coef_alpha, exp_beta, r_offset)
    elif edge_prune == 'PERCENT':
        filter_by_percent(G, percent_edges, num_st_dev_edges, distance_cutoff, coef_alpha, exp_beta, r_offset)
    else:
        raise PyeMapGraphException("Invalid choice of edge_prune. Must be set to 'DEGREE' or 'PERCENT'.")
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


def store_params(emap, params):
    params.pop('chains')
    params.pop('eta_moieties')
    params.pop('emap')
    params.pop('include_residues')
    emap._process_params = params


def process(emap,
            chains=None,
            eta_moieties=None,
            dist_def='COM',
            sdef='RSA',
            edge_prune='PERCENT',
            include_residues=["Y", "W"],
            custom="",
            distance_cutoff=20,
            max_degree=4,
            percent_edges=1.0,
            num_st_dev_edges=1.0,
            rd_thresh=3.03,
            rsa_thresh=0.2,
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
    dist_def: str, optional
        Definition of distance matrix. 'COM' for center of mass, 'CATM' for closest atom
    sdef: str, optional
        Algorithm to use for surface exposure. 'RD' for residue depth, 'RSA' for relative solvent accessibility
    edge_prune: str, optional
        Algorithm for pruning edges. 'DEGREE' for degree, 'PERCENT' for percent
    include_residues: list of str
        Included amino acids specified by 1 letter code
    custom: str, optional
        Custom atom string specified by user
    distance_cutoff: float
         Defines a pure distance threshold. PyeMap will only keep edges with distances less than or equal distance_cutoff.
    max_degree: int, optional
        Maximum degree of any vertex. Only used when edge_prune is set to 'DEGREE'.
    percent_edges: float, optional
        Percent of edges to keep for each node. Only used when edge_prune is set to 'PERCENT'.
    num_st_dev_edges: float, optional
        Number of standard deviations of edges to keep. Only used when edge_prune is set to 'PERCENT'.
    rd_thresh: float, optional
        Threshold for buried/surface exposed for residue depth
    rsa_thresh: float, optional
        Threshold for buried/surface exposed for relative solvent accessbility
    coef_alpha,exp_beta,r_offset: float, optional
        Penalty function parameters.
    Raises
    ------
    RuntimeError:
        Not enough residues to construct a graph

    """
    max_degree = int(max_degree)
    dist_def, edge_prune, sdef = validate_binary_params(dist_def, edge_prune, sdef)
    emap_params = locals().copy()
    emap._reset_process()
    pdb_file = emap.file_path
    if chains is None:
        chains = [emap.chains[0]]
    if eta_moieties is None:
        eta_moieties = []
        for resname, moiety in emap.eta_moieties.items():
            if moiety.get_full_id()[2] in chains:
                eta_moieties.append(resname)
    else:
        for resname in eta_moieties:
            if resname not in emap.eta_moieties:
                raise PyeMapGraphException("Error: " + str(resname) + " is not a valid residue name.")
    res_names = []
    res_chars = []
    for i, res in enumerate(include_residues):
        if res.upper() in char_to_res_name:
            res_names.append(char_to_res_name[res.upper()])
            res_chars.append(res.upper())
        elif res.upper() in res_name_to_char:
            res_names.append(res.upper())
            res_chars.append(res_name_to_char[res.upper()])
        else:
            raise PyeMapGraphException("Error: " + str(res) + " is not a valid 1-letter or 3-letter amino acid code.")
    model = emap._structure[0]
    all_residues = get_standard_residues(model.get_residues(), chains, res_names)
    for resname in eta_moieties:
        all_residues.append(emap.eta_moieties[resname])
    used_atoms = []
    for res in all_residues:
        for atm in res:
            used_atoms.append(atm.serial_number)
    user_residues = []
    if custom:
        selection = Bio.PDB.Selection.unfold_entities(emap._structure, target_level='A')
        serial_numbers = [atom.serial_number for atom in selection]
        serial_dict = dict(zip(serial_numbers, selection))
        user_residues = get_user_residues(custom, used_atoms, serial_dict)
    all_residues += user_residues
    if len(all_residues) < 2:
        raise PyeMapGraphException("Not enough residues to construct a graph.")
    node_labels = {}
    for i in range(0, len(all_residues)):
        node_labels[i] = all_residues[i].node_label
    if dist_def == 'COM':
        dmatrix = com_dmatrix(all_residues)
    elif dist_def == 'CATM':
        dmatrix = closest_atom_dmatrix(all_residues)
    else:
        raise PyeMapGraphException(
            "Invalid choice of dist_def. Must be set to 'COM' (center of mass) or 'CATM'(closest atom).")
    G = create_graph(dmatrix, node_labels, edge_prune, coef_alpha, exp_beta, r_offset, distance_cutoff, percent_edges,
                     num_st_dev_edges, max_degree, emap.eta_moieties.keys())
    G.graph['pdb_id'] = emap.pdb_id
    if len(G.edges()) == 0:
        raise PyeMapGraphException("Not enough edges to construct a graph.")
    # define surface exposed residues
    if sdef is None:
        warnings.warn("Protein surface will not be computed. All residues will be classified as buried...")
        surface_exposed_res = []
    else:
        try:
            if sdef == 'RD':
                surface_exposed_res = calculate_residue_depth(model, all_residues, rd_thresh)
            elif sdef == 'RSA':
                surface_exposed_res = calculate_rsa(pdb_file, model, node_labels.values(), rsa_thresh)
            else:
                surface_exposed_res = []
                warnings.warn("Invalid choice of surface definition. sdef must be set to 'RD' or 'RSA'.")
        except Exception:
            warnings.warn("Computing protein surface failed. All residues will be classified as buried...")
            surface_exposed_res = []
    finish_graph(G, surface_exposed_res)
    for res in all_residues:
        emap._add_residue(res)
    for res in user_residues:
        emap.user_residues[res.resname] = res
    emap._store_initial_graph(G)
    emap._include_residues = res_chars
    store_params(emap, emap_params)
    return emap
