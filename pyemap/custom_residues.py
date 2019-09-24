# This Python script is a part of
# eMAP: online mapping of electron transfer channels in biomolecules
# Copyright(C) 2017-2018 Ruslan Tazhigulov, James Gayvert, Melissa Wei, Ksenia Bravaya (Boston University, USA)
"""Takes a list of non-standard Bio.PDB Residue objects, and returns a list of custom residues.

Module used to identify and draw custom aromatic/conjugated residues. Aromatic residues are identified using
experimental bond lengths, and customized BioPython Residue objects are constructed for those aromatic
non protein residues. The residues are drawn using the RDKit Module and written out to file for visualization
on the front end.

"""
from collections import defaultdict
from heapq import heappop, heappush

import networkx as nx
import numpy as np
from rdkit import Chem
from .data import *


def is_pi_bonded(cur_atom, next_atom):
    """Determines whether two atoms are pi bonded based on experimental bond lengths.

    Parameters
    ----------
    cur_atom: BioPython Atom object
    next_atom: BioPython Atom object
        The two atoms being considered

    Returns
    -------
    True: if within 3 standard deviations of equilibrium pi bond length
    False: Otherwise

    References
    ----------
    Haynes, W. M. (Ed.). (2014). CRC Handbook of Chemistry and Physics. CRC Press.
        Reference for equilibrium pi bond lengths and standard deviations

    """
    x1 = cur_atom.coord[0]
    y1 = cur_atom.coord[1]
    z1 = cur_atom.coord[2]
    x2 = next_atom.coord[0]
    y2 = next_atom.coord[1]
    z2 = next_atom.coord[2]
    v1 = np.array((x1, y1, z1))
    v2 = np.array((x2, y2, z2))
    bond = str(cur_atom.element.upper()) + str(next_atom.element.upper())
    cutoff = SB_means.get(bond)
    if cutoff:
        cutoff -= 2 * SB_std_dev.get(bond)
        return dist(v1, v2) <= cutoff
    else:
        return False


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

def find_conjugated_systems(atoms, res_names):
    """Finds conjugated systems within a BioPython residue object, and returns them as individual customized BioPython
    Residue objects.

    Parameters
    ---------
    atoms: array-like
        List of atoms in the residue
    res_names: arary-like
        List of already used names for custom residues

    Returns
    -------
    custom_res_list: array-like
        List of customized BioPython residue objects

    Notes
    -----
    A chemical graph is constructed from the atoms in this residue using O, P, N, C, and S atoms. In this graph, the
    only edges that are included are those between atoms that considered to be pi-bonded (see is_pi_bonded). The
    chemical graph is now a forest, and so the next step is collect each connected component subgraph. Only subgraphs
    that contain cycles (aromatic) or are larger than 10 are considered (extended conjugated systems). For each subgraph,
    the atoms are collected and a customized Residue object is constructured and named.
    """
    # first let's create the chemical graph structure
    # only these elements will be considered in our search
    arom_atoms = ['O', 'P', 'N', 'C', 'S']
    cust_graph = nx.Graph()
    for i in range(len(atoms)):
        for k in range(i, len(atoms)):
            if (not i == k) and is_pi_bonded(atoms[i], atoms[k]):
                if atoms[i].element in arom_atoms and atoms[k].element in arom_atoms:
                    cust_graph.add_edge(i, k)
                    cust_graph.nodes[i]["element"] = atoms[i].element
                    cust_graph.nodes[i]["coords"] = atoms[i].coord
                    cust_graph.nodes[k]["element"] = atoms[k].element
                    cust_graph.nodes[k]["coords"] = atoms[k].coord

    # now we have a forest, let's get each individual tree (only take cycles or bigger than 10)
    subgraphs = []
    all_subgraphs = (cust_graph.subgraph(c)
                     for c in nx.connected_components(cust_graph))
    for graph in all_subgraphs:
        if nx.cycle_basis(graph) or len(graph.nodes()) >= 10:
            subgraphs.append(graph)
    custom_res_list = []
    count = 1

    # iterate over all of the subgraphs
    for graph in subgraphs:
        res_name = str(atoms[0].parent.get_resname()) + str(
            atoms[0].get_full_id()[3][1]) + "(" + str(
                atoms[0].get_full_id()[2]) + ")"
        if len(subgraphs) > 1:
            res_name += '-' + str(count)
            while res_name in res_names:
                count += 1
                res_name = res_name.split('-')[0] + '-' + str(count)
        atm_list = []
        for node in graph.nodes():
            atm_list.append(atoms[node])
        custom_res_list.append(create_custom_residue(atm_list, res_name))
        res_names.append(res_name)

    return custom_res_list


def create_custom_residue(atm_list, res_name):
    """Generates customized BioPython residue object corresponding to a conjugated system.

    Parameters
    ----------
    atm_list: array-like
        List of atoms to be included
    res_name: str
        Name of custom residue

    Returns
    -------
    custom_residue: BioPython Residue Object
        customized BioPython residue object corresponding to a conjugated system

    """
    # copy residue, assign resname
    atm_list[0].parent.get_full_id()  # initialize instance variables
    custom_res = atm_list[0].get_parent().copy()
    custom_res.resname = res_name
    custom_res.node_label = res_name

    # get serial numbers of aromatic atoms
    atm_serial_list = []
    for atm in atm_list:
        atm_serial_list.append(atm.serial_number)

    # remove atoms not in aromatic ring by serial number
    atm_list = list(custom_res.get_atoms())
    for k in range(0, len(atm_list)):
        if atm_list[k].serial_number not in atm_serial_list:
            custom_res.detach_child(atm_list[k].id)

    return custom_res


def process_custom_residues(non_standard_residue_list):
    """Main method of custom residues module.

    Executes all major functions of this module.

    Parameters
    ---------
    non_standard_residue_list: array-like
        List of non-protein residues in the structure

    Returns
    ------
    custom_res: array-like
        List of customized BioPython residue objects corresponding to aromatic/conjugated systems that are not part of
        standard protein residues

    """
    res_names = []
    custom_res = []

    for residue in non_standard_residue_list:
        atm_list = list(residue.get_atoms())
        conj_systems = find_conjugated_systems(
            atm_list, res_names)
        if conj_systems:
            for system in conj_systems:
                res_names.append(system.resname)
            custom_res += conj_systems

    # Clusters
    for residue in non_standard_residue_list:
        if residue.resname in clusters:
            atm_list = list(residue.get_atoms())
            res_name = str(atm_list[0].parent.get_resname()) + str(
                atm_list[0].get_full_id()[3][1]) + "(" + str(
                    atm_list[0].get_full_id()[2]) + ")"
            custom_res.append(create_custom_residue(atm_list, res_name))
    return custom_res
