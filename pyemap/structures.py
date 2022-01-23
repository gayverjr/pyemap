# PyeMap: A python package for automatic identification of electron and hole transfer pathways in proteins.
# Copyright(C) 2017-2020 Ruslan Tazhigulov, James Gayvert, Ksenia Bravaya (Boston University, USA)
import networkx as nx
import numpy as np
from .data import SB_means,SB_std_dev

def is_part_of_cycle(node, res_graph):
    for cycle in nx.cycle_basis(res_graph):
        if node in cycle:
            return True
    return False


def is_close(node, node2, res_graph):
    bond = str(res_graph.nodes[node]["element"].upper()) + str(res_graph.nodes[node2]["element"].upper())
    dist = np.sqrt(np.sum((np.array(res_graph.nodes[node]["coords"]) -\
     						np.array(res_graph.nodes[node2]["coords"]))**2))
    cutoff = SB_means.get(bond) + 3 * SB_std_dev.get(bond)
    return dist < cutoff


def cleanup_bonding(res_graph):
    '''Connects nodes that should be connected to fix broken aromaticity.
    
    Parameters
    -----------
    res_graph: :class:`networkx.Graph`
        residue graph
    '''
    for node in res_graph.nodes:
        if not is_part_of_cycle(node, res_graph) and len(list(res_graph.neighbors(node))) < 3:
            closest_neighbor = []
            min_dist = 10000
            for node2 in res_graph.nodes:
                if node2 != node and node2 not in res_graph.neighbors(node) and len(list(
                        res_graph.neighbors(node2))) < 3:
                    if is_close(node, node2, res_graph):
                        dist = np.sqrt(np.sum((np.array(res_graph.nodes[node]["coords"]) -\
                         					   np.array(res_graph.nodes[node2]["coords"]))**2))
                        if dist < min_dist:
                            min_dist = dist
                            closest_neighbor = node2
            if closest_neighbor:
                res_graph.add_edge(node, closest_neighbor)


def remove_atoms(prev, cur, remove_list, res_graph):
    if not is_part_of_cycle(cur, res_graph):
        for neighbor in res_graph.neighbors(cur):
            if neighbor != prev:
                remove_atoms(cur, neighbor, remove_list, res_graph)
        remove_list.append(cur)


def remove_side_chains(res_graph):
    ''' Removes non-aromatic sides chains on aromatic eta moieties.

    Parameters
    -----------
    res_graph: :class:`networkx.Graph`
        residue graph
    '''
    remove_list = []
    for node in res_graph.nodes:
        if not is_part_of_cycle(node,
                                res_graph) and res_graph.nodes[node]["element"] == "C" and node not in remove_list:
            remove_atoms(-1, node, remove_list, res_graph)
    for node in remove_list:
        res_graph.remove_node(node)
