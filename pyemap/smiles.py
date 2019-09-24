from heapq import heappop, heappush
import networkx as nx
import numpy as np
from .data import *
from collections import defaultdict


def buildSmiles(graph, cur, prev):
    """Two pass depth first search algorithm on chemical graph to generate smiles string.

    Parameters
    ----------
    graph: NetworkX chemical graph
    atoms: array-like
        List of atoms in structure
    cur: BioPython Atom object
    prev: BioPython Atom object
        Atoms being considered at this step of the iteration

    References
    ----------
    Varnek, A. Tutorials in Chemoinformatics; John Wiley & Sons, Inc.: Hoboken, NJ, 2017.
    """
    visited.add(cur)
    seq = ''
    seq += graph.nodes[cur]["element"].lower()
    for d in closingClosures[cur]:
        seq += d
        heappush(digits, d[-1])
    for a in openingClosures[cur]:
        d = str(heappop(digits))
        seq += d
        closingClosures[a].append(d)
    branches = []
    neighbors = list(graph.neighbors(cur))
    if prev in neighbors:
        neighbors.remove(prev)
    for neighbor in neighbors:
        if neighbor not in visited:
            branches.append(buildSmiles(graph, neighbor, cur))
    for branch in branches[:-1]:
        seq += "(" + branch + ")"
    if len(branches) > 0:
        seq += branches[-1]
    return seq


def getClosures(graph, cur, prev):
    """Two pass depth first search algorithm on chemical graph to generate smiles string.

    graph: NetworkX chemical graph
    cur: BioPython Atom object
    prev: BioPython Atom object
        Atoms being considered at this step of the iteration

    References
    ----------
    Varnek, A. Tutorials in Chemoinformatics; John Wiley & Sons, Inc.: Hoboken, NJ, 2017.

    """
    ancestor.add(cur)
    visited.add(cur)
    neighbors = list(graph.neighbors(cur))
    if prev in neighbors:
        neighbors.remove(prev)
    for neighbor in neighbors:
        if neighbor in ancestor:
            openingClosures[neighbor].append(cur)
        elif neighbor not in visited:
            getClosures(graph, neighbor, cur)
    ancestor.remove(cur)


def getSimpleSmiles(graph):
    """Two pass depth first search algorithm on chemical graph to generate smiles string.

    Parameters
    ----------
    graph: NetworkX chemical graph
    atoms: array-like
        List of atoms in structure

    References
    ----------
    Varnek, A. Tutorials in Chemoinformatics; John Wiley & Sons, Inc.: Hoboken, NJ, 2017.

    """
    root = list(graph.nodes())[0]
    global visited, ancestor, openingClosures, closingClosures, digits
    visited = set()
    ancestor = set()
    openingClosures = defaultdict(list)
    getClosures(graph, root, None)
    closingClosures = defaultdict(list)
    digits = [str(x) for x in range(1, 10)]
    visited = set()
    return buildSmiles(graph, root, None)

def is_part_of_cycle(node,res_graph):
    for cycle in nx.cycle_basis(res_graph):
        if node in cycle:
            return True
    return False

def is_close(node,node2,res_graph):
    v1= np.array(res_graph.nodes[node]["coords"])
    v2= np.array(res_graph.nodes[node2]["coords"])
    bond = str(res_graph.nodes[node]["element"].upper()) + str(res_graph.nodes[node2]["element"].upper())
    dist = np.sqrt(np.sum((v1 - v2)**2))
    cutoff = SB_means.get(bond) + 3 * SB_std_dev.get(bond)
    return dist < cutoff

def cleanup_bonding(res_graph):
    for node in res_graph.nodes:
        if not is_part_of_cycle(node,res_graph) and len(list(res_graph.neighbors(node)))<3:
            closest_neighbor = []
            min_dist = 10000
            for node2 in res_graph.nodes:
                if node2 != node and node2 not in res_graph.neighbors(node) and len(list(res_graph.neighbors(node2)))<3:
                    if is_close(node,node2,res_graph):
                        v1= np.array(res_graph.nodes[node]["coords"])
                        v2= np.array(res_graph.nodes[node2]["coords"])
                        dist = np.sqrt(np.sum((v1 - v2)**2))
                        if dist < min_dist:
                            min_dist = dist
                            closest_neighbor = node2
            if closest_neighbor:
                v1= np.array(res_graph.nodes[node]["coords"])
                v2= np.array(res_graph.nodes[closest_neighbor]["coords"])
                res_graph.add_edge(node,closest_neighbor)
        
def remove_atoms(prev,cur,remove_list,res_graph):
    if not is_part_of_cycle(cur,res_graph):
        for neighbor in res_graph.neighbors(cur):
            if neighbor != prev:
                remove_atoms(cur,neighbor,remove_list,res_graph)
        remove_list.append(cur)
            
def remove_side_chains(res_graph):
    remove_list =[]
    for node in res_graph.nodes:
        if not is_part_of_cycle(node,res_graph) and res_graph.nodes[node]["element"]=="C" and node not in remove_list:
            remove_atoms(-1,node,remove_list,res_graph)
    for node in remove_list:
        res_graph.remove_node(node)