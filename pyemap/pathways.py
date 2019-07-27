import itertools
import os
import string
import sys
import networkx as nx
import numpy as np
import pygraphviz as pg
from networkx.drawing.nx_agraph import from_agraph, to_agraph
from .dijkstras import yens_shortest_paths, dijkstras_shortest_paths


def draw_graph(G, original_shape_start, source):
    """Draws the graph with the shortest pathways highlighted.

    Parameters
    ----------
    G: NetworkX graph object
        Residue graph
    source: str
        name of source node
    original_shape_start: str
        shape of source node

    """
    G.node[source]['fillcolor'] = '#FFD700FF'
    G.node[source]['penwidth'] = 6.0
    G.node[source]['shape'] = original_shape_start
    # if not is not involved in pathways, make it less opaque on the graph.
    # change font color to slate gray, and change transparency of edges
    for name_node in G.nodes():
        if len(G.node[name_node]['fillcolor']) != 9:
            G.node[name_node]['fillcolor'] += '40'
            G.node[name_node]['fontcolor'] = '#708090'
    for edge in G.edges():
        name_node1, name_node2 = edge[0], edge[1]
        if G[name_node1][name_node2]['style'] == 'dashed':
            G[name_node1][name_node2]['color'] = '#7788994F'
    # draw graph
    A_new = to_agraph(G)
    A_new.graph_attr.update(ratio=1.0, overlap="ipsep", mode="ipsep", splines="true")
    A_new.layout(args="-n2")
    return A_new


def find_pathways(emap, source, target=None):
    """Function which calculates pathways from source to target or surface exposed residues.

    Performs shortest path analysis on source and (optionally) target residues. After analysis is completed, the pathways
    graph is drawn and saved to the passed emap object.

    Parameters
    ---------
    emap: emap object
        Object for storing state of emap analysis.
    source: str
        source node for analysis
    target: str, optional
        target node for analysis
    """
    # read in graph from file
    A = emap.init_agraph
    G = from_agraph(A)
    for u, v, d in G.edges(data=True):
        d['weight'] = np.float64(d['weight'])
    # process source and target
    source = source.strip()
    original_shape_start = G.node[source]['shape']
    G.node[source]['shape'] = 'oval'
    if target:
        target = target.strip()
        shortest_paths = yens_shortest_paths(G, source, target)
        # color target node blue
        G.node[target]['fillcolor'] = '#40e0d0FF'
        G.node[target]['penwidth'] = 6.0
    else:
        goals = []
        for n, d in G.nodes(data=True):
            if d['shape'] == "box":
                goals.append(n)
        shortest_paths = dijkstras_shortest_paths(G, source, goals)
    A = draw_graph(G, original_shape_start, source)
    emap.store_paths(shortest_paths)
    emap.store_paths_agraph(A)
