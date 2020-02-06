# PyeMap: A python package for automatic identification of electron and hole transfer pathways in proteins.
# Copyright(C) 2017-2020 Ruslan Tazhigulov, James Gayvert, Ksenia Bravaya (Boston University, USA)
import numpy as np
from .shortest_paths import yens_shortest_paths, dijkstras_shortest_paths

def _finish_graph(G, original_shape_start, source):
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
    G.nodes[source]['fillcolor'] = '#FFD700FF'
    G.nodes[source]['penwidth'] = 6.0
    G.nodes[source]['shape'] = original_shape_start
    # if not is not involved in pathways, make it less opaque on the graph.
    # change font color to slate gray, and change transparency of edges
    for name_node in G.nodes():
        if len(G.nodes[name_node]['fillcolor']) != 9:
            G.nodes[name_node]['fillcolor'] += '40'
            G.nodes[name_node]['fontcolor'] = '#708090'
    for edge in G.edges():
        name_node1, name_node2 = edge[0], edge[1]
        if G[name_node1][name_node2]['style'] == 'dashed':
            G[name_node1][name_node2]['color'] = '#7788994F'
    # draw graph


def find_paths(emap, source, target=None, max_paths=10):
    """Function which calculates pathways from source to target or surface exposed residues.

    Performs shortest path analysis on source and (optionally) target residues.

    Parameters
    ---------
    emap: :class:`~pyemap.emap` 
        Object for storing state of emap analysis.
    source: str
        source node for analysis
    target: str, optional
        target node for analysis
    max_paths: int, optional
        maximum number of paths to search for in yen's algorithm
    """
    # read in graph from file
    emap._reset_paths()
    G = emap.init_graph.copy()
    for u, v, d in G.edges(data=True):
        d['weight'] = np.float64(d['weight'])
    # process source and target
    source = source.strip()
    original_shape_start = G.nodes[source]['shape']
    G.nodes[source]['shape'] = 'oval'
    if target:
        target = target.strip()
        branches = yens_shortest_paths(G, source, target, max_paths=max_paths)
        # color target node blue
        G.nodes[target]['fillcolor'] = '#40e0d0FF'
        G.nodes[target]['penwidth'] = 6.0
    else:
        surface_exposed = emap.get_surface_exposed_residues()
        if source in surface_exposed:
            surface_exposed.remove(source)
        branches = dijkstras_shortest_paths(G, source, surface_exposed)
    _finish_graph(G, original_shape_start, source)
    emap._store_paths_graph(G)
    if target:
        emap._store_paths(branches, yens=True)
    else:
        emap._store_paths(branches)
