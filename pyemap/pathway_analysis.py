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
