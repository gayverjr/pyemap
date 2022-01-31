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
