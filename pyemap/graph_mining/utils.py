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

from ..data import char_to_res_name
from pysmiles import write_smiles, read_smiles
import networkx as nx
import math


def extract_chain(resname):
    try:
        return resname[resname.index('(') + 1:resname.index(")")]
    except Exception:
        return ''


def get_edge_label(G, edge, edge_thresholds):
    dist = G.edges[edge]['distance']
    try:
        label = 1
        for thresh in edge_thresholds:
            if dist < thresh:
                break
            else:
                label += 1
        return label
    except Exception:
        return 1


def write_graph_smiles(generic_subgraph):
    ''' Returns pseudo-SMILES string for supplied graph.

    Parameters
    -----------
    generic_subgraph: :class:`networkx.Graph`
        Graph to be transformed into string representation. The :attr:`label` attribute 
        of each node should be set to the 1-letter amino acid code or special character.

    Returns
    --------
    pseudosmiles: str
        String representation of graph of interest
    '''
    G = generic_subgraph.copy()
    element_dict = {}
    num_chars = 0
    for i, node in enumerate(G.nodes):
        G.nodes[node]['element'] = 'C' + str(i)
        element_dict['C' + str(i)] = G.nodes[node]['label']
        num_chars += len(G.nodes[node]['label'])
    proper_smiles = write_smiles(G)
    for key, val in element_dict.items():
        proper_smiles = proper_smiles.replace(key, val)
    # linear case
    if len(proper_smiles.replace('[', '').replace(']', '')) == num_chars:
        return proper_smiles.replace('[', '').replace(']', '')
    else:
        return proper_smiles


def get_numerical_node_label(u, res_to_num_label):
    if strip_res_number(u) in char_to_res_name and strip_res_number(u) in res_to_num_label:
        res_name = strip_res_number(u)
        result = res_to_num_label[res_name]
    elif strip_res_number(u) in char_to_res_name:
        result = res_to_num_label['X']
    else:
        result = res_to_num_label["#"]
    return result


def strip_res_number(u):
    for i in range(0, len(u)):
        if u[i].isdigit():
            return u[:i]


def node_match(node1, node2):
    return node1['num_label'] == node2['num_label']


def edge_match(edge1, edge2):
    return edge1['num_label'] == edge2['num_label']


def get_graph_matcher(protein_graph, generic_subgraph):
    return nx.algorithms.isomorphism.GraphMatcher(protein_graph,
                                                  generic_subgraph,
                                                  node_match=node_match,
                                                  edge_match=edge_match)


def set_defaults(kwargs):
    default = {
        'distance_cutoff': 20,
        'max_degree': 4,
        'dist_def': 'COM',
        'sdef': 'RSA',
        'edge_prune': 'PERCENT',
        'percent_edges': 1.0,
        'num_st_dev_edges': 1.0,
        'rd_thresh': 3.03,
        'rsa_thresh': 0.2,
        'coef_alpha': 1.0,
        'exp_beta': 2.3,
        'r_offset': 0.0
    }
    for arg in default:
        if arg not in kwargs:
            kwargs[arg] = default[arg]
    return kwargs


def make_pretty_subgraph(sg):
    for name_node in sg.nodes():
        sg.nodes[name_node]['style'] = 'filled'
        sg.nodes[name_node]['fontname'] = 'Helvetica-Bold'
        sg.nodes[name_node]['fontsize'] = 14
        sg.nodes[name_node]['margin'] = '0.04'
        sg.nodes[name_node]['fontcolor'] = "#000000"
        sg.nodes[name_node]['color'] = '#708090'
        sg.nodes[name_node]['penwidth'] = 2.0
        if (len(sg.nodes[name_node]['label']) == 1) or (len(sg.nodes[name_node]['label']) > 1
                                                        and sg.nodes[name_node]['label'][1].isdigit()):
            if 'Y' == sg.nodes[name_node]['label'][0]:
                sg.nodes[name_node]['fillcolor'] = '#96c8f0'
            elif 'W' == sg.nodes[name_node]['label'][0]:
                sg.nodes[name_node]['fillcolor'] = '#f07878'
            elif 'F' == sg.nodes[name_node]['label'][0]:
                sg.nodes[name_node]['fillcolor'] = '#f09664'
            elif 'H' == sg.nodes[name_node]['label'][0]:
                sg.nodes[name_node]['fillcolor'] = '#c8f0c8'
            else:
                sg.nodes[name_node]['fillcolor'] = '#FFC0CB'
        else:
            sg.nodes[name_node]['fillcolor'] = '#FFC0CB'
    for edge in sg.edges:
        try:
            dist = '{0:.2f}'.format(sg.edges[edge]['distance'])
            sg.edges[edge]['len'] = 1.0 + math.log10(float(dist))
            sg.edges[edge]['label'] = dist
        except Exception:
            pass
        sg.edges[edge]['fontname'] = 'Helvetica'
        sg.edges[edge]['color'] = '#778899'
        sg.edges[edge]['penwidth'] = 1.5
        sg.edges[edge]['style'] = 'dashed'
    return sg


def nodes_and_edges_from_smiles(smiles_str, edge_thresholds=[], residue_categories=[]):
    ''' Returns all possible combinations of nodes and edges based on graph string and edge thresholds and residue categories.

    Parameters
    ----------
    graph_str: str
        Specification of graph
    edge_thresholds: list of float
        Edge thresholds
    residue_categories: list of str
        List of 1 letter amino acid codes
    '''
    if '[' not in smiles_str:
        new_smiles = ""
        for char in smiles_str:
            new_smiles += '[{}]'.format(char)
        smiles_str = new_smiles
    # replace some problematic characters
    if 'H' in smiles_str:
        smiles_str = smiles_str.replace('H', 'He')
    if '#' in smiles_str:
        smiles_str = smiles_str.replace('#', 'Np')
    base_graph = read_smiles(smiles_str)
    node_list = []
    for node in base_graph.nodes:
        try:
            if base_graph.nodes[node]['element'] == 'He':
                node_list.append('H')
            elif base_graph.nodes[node]['element'] == 'Np':
                node_list.append('#')
            else:
                node_list.append(base_graph.nodes[node]['element'])
        except Exception:
            node_list.append("*")
    l1 = []
    for i in range(0, len(edge_thresholds) + 1):
        l1.append(i + 1)
    from itertools import product
    edge_combs = list(product(l1, repeat=len(base_graph.edges)))
    edges = list(base_graph.edges)
    indices = [i for i, x in enumerate(node_list) if x == "*"]
    if len(indices) == 0:
        node_combs = [node_list]
    else:
        node_combs = []
        wildcard_combs = list(product(residue_categories, repeat=len(indices)))
        for comb in wildcard_combs:
            for i, idx in enumerate(indices):
                node_list[idx] = comb[i]
            node_combs.append(node_list.copy())
    return node_combs, edge_combs, edges
