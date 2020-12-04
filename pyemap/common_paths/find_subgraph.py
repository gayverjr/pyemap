
import pyemap
import networkx as nx
from matplotlib import pyplot as plt

#buried labels, other is 7
res_labels = {
2:"W",
3:"Y",
4:"H",
5:"F",
6:"FAD",
7:"NP"
}

res_labels_inv = {
"W": 2,
"Y": 3,
"H": 4,
"F": 5,
"FAD": 6,
"FMN:":6,
"NP": 7
}

def get_numerical_label(G,u):
    res_name = strip_res_number(u)
    if res_name in res_labels_inv:
        return res_labels_inv[res_name]
    else:
        return 7

def strip_res_number(u):
    digit_idx=-1
    for i in range(0,len(u)):
        if u[i].isdigit():
            return u[:i]

# defines edge labels in the graphs for gspan, based on distance
def get_edge_label(G,edge):
    #dist = G.edges[edge[0], edge[1]]['distance']
    dist = G.edges[edge]['distance']
    if dist <= 10.0:
        return 2
    elif dist <=15.0:
        return 3
    else:
        return 4

def strip_res_number(u):
    digit_idx=-1
    for i in range(0,len(u)):
        if u[i].isdigit():
            return u[:i]


def generate_candidate_subgraph(graph):
    G = nx.Graph()
    has_FAD= False
    for i,node in enumerate(graph.nodes):
        G.add_node(i)
        num_label = get_numerical_label(graph,node)
        G.nodes[i]['label']= res_labels[num_label]
    for i,edge in enumerate(graph.edges):
        node1_idx = list(graph.nodes()).index(edge[0])
        node2_idx = list(graph.nodes()).index(edge[1])
        edge_label = graph.edges[edge]['label']
        G.add_edge(node1_idx,node2_idx,label=edge_label)
    return G

def get_candidates(G,source):
    # returns networkx node
    candidates = []
    for node in G.nodes:
        if strip_res_number(node)==source['label']:
            candidates.append(node)
    return candidates

def node_match(node1,node2):
    return node1['label'] == node2['label']

def edge_match(edge1,edge2):
    return edge1['label'] == edge2['label']

#non-recursive
def dfs_nr(full_graph,target_subgraph,src):
    found_subgraphs = []
    target_num_nodes = len(list(target_subgraph.nodes))
    target_num_edges = len(list(target_subgraph.edges))
    G = nx.Graph()
    G.add_node(src,shape=full_graph.nodes[src]['shape'])
    graph_stack = []
    stack = []
    stack.append(src)
    graph_stack.append(G)
    while(len(stack)):
        prev_node = stack.pop()
        G = graph_stack.pop()
        for cur_node in full_graph.neighbors(prev_node):
            if not G.has_edge(prev_node,cur_node):
                cur_G = G.copy()
                edge = (prev_node,cur_node)
                edge_label = get_edge_label(full_graph,edge)
                cur_G.add_node(cur_node)
                cur_G.add_edge(prev_node,cur_node,label=edge_label,distance=full_graph.edges[edge]['distance'])
                cur_num_nodes = len(list(cur_G.nodes))
                cur_num_edges = len(list(cur_G.edges))
                if target_num_nodes == cur_num_nodes and target_num_edges == cur_num_edges:
                    test_graph = generate_candidate_subgraph(cur_G)
                    if nx.is_isomorphic(test_graph,target_subgraph,edge_match=edge_match,node_match=node_match):
                        found_subgraphs.append(cur_G)
                elif target_num_nodes > cur_num_nodes:
                    stack.append(cur_node)
                    graph_stack.append(cur_G)
    return found_subgraphs

def find_sg(graph,subgraph,source):
    # graph is a networkx object generated by emap
    # subgraph is networkx object with node labels in the same format as our gspan stuff
    # source is a numerical label of the starting node
    # residue_map is a dict mapping node labels to residue types
    source_candidates = get_candidates(graph,source)
    sgs = []
    for src in source_candidates:
        sgs += dfs_nr(graph,subgraph,src)
    return sgs




