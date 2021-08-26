from ..data import char_to_res_name
from networkx.algorithms import isomorphism

def get_edge_label(G, edge, edge_thresholds):
    dist = G.edges[edge]['distance']
    label = 2
    for thresh in edge_thresholds:
        if dist < thresh:
            break
        else:
            label += 1
    return label

# _exp denotes surface exposed residue
def get_numerical_node_label(u, pdb_id, res_to_num_label):
    if strip_res_number(u) in char_to_res_name:
        res_name = strip_res_number(u)
        result = res_to_num_label[res_name]
    elif (pdb_id + "_" + str(u)) in res_to_num_label:
        res_label = pdb_id + "_" + str(u)
        result = res_to_num_label[res_label]
    else:
        result = res_to_num_label["X"]
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
    return isomorphism.GraphMatcher(protein_graph, generic_subgraph, node_match=node_match, edge_match=edge_match)

