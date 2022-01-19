from ..data import char_to_res_name
from networkx.algorithms import isomorphism
import math

def extract_chain(resname):
    try:
        return resname[resname.index('(')+1:resname.index(")")]
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

#adapted from: https://github.com/pckroon/pysmiles/
def write_graph_smiles(generic_subgraph):
    import networkx as nx
    from collections import defaultdict
    G = generic_subgraph.copy()
    start = min(G.nodes, key=lambda x: G.degree(x))
    dfs_successors = nx.dfs_successors(G, source=start)
    predecessors = defaultdict(list)
    for node_key, successors in dfs_successors.items():
        for successor in successors:
            predecessors[successor].append(node_key)
    predecessors = dict(predecessors)
    # We need to figure out which edges we won't cross when doing the dfs.
    # These are the edges we'll need to add to the smiles using ring markers.
    edges = set()
    for n_idx, n_jdxs in dfs_successors.items():
        for n_jdx in n_jdxs:
            edges.add(frozenset((n_idx, n_jdx)))
    #total_edges = set(map(frozenset, G.edges))
    branch_depth = 0
    branches = set()
    to_visit = [start]
    smiles = ''
    while to_visit:
        current = to_visit.pop()
        if current in branches:
            branch_depth += 1
            smiles += '('
            branches.remove(current)
        smiles += G.nodes[current]['label']
        if current in dfs_successors:
            # Proceed to the next node in this branch
            next_nodes = dfs_successors[current]
            # ... and if needed, remember to return here later
            branches.update(next_nodes[1:])
            to_visit.extend(next_nodes)
        elif branch_depth:
            # We're finished with this branch.
            smiles += ')'
            branch_depth -= 1
    smiles += ')' * branch_depth
    return smiles

# _exp denotes surface exposed residue
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
    return isomorphism.GraphMatcher(protein_graph, generic_subgraph, node_match=node_match, edge_match=edge_match)

def set_defaults(kwargs):
    default = {'distance_cutoff':20,
               'max_degree' : 4,
                'dist_def':'COM',
                'sdef':'RSA',
                'edge_prune' : 'PERCENT',
                'percent_edges':1.0,
               'num_st_dev_edges':1.0,
               'rd_thresh':3.03,
               'rsa_thresh':0.2,
               'coef_alpha':1.0,
               'exp_beta':2.3,
               'r_offset':0.0}
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
        if (len(sg.nodes[name_node]['label']) == 1) or (len(sg.nodes[name_node]['label']) > 1 and sg.nodes[name_node]['label'][1].isdigit()):
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