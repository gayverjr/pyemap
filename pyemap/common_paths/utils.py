from ..data import char_to_res_name
from networkx.algorithms import isomorphism

def get_edge_label(G, edge, edge_thresholds):
    dist = G.edges[edge]['distance']
    try:
        label = 2
        for thresh in edge_thresholds:
            if dist < thresh:
                break
            else:
                label += 1
        return label
    except:
        return 2

def strip_insertion_code(u):
    for i in range(0, len(u)):
        if u[i].isalpha():
            return u[:i-1]
    return u

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
    total_edges = set(map(frozenset, G.edges))
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
        if current in predecessors:
            # It's not the first atom we're visiting, so we want to see if the
            # edge we last crossed to get here is interesting.
            previous = predecessors[current]
            assert len(previous) == 1
            previous = previous[0]
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

