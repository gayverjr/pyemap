import numpy as np
from ..data import char_to_res_name, res_name_to_char
from numpy import linalg as LA
from Bio.PDB import Superimposer
import matplotlib.pyplot as plt
import warnings
from .utils import get_edge_label, get_numerical_node_label, strip_res_number, get_graph_matcher, write_graph_smiles
import networkx as nx


def _do_fiedler_clustering(D,A,all_graphs):
    L = D - A
    eigv, eigvc = LA.eig(L)
    eigv = np.real(eigv)
    eigvc = np.real(eigvc)
    idx = eigv.argsort()
    eigv = eigv[idx]
    eigvc = eigvc[:, idx]
    # second lowest eigenvector
    eigvc2 = eigvc[:, 1][:-1]
    groups = {}
    protein_subgraphs = {}
    for i, val in enumerate(eigvc2):
        rounded_val = np.round(val, decimals=4)
        if rounded_val not in groups:
            groups[rounded_val] = [all_graphs[i]]
        else:
            graphs = groups[rounded_val]
            graphs.append(all_graphs[i])
            groups[rounded_val] = graphs
    tuples = []
    for key, val in groups.items():
        tuples.append((key, val))
    # largest groups first
    groups = {}
    tuples.sort(key=lambda x: len(x[1]), reverse=True)
    for idx, tuple1 in enumerate(tuples):
        key, group = tuple1
        pdb_list = []
        for graph in group:
            pdb_list.append(graph.graph['pdb_id'])
            id = graph.graph['pdb_id'] + "(" + str(idx + 1) + ")-" + str(
                pdb_list.count(graph.graph['pdb_id']))
            protein_subgraphs[id] = graph
        groups[idx + 1] = group
    return groups,protein_subgraphs

class FrequentSubgraph():
    '''
    Stores all information regarding a subgraph pattern identified by the gSpan algorithm.

    Attributes
    ----------
    id: str
        Unique identifier for subgraph pattern. 
    generic_subgraph: :class:`networkx.Graph`
        Graph representation of subgraph pattern found by gSpan algorithm
    support: list of str
        List of PDB IDs which contain this subgraph
    protein_subgraphs: dict of protein subgraph id (str): :class:`networkx.Graph`
        Dict which contains protein subgraphs which match this subgraph pattern. Each entry has a unique identifier and a
        :class:`networkx.Graph` derived from the graphs generated by the :class:`~pyemap.emap` class which match 
        the pattern of this subgraph pattern.
    support_number: int
        Number of PDBs this subgraph pattern was identified in 
    
    '''
    def __init__(self, G, graph_number, support, res_to_num_label, edge_thresholds):
        '''Initializes SubgraphPattern object.

        Parameters
        ----------
        G: :class:`networkx.Graph`
            Graph representation of subgraph pattern found by gSpan algorithm
        graph_number: int
            Unique numerical ID of subgraph pattern
        support: list of str
            List of PDB IDs which contain this subgraph
        '''
        self.generic_subgraph = G.copy()
        self.support = support
        self.protein_subgraphs = {}
        self.groups = {}
        self.L = []
        self.res_to_num_label = res_to_num_label
        self.edge_thresholds = edge_thresholds
        self.support_number = len(support)
        self.id = str(graph_number+1) + "_" + str(write_graph_smiles(self.generic_subgraph)) + "_" + str(self.support_number)
        if "#" in self.id:
            self.file_id = self.id.replace("#","NP")
        else:
            self.file_id = self.id
        for node in self.generic_subgraph.nodes:
            if self.generic_subgraph.nodes[node]['label'] == "#":
                self.generic_subgraph.nodes[node]['label'] = "NP"

    def general_report(self):
        ''' Generates general report which describes this subgraph pattern.
        '''
        full_str = ""
        full_str += "ID:" + str(self.id) + "\n"
        full_str += "Support:" + str(self.support_number) + "\n"
        full_str += "Where:" + str(list(self.support.keys())) + "\n"
        full_str += "Adjacency list:\n"
        G = self.generic_subgraph
        for node in G.nodes:
            full_str += G.nodes[node]['label'] + str(node) + ":["
            for neighbor in G.neighbors(node):
                full_str += G.nodes[neighbor]['label'] + str(neighbor) + "(" + str(
                    G.edges[(node, neighbor)]['num_label']) + "), "
            full_str = full_str[:-2]
            full_str += "]\n"
        return full_str

    def full_report(self):
        full_str = self.general_report()
        if len(self.protein_subgraphs) == 0:
            full_str += "Please run `find_protein_subgraphs' to get a full report.\n"
            return full_str
        full_str += str(len(self.protein_subgraphs)) + " subgraphs matching this pattern were found.\n"
        full_str += "Graphs are classified using " + self.clustering_option + " clustering.\n\n"
        for key in self.groups:
            graphs = self.groups[key]
            full_str += "Group " + str(key) + ": " + str(len(graphs)) + " members\n---------------\n"
            for graph in graphs:
                full_str += self._report_for_graph(graph) + "\n"
        return full_str

    def _report_for_graph(self, G):
        full_str = ""
        full_str += str(G.graph['id']) + "\n"
        full_str += "Nodes\n"
        for node in G.nodes:
            full_str += G.nodes[node]['label']
            full_str += " Position in alignment:" + str(G.nodes[node]['aligned_resnum'])
            full_str += "\n"
        full_str += "Adjacency list:\n"
        for node in G.nodes:
            full_str += G.nodes[node]['label'] + ":["
            for neighbor in G.neighbors(node):
                dist = '{0:.2f}'.format(G.edges[(node, neighbor)]['distance'])
                full_str += G.nodes[neighbor]['label'] + "(" + str(dist) + "), "
            full_str = full_str[:-2]
            full_str += "]\n"
        return full_str

    def _gen_node_rep(self):
        node_rep = ""
        for node, node_data in self.generic_subgraph.nodes(data=True):
            node_rep = ''.join([node_rep, node_data['label']])
        return node_rep

    def visualize_subgraph_in_ngl(self, emap, G):
        ''' Gets visualization of subgraph in ngl viewer

        Parameters
        ----------
            emap: :class:`~pyemap.emap`
                :class:`~pyemap.emap` object containing the protein subgraph
            idx: int
                Index of protein subgraph to be visualized
        '''
        colors = {"F": "orange", "Y": "blue", "W": "red", "H": "green"}
        label_texts = []
        labeled_atoms = []
        color_list = []
        selection_strs = []
        for res in G.nodes:
            label_texts.append(res)
            try:
                if res not in emap.eta_moieties:
                    color_list.append(colors[res[0]])
                    labeled_atoms.append(".CA")
                else:
                    color_list.append("pink")
                    labeled_atoms.append(next(emap.residues[res].get_atoms()).name)
            except KeyError:
                color_list.append("pink")
                labeled_atoms.append(next(emap.residues[res].get_atoms()).name)
            selection_strs.append(emap.residues[res].ngl_string)
        return label_texts, labeled_atoms, color_list, selection_strs

    def find_protein_subgraphs(self,clustering_option="structural"):
        self.groups = {}
        self.protein_subgraphs = {}
        all_graphs = []
        for pdb_id in self.support:
            all_graphs += self._find_subgraph_in_pdb(pdb_id)
        if len(all_graphs) > 1:
            D,A = self._structural_clustering(all_graphs)
            self._structural_groups,self._structural_ids = _do_fiedler_clustering(D,A,all_graphs)
            D,A = self._sequence_clustering(all_graphs)
            self._sequence_groups,self._sequence_ids = _do_fiedler_clustering(D,A,all_graphs)
            self.set_clustering(clustering_option)
        else:
            graph = all_graphs[0]
            graph.graph['id'] = graph.graph['pdb_id'] + "(" + str(1) + ")-" + str(1)
            graph.graph['group_val'] = 0.0
            self.protein_subgraphs[graph.graph['id']] = graph
            self.groups[1] = all_graphs
            self.clustering_option=clustering_option

    def set_clustering(self,clustering_option):
        if clustering_option=="structural":
            self.groups,self.protein_subgraphs = (self._structural_groups,self._structural_ids)
        elif clustering_option=="sequence":
            self.groups,self.protein_subgraphs = (self._sequence_groups,self._sequence_ids)
        else:
            raise Exception("Either structural or sequence.")
        self.clustering_option = clustering_option
        for id,graph in self.protein_subgraphs.items():
            graph.graph['id']=id

    def _sequence_clustering(self, all_graphs):
        num_graphs = len(all_graphs)
        dims = (num_graphs+1, num_graphs+1)
        D = np.zeros(dims)
        A = np.ones(dims)*0.001
        for i in range(0,num_graphs):
            A[i][i] = 0.0
        for i in range(0, len(all_graphs)):
            for j in range(i + 1, len(all_graphs)):
                G1 = all_graphs[i]
                G2 = all_graphs[j]
                G2nodes = list(G2.nodes())
                distance = 0
                for k, node1 in enumerate(G1.nodes):
                    # only count for standard amino acid residues
                    if strip_res_number(node1) in char_to_res_name:
                        node2 = G2nodes[k]
                        if str(G1.nodes[node1]['aligned_resnum']).isdigit() and str(
                                G2.nodes[node2]['aligned_resnum']).isdigit():
                            distance += np.absolute(G1.nodes[node1]['aligned_resnum'] -
                                                    G2.nodes[node2]['aligned_resnum'])
                if distance < len(G2nodes)+1:
                    A[i][j] = 1 / (distance + 1)
                    A[j][i] = 1 / (distance + 1)
            D[i][i] = np.sum(A[i])
        D[-1][-1] = np.sum(A[-1])
        return D, A
    # J. Mol. Biol. (1999) 292, 441-464
    # This is known as Fiedler eigenvalue, or Spectral Graph Partitioning
    def _structural_clustering(self, all_graphs):
        num_graphs = len(all_graphs)
        dims = (num_graphs+1, num_graphs+1)
        D = np.zeros(dims)
        A = np.ones(dims)*0.001
        for i in range(0,num_graphs):
            A[i][i] = 0.0
        for i in range(0, num_graphs):
            for j in range(i + 1, num_graphs):
                distance = self._subgraph_rmsd(all_graphs[i], all_graphs[j])
                if distance < 1:
                    A[i][j] = 1 / (distance + 1)
                    A[j][i] = 1 / (distance + 1)
            D[i][i] = np.sum(A[i])
        D[-1][-1] = np.sum(A[-1])
        return D, A

    def subgraph_rmsd(self,sg1,sg2):
        return self._subgraph_rmsd(self.protein_subgraphs[sg1],self.protein_subgraphs[sg2])

    def _subgraph_rmsd(self, sg1, sg2):
        emap1 = self.support[sg1.graph['pdb_id']]
        emap2 = self.support[sg2.graph['pdb_id']]
        atoms1 = []
        atoms2 = []
        nodes1 = list(sg1.nodes)
        nodes2 = list(sg2.nodes)
        for i in range(0, len(nodes1)):
            res1 = emap1.residues[nodes1[i]]
            res2 = emap2.residues[nodes2[i]]
            if 'CA' in res1 and 'CA' in res2:
                atoms1.append(res1['CA'])
                atoms2.append(res2['CA'])
            else:
                shared_id = None
                for atm in res1:
                    if atm.id in res2:
                        shared_id = atm.id
                        break
                if shared_id is not None:
                    atoms1.append(res1[shared_id])
                    atoms2.append(res2[shared_id])
                else:
                    return 100
        if len(atoms1) == len(nodes1):
            si = Superimposer()
            si.set_atoms(atoms1, atoms2)
            return si.rms
        else:
            return 100

    def _find_subgraph_in_pdb(self,pdb_id):
        GM = get_graph_matcher(self.support[pdb_id].init_graph, self.generic_subgraph)
        subgraph_isos = GM.subgraph_monomorphisms_iter()
        sgs = []
        degree_dicts = []
        for mapping in subgraph_isos:
            sg = self._generate_protein_subgraph(mapping, self.support[pdb_id].init_graph, self.generic_subgraph, self.support[pdb_id])
            degree_dict = dict(sg.degree)
            if degree_dict not in degree_dicts:
                degree_dicts.append(degree_dict)
                sgs.append(sg)
        return sgs

    def _generate_protein_subgraph(self, mapping, protein_graph, generic_subgraph, emap_obj):
        mapping = dict((v, k) for k, v in mapping.items())
        protein_subgraph = generic_subgraph.copy()
        protein_subgraph = nx.relabel_nodes(protein_subgraph, mapping)
        for node in protein_subgraph.nodes():
            protein_subgraph.nodes[node]['shape'] = protein_graph.nodes[node]['shape']
            protein_subgraph.nodes[node]['label'] = str(node)
            protein_subgraph.nodes[node]['aligned_resnum'] = emap_obj.residues[node].aligned_residue_number
            protein_subgraph.graph['pdb_id'] = protein_graph.graph['pdb_id']
        for edge in protein_subgraph.edges():
            for key in protein_graph.edges[edge]:
                protein_subgraph.edges[edge][key] = protein_graph.edges[edge][key]
        return protein_subgraph


        

