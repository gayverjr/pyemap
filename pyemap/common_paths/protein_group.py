import numpy as np
from io import StringIO
import os
from ..process_data import process
from collections import OrderedDict
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import isomorphism
import time
from ..data import char_to_res_name


def strip_res_number(u):
    for i in range(0, len(u)):
        if u[i].isdigit():
            return u[:i]


def node_match(node1, node2):
    return node1['num_label'] == node2['num_label']


def edge_match(edge1, edge2):
    return edge1['num_label'] == edge2['num_label']


class FrequentSubgraph():
    def __init__(self, G, graph_number, support):
        self.generic_subgraph = G
        self.support = support
        self.specific_subgraphs = {}
        self.node_rep = self.gen_node_rep()
        self.support_number = len(support)
        self.id = str(graph_number) + "_" + str(self.node_rep) + "_" + str(self.support_number)

    def generate_generic_report(self):
        full_str = ""
        a_str=""
        a_str+="Generic Adjacency List: \n"
        for node in self.G.nodes:
            main_node=(self.G.nodes[node]['label'])
            neighborhood=(list(self.G.neighbors(node)))
            a_str+=main_node + str(node)+ "["
            for l in range(0,len(neighborhood)):
                if l>0:
                    a_str+= ","  + str(neighborhood[l])
                    if l==0:
                        a_str+= str(neighborhood[l])
                a_str+="]"  + "\n"
        a_str+="\n"
        a_str+= "Support=" + str(len(self.occurences))
        a_str+="\n"
        for i in range(len(self.occurences)):
            full_str+="\n"
            # print(i)
            full_str+= str(self.occurences[i]) +"\n" #first pdb id
            a_str+= "PDB ID: " + str(self.occurences[i]) + "\n"
        # grab all specific graphs for first pdb
        return  a_str

    def specific_subgraph_example(self):
        full_str = ""
        a_str=""
        # print(len(self.occurences))
        a_str+= "Support=" + str(len(self.occurences))
        a_str+="\n"
        for i in range(len(self.occurences)):
            full_str+="\n"
            full_str+= str(self.occurences[i]) +"\n" #first pdb id
            a_str+="\n"
            a_str+= "PDB ID: " + str(self.occurences[i]) + "\n"
            # grab all specific graphs for first pdb
            specific_graphs_for_first_pdb = self.specific_graphs[i]
            for j in range (len(specific_graphs_for_first_pdb)):
                full_str+="\n"
                a_str+="\n"
                a_str+= "Adjacency List: "
                a_str+= "PDB Specific Subgraph " + str(j) +"\n"
                G = specific_graphs_for_first_pdb[j] #first specific graph
                for edge in G.edges:
                    node1_label = G.nodes[edge[0]]['label']
                    node2_label = G.nodes[edge[1]]['label']
                for node in G.nodes:
                    main_node=(G.nodes[node]['label'])
                    neighborhood=(list(G.neighbors(node)))
                    a_str+=main_node +"[ "
                    for l in range(0,len(neighborhood)):
                        a_str+= neighborhood[l] + ":" + str(round((G.edges[(node, neighborhood[l])]['distance']),4))+","
                    a_str+="]"  + "\n"
                    full_str+="Node1:"+str(node1_label)+", Node2:"+str(node2_label)+": Distance:" + str(G.edges[edge]['distance']) + "\n"
            return  a_str

    def gen_node_rep(self):
        node_rep = ""
        for node, node_data in self.generic_subgraph.nodes(data=True):
            node_rep = ''.join([node_rep, node_data['label']])
        return node_rep

    def visualize_subgraph_in_ngl(self, emap, idx):
        colors = {"F": "orange", "Y": "blue", "W": "red", "H": "green"}
        label_texts = []
        labeled_atoms = []
        color_list = []
        selection_strs = []
        G = self.specific_subgraphs[emap.pdb_id][idx]
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


class PDBGroup():
    def __init__(self, title, temp_dir=""):
        self.title = title
        self.emaps = OrderedDict()
        self.temp_dir = temp_dir
        self.frequent_subgraphs = {}
        self.res_to_num_label = {}
        self.num_label_to_res = {}
        self.edge_thresholds = []

    # chains and eta_moieties should be dictionaries
    def process_emaps(self, chains=None, eta_moieties=None, **kwargs):
        if not chains:
            chains = {}
            for pdb_id in self.emaps:
                chains[pdb_id] = [self.emaps[pdb_id].chains[0]]
        for pdb_id in self.emaps:
            if eta_moieties:
                cur_eta_moieties = eta_moieties[pdb_id]
            else:
                cur_eta_moieties = []
            process(self.emaps[pdb_id], chains=chains[pdb_id], eta_moieties=cur_eta_moieties, **kwargs)

    def _set_edge_labels(self, edge_thresholds):
        if edge_thresholds == None:
            edge_thresholds = [8.0, 12.0]
        self.edge_thresholds = edge_thresholds

    def _set_node_labels(self, node_labels, categories):
        if node_labels == None:
            self.res_to_num_label = {"W": 2, "Y": 3, "H": 4, "F": 5, "NP": 6}
            self.num_label_to_res = {2: "W", 3: "Y", 4: "H", 5: "F", 6: "NP"}
        else:
            self.res_to_num_label = node_labels
            self.num_label_to_res = {2: "W", 3: "Y", 4: "H", 5: "F"}
            num_label = 5
            for category in categories:
                num_label += 1
                self.num_label_to_res[num_label] = category
            self.num_label_to_res[num_label + 1] = "NP"
            self.res_to_num_label["NP"] = num_label + 1

    def get_edge_label(self, G, edge):
        dist = G.edges[edge]['distance']
        label = 2
        for thresh in self.edge_thresholds:
            if dist < thresh:
                return label
            else:
                label += 1
        return label

    def get_numerical_node_label(self, u):
        if u in self.res_to_num_label:
            result = self.res_to_num_label[u]
        elif strip_res_number(u) in char_to_res_name:
            res_name = strip_res_number(u)
            result = self.res_to_num_label[res_name]
        else:
            result = len(self.num_label_to_res) + 1
        return result

    def add_emap(self, emap_obj):
        if emap_obj.pdb_id not in self.emaps:
            self.emaps[emap_obj.pdb_id] = emap_obj
            print("Added emap object with PDB ID: " + emap_obj.pdb_id)
        else:
            print("An emap object with PDB ID:" + str(emap_obj.pdb_id) + " is already in the data set. Skipping...")

    def generate_graph_database(self, node_labels=None, categories=None, edge_thresholds=None):
        # remove all old data
        self.frequent_subgraphs = {}
        self.res_to_num_label = {}
        self.num_label_to_res = {}
        self.edge_thresholds = []
        #-------------------------
        self._set_node_labels(node_labels, categories)
        self._set_edge_labels(edge_thresholds)
        f = open(os.path.join(self.temp_dir, 'graphdatabase.txt'), "w")
        for i, key in enumerate(self.emaps):
            G = self.emaps[key].init_graph
            f.write("t # " + str(i) + "\n")
            for i, node in enumerate(G.nodes):
                f.write("v " + str(i) + " " + str(self.get_numerical_node_label(node)) + "\n")
            for i, edge in enumerate(G.edges):
                f.write("e " + str(list(G.nodes()).index(edge[0])) + " " + str(list(G.nodes()).index(edge[1])) + " " +
                        str(self.get_edge_label(G, edge)) + "\n")
        f.write("t # -1")

    def run_gspan(self, support=10, lower_bound=4):
        import sys
        old_stdout = sys.stdout
        f = open(os.path.join(self.temp_dir, 'gspan_results.out'), "w")
        sys.stdout = f
        from gspan_mining.config import parser
        from gspan_mining.main import main
        args_str = '-s ' + str(support) + ' -d False -l ' + str(lower_bound) + ' -p False -w True ' + str(
            os.path.join(self.temp_dir, 'graphdatabase.txt'))
        FLAGS, _ = parser.parse_known_args(args=args_str.split())
        gs = main(FLAGS)
        # give us our old standard output back
        sys.stdout = old_stdout
        f.close()
        self.generate_frequent_subgraphs()

    def generate_frequent_subgraphs(self):
        buff = open(os.path.join(self.temp_dir, 'gspan_results.out'), "r")
        subgraphs = []
        lines = buff.readlines()
        line_idx = 0
        while line_idx < len(lines):
            line = lines[line_idx]
            if len(line.split()) == 3 and line.split()[0] == "t" and line.split()[1] == "#":
                graph_number = line.split()[2]
                line_idx += 1
                start_idx = line_idx - 1
                G = nx.Graph()
                line = lines[line_idx]
                while "---" not in line:
                    line = lines[line_idx]
                    if len(line.split()) > 1 and line.split()[0] == "v":
                        node_idx = int(line.split()[1])
                        node_label = int(line.split()[2])
                        G.add_node(node_idx)
                        G.nodes[node_idx]['label'] = self.num_label_to_res[node_label]
                        G.nodes[node_idx]['num_label'] = node_label
                    if len(line.split()) > 1 and line.split()[0] == "e":
                        idx1 = int(line.split()[1])
                        idx2 = int(line.split()[2])
                        edge_label = int(line.split()[3])
                        G.add_edge(idx1, idx2, label=edge_label)
                        G.edges[(idx1, idx2)]['num_label'] = edge_label
                    if "where" in line:
                        pdb_list_by_index = line[7:-2].strip('][').split(', ')
                        pdb_list_by_index = list(np.array(pdb_list_by_index, dtype=int))
                        pdb_list = list(self.emaps.keys())
                        support = []
                        for idx in pdb_list_by_index:
                            support.append(pdb_list[idx])
                        subgraphs.append(FrequentSubgraph(G, graph_number, support))
                    line_idx += 1
            line_idx += 1
        buff.close()
        subgraphs.sort(key=lambda x: x.support_number, reverse=True)
        for sg in subgraphs:
            for pdb_id in sg.support:
                sg.specific_subgraphs[pdb_id] = self.find_subgraph_in_pdb(sg, pdb_id)
            self.frequent_subgraphs[sg.id] = sg

    def generate_specific_subgraph(self, mapping, protein_graph, generic_subgraph):
        mapping = dict((v, k) for k, v in mapping.items())
        specific_subgraph = generic_subgraph.copy()
        specific_subgraph = nx.relabel_nodes(specific_subgraph, mapping)
        for node in specific_subgraph.nodes():
            specific_subgraph.nodes[node]['shape'] = protein_graph.nodes[node]['shape']
            specific_subgraph.nodes[node]['label'] = str(node)
        for edge in specific_subgraph.edges():
            for key in protein_graph.edges[edge]:
                specific_subgraph.edges[edge][key] = protein_graph.edges[edge][key]
        return specific_subgraph

    def find_subgraph_in_pdb(self, subgraph, pdb_id):
        generic_subgraph = subgraph.generic_subgraph
        if pdb_id not in subgraph.support:
            return []
        protein_graph = self.emaps[pdb_id].init_graph
        for node in protein_graph.nodes:
            protein_graph.nodes[node]['num_label'] = self.get_numerical_node_label(node)
        for edge in protein_graph.edges:
            protein_graph.edges[edge]['num_label'] = self.get_edge_label(protein_graph, edge)
        GM = isomorphism.GraphMatcher(protein_graph, generic_subgraph, node_match=node_match, edge_match=edge_match)
        subgraph_isos = GM.subgraph_monomorphisms_iter()
        sgs = []
        for mapping in subgraph_isos:
            sg = self.generate_specific_subgraph(mapping, protein_graph, generic_subgraph)
            sgs.append(sg)
        return sgs
