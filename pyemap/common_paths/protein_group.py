import numpy as np
from io import StringIO
import os
from ..process_data import process
from collections import OrderedDict
import networkx as nx
from .find_subgraph import find_sg
import matplotlib.pyplot as plt


def strip_res_number(u):
    digit_idx=-1
    for i in range(0,len(u)):
        if u[i].isdigit():
            return u[:i]

def node_match(node1,node2):
    return node1['label'] == node2['label']
    
def edge_match(edge1,edge2):
    return edge1['label'] == edge2['label']

class Subgraph():
    def __init__(self,subgraph,graph_id,occurences):
        self.G = subgraph
        self.occurences = occurences
        self.node_rep = self.gen_node_rep()
        self.support = len(occurences)
        #self.id = id
        self.id = str(graph_id) + "_"+str(self.node_rep) + "_" + str(self.support)

    def gen_node_rep(self):
        node_rep = ""
        for node,node_data in self.G.nodes(data=True):
            node_rep = ''.join([node_rep, node_data['label']])
        return node_rep
        
class protein_group():

    def _set_edge_labels(self,edge_thresholds):
        if edge_thresholds == None:
            edge_thresholds = [8,12]
        self.edge_thresholds = edge_thresholds
    
    def _set_node_labels(self,node_labels,categories,surface_exposed):
        if node_labels == None:
            self.res_to_num_label = { "W": 2, "Y": 3, "H": 4, "F": 5}
            self.num_label_to_res = { 2:"W", 3:"Y", 4:"H", 5:"F"}
        else:
            self.res_to_num_label = node_labels
            self.num_label_to_res = { 2:"W", 3:"Y", 4:"H", 5:"F"}
            num_label = 5
            for category in categories:
                num_label+=1
                self.num_label_to_res[num_label]=category
            self.num_label_to_res[num_label+1] = "NP"

    def get_edge_label(self,G,edge):
        dist = G.edges[edge]['distance']
        label = 2
        for thresh in self.edge_thresholds:
            if dist < thresh:
                return label
            else:
                label+=1
        return label

    def get_numerical_node_label(self,u):
        if u in self.res_to_num_label:
            result = self.res_to_num_label[u]
        elif strip_res_number(u) in ["W","Y","H","F"]:
            res_name = strip_res_number(u)
            result = self.res_to_num_label[res_name]
        else:
            result = len(self.num_label_to_res)+1
        return result

    def __init__(self,title,temp_dir=""):
        self.title = title
        self.emaps = OrderedDict()
        self.parameters = {}
        self.raw_gspan_output = StringIO()
        self.temp_dir = temp_dir
        self.subgraphs = {}
    
    def add_emap(self,emap_obj):
        if emap_obj.pdb_id not in self.emaps:
            self.emaps[emap_obj.pdb_id] = emap_obj
            print("Added "+ emap_obj.pdb_id)
        else:
            print("Duplicate!")
    
    def set_parameters(self,parameters):
        self.parameters = parameters

    def find_subgraph(self,graph_id,emap_id):
        graph_id = str(graph_id)
        subgraph = self.subgraphs[graph_id]
        if emap_id in subgraph.occurences:
            G = self.subgraphs[graph_id].G
            src = G.nodes[0]
            sgs = self.find_sg(self.emaps[emap_id].init_graph,G,src)
        return sgs

    def process_emap(self,emap_name,eta_moieties="All",custom_str=""):
        my_emap = self.emaps[emap_name]
        print("Processing:" + emap_name)
        process(my_emap,dist_def=1)

    def generate_graph_db(self,node_labels=None,categories=None,edge_thresholds=None,surface_exposed=False):
        self._set_node_labels(node_labels,categories,surface_exposed)
        self._set_edge_labels(edge_thresholds)
        f = open(os.path.join(self.temp_dir,'graphdatabase.txt'), "w")
        print(os.path.join(self.temp_dir,'graphdatabase.txt'))
        for i,key in enumerate(self.emaps):
            G = self.emaps[key].init_graph
            f.write("t # "+str(i)+"\n")
            for i,node in enumerate(G.nodes):
                f.write("v "+str(i) + " " + str(self.get_numerical_node_label(node))+"\n")
            for i,edge in enumerate(G.edges):
                f.write("e " + str(list(G.nodes()).index(edge[0]))+ " " + str(list(G.nodes()).index(edge[1]))+ " " + str(self.get_edge_label(G,edge))+ "\n")
        f.write("t # -1")

    def run_gspan(self,support=10,lower_bound=4):
        import sys
        # change standard output temporarily
        old_stdout = sys.stdout
        #sys.stdout = self.raw_gspan_output
        f = open(os.path.join(self.temp_dir,'gspan_results.out'), "w")
        sys.stdout = f
        from gspan_mining.config import parser
        from gspan_mining.main import main
        args_str = '-s '+ str(support) + ' -d False -l '+str(lower_bound)+ ' -p False -w True ' + str(os.path.join(self.temp_dir,'graphdatabase.txt'))
        FLAGS, _ = parser.parse_known_args(args=args_str.split())
        gs = main(FLAGS)
        # give us our old standard output back
        sys.stdout = old_stdout
        f.close()
        self.prune_gspan()

    def prune_gspan(self):
        buff = open(os.path.join(self.temp_dir,'gspan_results.out'), "r")
        f = open(os.path.join(self.temp_dir,'pruned_results.out'), "w")
        subgraphs = []
        buff.seek(0)
        lines = buff.readlines()
        line_idx = 0
        while line_idx < len(lines):
            line = lines[line_idx]
            if len(line.split())==3 and line.split()[0]=="t" and line.split()[1]=="#":
                graph_id = line.split()[2]
                line_idx+=1
                start_idx = line_idx-1
                G = nx.Graph()
                line = lines[line_idx]
                while "---" not in line:
                    line = lines[line_idx]
                    if len(line.split())>1 and line.split()[0]=="v":
                        node_idx = int(line.split()[1])
                        node_label = int(line.split()[2])
                        G.add_node(node_idx)
                        G.nodes[node_idx]['label']= self.num_label_to_res[node_label]
                    if len(line.split())>1 and line.split()[0]=="e":
                        idx1 = int(line.split()[1])
                        idx2 = int(line.split()[2])
                        edge_label = int(line.split()[3])
                        G.add_edge(idx1,idx2,label=edge_label)
                    if "where" in line: 
                        for i in range(start_idx,line_idx+1):
                            f.write(lines[i])
                        f.write("----------\n")
                        where_list = line[7:-2].strip('][').split(', ')
                        where_list = list(np.array(where_list,dtype=int))
                        emap_list = list(self.emaps.keys())
                        occurences = []
                        for where_idx in where_list:
                            occurences.append(emap_list[where_idx])
                        subgraphs.append(Subgraph(G,graph_id,occurences))
                        #self.subgraphs[sg.id] = sg
                    line_idx+=1
            line_idx+=1
        buff.close()
        subgraphs.sort(key=lambda x: x.support, reverse=True)
        for sg in subgraphs:
            self.subgraphs[sg.id] = sg

    def generate_candidate_subgraph(self,graph):
        G = nx.Graph()
        has_FAD= False
        for i,node in enumerate(graph.nodes):
            G.add_node(i)
            num_label = self.get_numerical_node_label(node)
            G.nodes[i]['label']= self.num_label_to_res[num_label]
        for i,edge in enumerate(graph.edges):
            node1_idx = list(graph.nodes()).index(edge[0])
            node2_idx = list(graph.nodes()).index(edge[1])
            edge_label = graph.edges[edge]['label']
            G.add_edge(node1_idx,node2_idx,label=edge_label)
        return G

    def get_candidates(self,G,source):
        # returns networkx node
        candidates = []
        for node in G.nodes:
            if strip_res_number(node)==source['label']:
                candidates.append(node)
        return candidates

#non-recursive
    def dfs_nr(self,full_graph,target_subgraph,src):
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
                    edge_label = self.get_edge_label(full_graph,edge)
                    cur_G.add_node(cur_node)
                    cur_G.add_edge(prev_node,cur_node,label=edge_label,distance=full_graph.edges[edge]['distance'])
                    cur_num_nodes = len(list(cur_G.nodes))
                    cur_num_edges = len(list(cur_G.edges))
                    if target_num_nodes == cur_num_nodes and target_num_edges == cur_num_edges:
                        test_graph = self.generate_candidate_subgraph(cur_G)
                        if nx.is_isomorphic(test_graph,target_subgraph,edge_match=edge_match,node_match=node_match):
                            found_subgraphs.append(cur_G)
                    elif target_num_nodes > cur_num_nodes:
                        stack.append(cur_node)
                        graph_stack.append(cur_G)
        return found_subgraphs

    def find_sg(self,graph,subgraph,source):
        # graph is a networkx object generated by emap
        # subgraph is networkx object with node labels in the same format as our gspan stuff
        # source is a numerical label of the starting node
        # residue_map is a dict mapping node labels to residue types
        source_candidates = self.get_candidates(graph,source)
        sgs = []
        for src in source_candidates:
            sgs += self.dfs_nr(graph,subgraph,src)
        return sgs










