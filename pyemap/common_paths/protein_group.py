import numpy as np
from io import StringIO
import os
from ..process_data import process
from collections import OrderedDict
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import isomorphism
import time
from fpdf import FPDF



  
# save the pdf with name .pdf 
  
# imagelist is the list with all image filenames



def strip_res_number(u):
    for i in range(0,len(u)):
        if u[i].isdigit():
            return u[:i]

def node_match(node1,node2):
    return node1['num_label'] == node2['num_label']
    
def edge_match(edge1,edge2):
    return edge1['num_label'] == edge2['num_label']

class Subgraph():
    def __init__(self,subgraph,graph_id,occurences):
        self.G = subgraph
        self.occurences = occurences
        self.specific_graphs = []
        self.node_rep = self.gen_node_rep()
        self.support = len(occurences)
        #self.id = id
        self.id = str(graph_id) + "_"+str(self.node_rep) + "_" + str(self.support)
 
    def generate_generic_report(self):
        full_str = ""
        full_str+= "ID:" + str(self.id) + "\n"
        full_str+= "Support:" + str(self.support) + "\n"
        full_str+= "Nodes:\n"
        for node in self.G.nodes:
            full_str+= self.G.nodes[node]['label']+","
            full_str+= "\nEdges:\n"
        for edge in self.G.edges:
            node1_label = self.G.nodes[edge[0]]['label']
            node2_label = self.G.nodes[edge[1]]['label']
            full_str+="Node1:"+str(node1_label)+", Node2:"+str(node2_label)+": Edge label:" + str(self.G.edges[edge]['num_label']) +"\n"
        return full_str 


    def specific_subgraph_example(self):
        full_str = ""
        a_str=""
       # print(len(self.occurences))
        a_str+= "Support=" + str(len(self.occurences))
        a_str+="\n"
        for i in range(len(self.occurences)):
            full_str+="\n"  
            
           # print(i)
            full_str+= str(self.occurences[i]) +"\n" #first pdb id
            a_str+="\n"
            a_str+= "PDB ID: " + str(self.occurences[i]) + "\n" 
            
    # grab all specific graphs for first pdb
            specific_graphs_for_first_pdb = self.specific_graphs[i]

    # Important: this G is different than the self.G in the above example,
    #  which is the generic version
    #  this guy is the networkx graph for the specific instance of the subgraph found in a pdb
    #  which contains all the info about the nodes in the specific pdb
            #print(len(specific_graphs_for_first_pdb))
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
                        if l==len(neighborhood) -1:
                            a_str+= neighborhood[l] + ":" + str(round((G.edges[(node, neighborhood[l])]['distance']),4))
                        else:
                            a_str+= neighborhood[l] + ":" + str(round((G.edges[(node, neighborhood[l])]['distance']),4)) +","
                        a_str+="]"  + "\n"  



                    full_str+="Node1:"+str(node1_label)+", Node2:"+str(node2_label)+": Distance:" + str(G.edges[edge]['distance']) + "\n"
            
        return  a_str


        

    def gen_node_rep(self):
        node_rep = ""
        for node,node_data in self.G.nodes(data=True):
            node_rep = ''.join([node_rep, node_data['label']])
        return node_rep
        print(node_rep)
    
    def visualize_subgraph_in_ngl(self,emap,idx):
        colors = {"F": "orange", "Y": "blue", "W": "red", "H": "green"}
        label_texts = []
        labeled_atoms = []
        color_list = []
        selection_strs = []
        G = self.specific_graphs[self.occurences.index(emap.pdb_id)][idx]
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
        return label_texts,labeled_atoms,color_list,selection_strs


    



class protein_group():

    def _set_edge_labels(self,edge_thresholds):
        if edge_thresholds == None:
            edge_thresholds = [8.0,12.0]
        self.edge_thresholds = edge_thresholds



        

       
    
    def _set_node_labels(self,node_labels,categories,surface_exposed):
        if node_labels == None:
            self.res_to_num_label = {"W":2, "Y": 3, "H": 4, "F": 5, "NP":6}
            self.num_label_to_res = {2:"W", 3:"Y", 4:"H", 5:"F", 6:"NP"}
        else:
            self.res_to_num_label = node_labels
            self.num_label_to_res = { 2:"W", 3:"Y", 4:"H", 5:"F"}
            num_label = 5
            for category in categories:
                num_label+=1
                self.num_label_to_res[num_label]=category
            self.num_label_to_res[num_label+1] = "NP"
            self.res_to_num_label["NP"] = num_label+1


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


    def find_subgraph(self,graph_id,emap_id):
        graph_id = str(graph_id)
        subgraph = self.subgraphs[graph_id]
        sgs = []
        if emap_id in subgraph.occurences:
            sgs = self.find_subgraph_in_pdb(self.emaps[emap_id].init_graph,subgraph,emap_id)
        return sgs


    def generate_graph_db(self,node_labels=None,categories=None,edge_thresholds=None,surface_exposed=False):
        self._set_node_labels(node_labels,categories,surface_exposed)
        self._set_edge_labels(edge_thresholds)
        f = open(os.path.join(self.temp_dir,'graphdatabase.txt'), "w")
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
        
        self.generate_subgraphs()
        

    def generate_subgraphs(self):
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
                        G.nodes[node_idx]['num_label'] = node_label
                    if len(line.split())>1 and line.split()[0]=="e":
                        idx1 = int(line.split()[1])
                        idx2 = int(line.split()[2])
                        edge_label = int(line.split()[3])
                        G.add_edge(idx1,idx2,label=edge_label)
                        G.edges[(idx1,idx2)]['num_label']= edge_label
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
                    line_idx+=1
            line_idx+=1
        buff.close()
        subgraphs.sort(key=lambda x: x.support, reverse=True)
        for sg in subgraphs:
            self.subgraphs[sg.id] = sg
            occurences = sg.occurences
            for emap_id in occurences:
                specific_subgraphs = self.find_subgraph(sg.id,emap_id)
                self.subgraphs[sg.id].specific_graphs.append(specific_subgraphs)

            

    def generate_specific_subgraph(self,mapping,graph,subgraph):
        mapping = dict((v, k) for k, v in mapping.items())
        specific_graph = subgraph.copy()
        specific_graph = nx.relabel_nodes(specific_graph,mapping)
        for node in specific_graph.nodes():
            for key in graph.nodes[node]:
                specific_graph.nodes[node]['shape'] = graph.nodes[node]['shape']
                specific_graph.nodes[node]['label'] = graph.nodes[node]['label']
                specific_graph.nodes[node]['num_label'] = self.get_numerical_node_label(node)
        for edge in specific_graph.edges():
            for key in graph.edges[edge]:
                specific_graph.edges[edge][key] = graph.edges[edge][key]
            specific_graph.edges[edge]['num_label'] = self.get_edge_label(specific_graph,edge) 
        return specific_graph



    def find_subgraph_in_pdb(self,graph,subgraph,emap_id):
        for node in graph.nodes:
            num_label = self.get_numerical_node_label(node)
            graph.nodes[node]['num_label'] = num_label
            graph.nodes[node]['label'] = str(node)
        for node in subgraph.G.nodes():
            node_label = subgraph.G.nodes[node]['label']
            num_label = self.res_to_num_label[node_label]
            subgraph.G.nodes[node]['num_label'] = num_label
        for edge in graph.edges:
            graph.edges[edge]['num_label'] = self.get_edge_label(graph,edge)
        GM = isomorphism.GraphMatcher(graph, subgraph.G,node_match=node_match,edge_match=edge_match)
        subgraph_isos = GM.subgraph_monomorphisms_iter()
        sgs = []
        for mapping in subgraph_isos:
            sg = self.generate_specific_subgraph(mapping,graph,subgraph.G)
            sgs.append(sg)
        return sgs     

    
    
    
        
     






