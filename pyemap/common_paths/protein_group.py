import numpy as np
import os
import shutil
from ..process_data import process
import networkx as nx
import time
import datetime
from ..data import res_name_to_char, char_to_res_name
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1
from .frequent_subgraph import FrequentSubgraph
import warnings
from .utils import get_edge_label, get_numerical_node_label, get_graph_matcher, strip_insertion_code

def nodes_and_edges_from_string(graph_str, edge_thresholds, residue_categories):
    graph_str = graph_str.replace(" ", "")
    graph_str = list(graph_str)
    node_list = []
    idx = 0
    while idx < len(graph_str):
        node_list.append(str(graph_str[idx]))
        idx+=1
    l1 = []
    for i in range(0, len(edge_thresholds)):
        l1.append(i + 2)
    from itertools import product
    indices = [i for i,x in enumerate(node_list) if x=="*"]
    if len(indices) == 0:
        node_combs = [node_list]
    else:
        node_combs = []
        wildcard_combs = list(product(residue_categories,repeat=len(indices)))
        for comb in wildcard_combs:
            for i,idx in enumerate(indices):
                node_list[idx] = comb[i]
            node_combs.append(node_list.copy())
    edge_combs = list(product(l1, repeat=len(node_list) - 1))
    return node_combs, edge_combs

class PDBGroup():
    '''
    Contains all information regarding the group of proteins being analyzed, and all of the 
    the subgraph patterns identified by the gSpan algorithm.

    Attributes
    ----------
    title: str
        Title of PDB group 
    emaps: dict of str: :class:`~pyemap.emap`
        Dict of PDBs being analyzed by PyeMap. The keys are PDB IDs, meaning that only one :class:`~pyemap.emap` object per PDB ID is allowed.
    temp_dir: str
        Path where temporary files are to be stored   
    subgraph_patterns: dict of str: :class:`~pyemap.common_paths.FrequentSubgraph`
        Dict of subgraph patterns found by GSpan. Keys are the unique IDs of the :class:`~pyemap.common_paths.FrequentSubgraph` objects.

    '''
    def __init__(self, title, temp_dir="."):
        ''' Initializes PDBGroup object
        
        Parameters
        ----------
        title: str
            Title of PDB group 
        temp_dir: str
            Path where temporary files are to be stored 
        '''
        self.title = title
        self.emaps = {}
        self.temp_dir = temp_dir
        self.subgraph_patterns = {}
        self.node_labels = {}
        self.residue_categories = {}
        self.edge_thresholds = []
        self.emap_parameters = {}
        self.gspan_parameters = {}
        self.graph_database_parameters = {}
        self.included_eta_moieties = {}
        self.included_chains = {}
        self.included_standard_residues = []
        self.sequences = {}
        self.aligned_sequences = {}

    def _clean_graph_database(self):
        self.graph_database_parameters = {}
        self.node_labels = {}
        self.residue_categories = {}
        self.edge_thresholds = []
        self.sep_buried_exposed = False
        self._clean_subgraphs()
    
    def _clean_subgraphs(self):
        self.gspan_parameters = {}
        self.subgraph_patterns = {}

    def _reset_process(self):
        self.emap_parameters = {}
        self.included_eta_moieties = {}
        self.included_chains = {}
        self.included_standard_residues = []
        self.sequences = {}
        self.aligned_sequences = {}
        self._clean_graph_database()

    def _align_sequences(self, chains):
        records = []
        valid_ids = []
        inp = os.path.join(self.temp_dir, "data.fasta")
        with open(inp, 'w') as handle:
            for pdb_id, emap in self.emaps.items():
                for chain in emap.active_chains:
                    handle.write(emap.sequences[chain] + "\n")
        out = os.path.join(self.temp_dir, "data_aligned.fasta")
        log = os.path.join(self.temp_dir, "log.txt")
        try:
            muscle_cline = MuscleCommandline(input=inp, out=out, log=log)
            muscle_cline()
        except Exception as e:
            shutil.copyfile(inp, out)
        seqIO = SeqIO.parse(out, "fasta")
        for record in seqIO:
            self.aligned_sequences[record.id.upper()] = record.seq
        # now lets save the updated sequence numbers
        for pdb_id, emap in self.emaps.items():
            for chain, residues in emap.active_chains.items():
                original_idx = emap.chain_start[chain]
                aligned_idx = 1
                seq_map = {}
                if pdb_id + ":" + chain in self.aligned_sequences:
                    aligned_seq = self.aligned_sequences[pdb_id + ":" + chain]
                    for res in aligned_seq:
                        if not res == "-":
                            seq_map[original_idx] = aligned_idx
                            original_idx += 1
                        aligned_idx += 1
                for residue in residues:
                    resnum = residue.id[1]
                    if int(resnum) in seq_map:
                        residue.aligned_residue_number = seq_map[int(resnum)]
                    else:
                        residue.aligned_residue_number = 'X'

    def process_emaps(self, chains={}, eta_moieties={}, sdef=None, include_residues=["TYR", "TRP"], **kwargs):
        ''' Processes :class:`~pyemap.emap` objects in order to generate protein graphs. 
        
        For a list of accepted kwargs, see the documentation for :func:`~pyemap.process_data.process`.

        Parameters
        -----------
        chains: dict of str: list of str, optional
            Chains to include for each PDB
        eta_moieties: dict of str: list of str, optional
            Dict containing list of ETA moieties(specified by their residue label) to include for each PDB
        
        Examples
        --------
        >>> my_pg = pyemap.common_paths.PDBGroup()
        >>> # Add pdbs 1u3d,1u3c,6PU0,4I6G,2J4D ...
        >>> eta_moieties = {'1u3d': ['FAD510(A)-2'], '1u3c': ['FAD510(A)-2'], '6PU0': ['FAD501(A)-2'], '4I6G': ['FAD900(A)-2'], '2J4D': ['FAD1498(A)-2']}
        >>> chains = {'1u3d': ['A'], '1u3c': ['A'], '6PU0': ['A'], '4I6G': ['A'], '2J4D': ['A']}
        >>> my_pg.process_emaps(chains=chains,eta_moieties=eta_moieties)

        '''
        start_time = time.time()
        self._reset_process()
        remove_pdbs = []
        for pdb_id in self.emaps:
            if pdb_id not in eta_moieties:
                eta_moieties[pdb_id] = list(self.emaps[pdb_id].eta_moieties.keys())
            if pdb_id not in chains:
                chains[pdb_id] = [self.emaps[pdb_id].chains[0]]
            try:
                process(self.emaps[pdb_id],
                        chains=chains[pdb_id],
                        eta_moieties=eta_moieties[pdb_id],
                        include_residues=include_residues,
                        sdef = sdef,
                        **kwargs) 
                print("Finished:"+str(pdb_id))
            except Exception as e:
                print(e)
                remove_pdbs.append(pdb_id)
                warnings.warn("Could not generate graph for: "+ pdb_id + ". It will not be included in the analysis.")
        for pdb_id in remove_pdbs:
            self.emaps.pop(pdb_id)
        self._align_sequences(chains)
        self.emap_parameters = kwargs
        self.included_chains = chains
        self.included_eta_moieties = eta_moieties
        self.included_standard_residues = include_residues
        print("Processing took:" + str(time.time() - start_time) + " seconds.")

    def _apply_num_labels(self):
        for emap in self.emaps.values():
            protein_graph = emap.init_graph
            for node in protein_graph.nodes:
                protein_graph.nodes[node]['num_label'] = get_numerical_node_label(node, self.node_labels)
            for edge in protein_graph.edges:
                protein_graph.edges[edge]['num_label'] = get_edge_label(protein_graph, edge, self.edge_thresholds)

    def report_header(self):
        full_str = ""
        full_str += "Generated:\n" + str(datetime.datetime.now()) + "\n"
        full_str += "Parameters:\n"
        if not self.emap_parameters:
            full_str += "Custom.\n"
        else:
            full_str += str(self.emap_parameters)
            full_str += "\n"
        full_str += "Chains:\n"
        if not self.included_chains:
            full_str += "Custom.\n"
        else:
            full_str += str(self.included_chains)
            full_str += "\n"
        full_str += "Included non protein moieties:\n"
        if not self.included_eta_moieties:
            full_str += "Custom.\n"
        else:
            full_str += str(self.included_eta_moieties)
            full_str += "\n"
        full_str += "Edge thresholds:\n" + str(self.edge_thresholds) + "\n"
        full_str += "Node labels:\n" + str(self.node_labels) + "\n"
        full_str += "Residue categories:\n" + str(self.residue_categories) + "\n"
        return full_str

    def general_report(self, dest=None):
        ''' Generates general report of all subgraph patterns found in the analysis.

        Returns
        -------
        report: str
            General report of all subgraph patterns found in the analysis.
        '''
        full_str = "Overview of all subgraphs:\n"
        full_str += self.report_header()
        full_str += "\nSubgraphs found:\n\n"
        for fsg in self.subgraph_patterns:
            full_str += self.subgraph_patterns[fsg].general_report() + "\n"
        if dest:
            fi = open(dest, "w")
            fi.write(full_str)
            fi.close()
        return full_str

    def subgraph_report(self, sg_id, dest=None):
        ''' Generates detailed report for a given subgraph pattern.

        Parameters
        -----------
        sg_id: str
            ID corresponding to a :class:`~pyemap.common_paths.FrequentSubgraph` object 
        dest: str, optional
            Destination to write report to file

        Returns
        -------
        report: str
            Detailed report for a particular subgraph pattern.
        '''
        sg = self.subgraph_patterns[sg_id]
        full_str = "Full report for subgraph:" + str(sg_id) + "\n"
        full_str += self.report_header() + "\n"
        full_str += sg.full_report()
        if dest:
            fi = open(dest, "w")
            fi.write(full_str)
            fi.close()
        return full_str

    def _set_node_labels(self,nodes):
        '''
        assert (categories is None and labels is None) or  (categories is not None and labels is not None)
        if categories is not None:
            if "X" in labels or "X" in categories.values():
                raise KeyError("X is reserved for unspecified residue type. Do not use X as a key.")
            if not nodes==['H','W','F',"Y"]:
                warnings.warn("Warning: 'nodes' keyword is incompatible and will be ignored.")
            self.residue_categories = categories
            self.node_labels = labels
        '''
        if nodes is None:
            nodes = [res_name_to_char[x] for x in self.included_standard_residues]
        assert all(x in char_to_res_name for x in nodes)
        num_label = 2
        for res in nodes:
            self.node_labels[res] = num_label
            self.residue_categories[num_label] = res
            num_label+=1
        num_label = max(self.node_labels.values()) + 1
        self.residue_categories[num_label] = "X"
        self.node_labels["X"] = num_label
        self.residue_categories[num_label+1] = "#"
        self.node_labels["#"] = num_label+1

    def add_emap(self, emap_obj):
        ''' Adds a parsed :class:`~pyemap.emap` object to the PDB group.

        Parameters
        ----------
            emap_obj: :class:`~pyemap.emap` object 
                Parsed PDB generated by :func:`~pyemap.parser.parse` or :func:`~pyemap.parser.fetch_and_parse`
        '''
        if emap_obj.pdb_id not in self.emaps:
            self.emaps[emap_obj.pdb_id] = emap_obj
        else:
            print("An emap object with PDB ID:" + str(emap_obj.pdb_id) + " is already in the data set. Skipping...")

    def generate_graph_database(self,node_categories=None,edge_thresholds=[10,15]):
        ''' Generates graph database for analysis by GSpan using specified node labels, node categories, and edge thresholds.

        Parameters
        ----------
        nodes: list of str
            List of one character amino acid codes to be given their own category. The remaining
            AA will be labeled as "X". Default is all standard amino acids receive their own category.
        edge_thresholds: list of float
            List of edge thresholds. Edges with weight below the first value will be given the label 2, edges
            between the 1st and second values will be labeled as 3, and so on. 
            Default is [5,9,13...] for closest atom, and [10,14,18...] for center of mass distance.
        Examples
        ---------
        pg.generate_graph_database(nodes=['W','Y'],edge_thresholds=[5,15])
        '''
        # check if we need to regenerate database
        self._clean_graph_database()
        assert (float(x) for x in edge_thresholds)
        assert all(edge_thresholds[i] <= edge_thresholds[i+1] for i in range(len(edge_thresholds)-1))
        self.edge_thresholds = edge_thresholds.copy()
        self._set_node_labels(node_categories)
        f = open(os.path.join(self.temp_dir, 'graphdatabase.txt'), "w")
        for i, key in enumerate(self.emaps):
            G = self.emaps[key].init_graph
            f.write("t # " + str(i) + "\n")
            for i, node in enumerate(G.nodes):
                f.write("v " + str(i) + " " + str(get_numerical_node_label(node,self.node_labels)) + "\n")
            for i, edge in enumerate(G.edges):
                f.write("e " + str(list(G.nodes()).index(edge[0])) + " " + str(list(G.nodes()).index(edge[1])) + " " +
                        str(get_edge_label(G, edge,self.edge_thresholds)) + "\n")
        f.write("t # -1")
        f.close()
        self._apply_num_labels()

    def run_gspan(self, support, lower_bound=4):
        ''' Runs gSpan algorithm to mine for subgraph patterns, and then identifies each occurence of each subgraph pattern in each PDB which supports it.
        
        Parameters
        ----------
        support: int, optional
            Minimum support number of subgraphs in the search space 
        lower_bound: int, optional
            Minimum number of nodes for subgraphs in the search space
        '''
        start_time = time.time()
        self._clean_subgraphs()
        self.gspan_parameters["support"] = support
        self.gspan_parameters["lower_bound"] = lower_bound
        self.gspan_parameters["graph_specification"] = ""
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
        print("GSpan took:" + str(time.time() - start_time) + " seconds.")
        self._generate_subgraph_patterns()

    def _generate_subgraph_patterns(self):
        buff = open(os.path.join(self.temp_dir, 'gspan_results.out'), "r")
        subgraphs = []
        lines = buff.readlines()
        line_idx = 0
        while line_idx < len(lines):
            line = lines[line_idx]
            if len(line.split()) == 3 and line.split()[0] == "t" and line.split()[1] == "#":
                graph_number = int(line.split()[2])
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
                        G.nodes[node_idx]['label'] = self.residue_categories[node_label]
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
                        support = {}
                        for idx in pdb_list_by_index:
                            support[pdb_list[idx]] = self.emaps[pdb_list[idx]]
                        subgraphs.append(FrequentSubgraph(G,graph_number,support,self.node_labels,self.edge_thresholds))
                    line_idx += 1
            line_idx += 1
        buff.close()
        subgraphs.sort(key=lambda x: x.support_number, reverse=True)
        for sg in subgraphs:
            self.subgraph_patterns[sg.id] = sg

    def find_subgraph(self, graph_specification):
        self._clean_subgraphs()
        self.gspan_parameters["support"] = None
        self.gspan_parameters["lower_bound"] = None
        self.gspan_parameters["graph_specification"] = graph_specification
        node_combs, edge_combs = nodes_and_edges_from_string(graph_specification, self.edge_thresholds, list(self.residue_categories.values()))
        subgraph_patterns = []
        for node_list in node_combs:
            G = nx.Graph()
            for node_idx, node in enumerate(node_list):
                G.add_node(node_idx)
                G.nodes[node_idx]['label'] = node
                G.nodes[node_idx]['num_label'] = self.node_labels[node]
                if node_idx > 0:
                    G.add_edge(node_idx - 1, node_idx)
            for edge_comb in edge_combs:
                for j, edge in enumerate(G.edges):
                    G.edges[edge]['num_label'] = edge_comb[j]
                    G.edges[edge]['label'] = edge_comb[j]
                protein_subgraphs = []
                support = {}
                for pdb_id in self.emaps:
                    GM = get_graph_matcher(self.emaps[pdb_id].init_graph, G)
                    if GM.subgraph_is_monomorphic():
                        support[pdb_id] = self.emaps[pdb_id]
                if len(support) > 0:
                    subgraph_patterns.append(FrequentSubgraph(G,len(subgraph_patterns),support,self.node_labels,self.edge_thresholds))
        subgraph_patterns.sort(key=lambda x: x.support_number, reverse=True)
        for fs in subgraph_patterns:
            self.subgraph_patterns[fs.id] = fs