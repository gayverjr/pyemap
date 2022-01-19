import numpy as np
import os
import shutil
from ..process_data import process
import networkx as nx
import datetime
from ..data import char_to_res_name
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from .frequent_subgraph import SubgraphPattern
import warnings
from .utils import get_edge_label, get_numerical_node_label, get_graph_matcher, extract_chain, set_defaults
from ..pyemap_exceptions import *
from gspan_mining import gSpan
from tempfile import NamedTemporaryFile
import sys
from io import StringIO


def moieties_on_chains(chains, moieties):
    ''' Extract eta moieties which belong to chains

    Parameters
    -----------
    chains: list of str
        Chains included in the analysis
    moieties: list of str
        Non-protein moieties included in the analysis

    '''
    remove = False
    for itm in moieties:
        if extract_chain(itm) not in chains:
            remove = True
            warnings.warn(itm + " is not part of an included chain, and will not be included in the analysis.")
    if remove:
        return [item for item in moieties if extract_chain(item) in chains]
    else:
        return moieties


def nodes_and_edges_from_string(graph_str, edge_thresholds, residue_categories):
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
    graph_str = graph_str.replace(" ", "")
    graph_str = list(graph_str)
    node_list = []
    idx = 0
    while idx < len(graph_str):
        node_list.append(str(graph_str[idx]))
        idx += 1
    l1 = []
    for i in range(0, len(edge_thresholds)+1):
        l1.append(i + 1)
    from itertools import product
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
    edge_combs = list(product(l1, repeat=len(node_list) - 1))
    if len(edge_combs) == 0:
        edge_combs = [tuple([1 for x in range(0,len(graph_str)-1)])]
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
    subgraph_patterns: dict of str: :class:`~pyemap.graph_mining.SubgraphPattern`
        Dict of subgraph patterns found by GSpan. Keys are the unique IDs of the :class:`~pyemap.graph_mining.SubgraphPattern` objects.
    fasta: str
        List of unaligned sequences in FASTA format
    aligned_fasta: str
        List of sequences in FASTA format after multiple sequence alignment

    '''
    def __init__(self, title):
        ''' Initializes PDBGroup object
        
        Parameters
        ----------
        title: str
            Title of PDB group 

        '''
        self.title = title
        self.emaps = {}
        self.subgraph_patterns = {}
        self.fasta=""
        self.aligned_fasta=""
        self._emaps = {}
        self._node_labels = {}
        self._residue_categories = {}
        self._edge_thresholds = []
        self._substitutions = []
        self._emap_parameters = {}
        self._gspan_parameters = {}
        self._included_eta_moieties = {}
        self._included_chains = {}
        self._include_residues = []
        self._sequences = {}
        self._aligned_sequences = {}
        self._graph_database = ""

    def _clean_graph_database(self):
        ''' Cleans graph database and associated parameters.

        '''
        self._node_labels = {}
        self._residue_categories = {}
        self._edge_thresholds = []
        self._substitutions = []
        self._graph_database = ""
        self._clean_subgraphs()

    def _clean_subgraphs(self):
        ''' Cleans results of subgraph minings.

        '''
        self._gspan_parameters = {}
        self.subgraph_patterns = {}

    def _reset_process(self):
        ''' Cleans results of processing :class:`~pyemap.emap` objects.

        '''
        self._emap_parameters = {}
        self._included_eta_moieties = {}
        self._included_chains = {}
        self._include_residues = []
        self._sequences = {}
        self._aligned_sequences = {}
        self.fasta=""
        self.aligned_fasta=""
        self._clean_graph_database()

    def _align_sequences(self):
        ''' Performs multiple sequence alignment using MUSCLE.

        Should be performed after all :class:`~pyemap.emap` objects have been processed.
        :class:`Bio.PDB.Residue.Residue` objects are updated with an additional attribute 
        `aligned_residue_number` which contains the residue number w.r.t. the sequence alignment. 
        Non-protein residues will simply have the original residue number placed in this attribute for consistency.

        If the alignment fails, the original residue number is copied into the `aligned_residue_number` attribute for consistency.
        

        Notes
        ------
        Successful alignment requires installing MUSCLE (https://www.ebi.ac.uk/Tools/msa/muscle/) and having the exectuable in the path.

        '''
        inp = NamedTemporaryFile(mode='w', delete=False)
        out = NamedTemporaryFile(mode='w+', delete=False)
        orig_fasta = ""
        try:
            for pdb_id, emap in self.emaps.items():
                for chain in emap.active_chains:
                    inp.write(emap.sequences[chain] + "\n")
                    orig_fasta+=emap.sequences[chain] + "\n"
            inp.close()
            try:
                muscle_cline = MuscleCommandline(input=inp.name, out=out.name)
                muscle_cline()
                seqIO = SeqIO.parse(out, "fasta")
                aligned_fasta=""
                for record in seqIO:
                    self._aligned_sequences[record.id] = record.seq
                    aligned_fasta+='>'+str(record.id)+'\n'
                    aligned_fasta+=str(record.seq)+'\n'
                self.fasta = orig_fasta
                self.aligned_fasta = aligned_fasta
                # now lets save the updated sequence numbers
                for pdb_id, emap in self.emaps.items():
                    for chain, residues in emap.active_chains.items():
                        original_idx = 1
                        aligned_idx = 1
                        seq_map = {}
                        if pdb_id + ":" + chain in self._aligned_sequences:
                            aligned_seq = self._aligned_sequences[pdb_id + ":" + chain]
                            for res in aligned_seq:
                                if not res == "-":
                                    seq_map[original_idx] = aligned_idx
                                    original_idx += 1
                                aligned_idx += 1
                        for residue in residues:
                            try:
                                resnum = residue.sequence_index
                                residue.aligned_residue_number = seq_map[int(resnum)]
                            except:
                                residue.aligned_residue_number = 'X'
            except Exception:
                for pdb_id, emap in self.emaps.items():
                    for chain, residues in emap.active_chains.items():
                        for residue in residues:
                            residue.aligned_residue_number = residue.id[1]
                warnings.warn("Warning: could not align sequences. Make sure that MUSCLE (https://www.drive5.com/muscle/manual/) is installed and "\
                "accessible in the current path. Original residue numbers will be used for sequence similarity.")
        finally:
            os.remove(inp.name)
            os.remove(out.name)

    def _setup_process(self, chains={}, eta_moieties={}, include_residues=["Y", "W"], **kwargs):
        ''' Sets up parameters to begin emap graph generation.

        Parameters
        -----------
        chains: dict of str: list of str, optional
            Chains to include for each PDB
        eta_moieties: dict of str: list of str, optional
            Dict containing list of ETA moieties(specified by their residue label to include for each PDB
        include_residues: list of str, optional
            List of 1-letter standard AA codes to include in the graph
        **kwargs
            For a list of accepted kwargs, see the documentation for :func:`~pyemap.process_data.process`.

        This function should be used to set up :func:`_process_emap`.   
          
        '''
        self._reset_process()
        for pdb_id in self.emaps:
            if pdb_id not in chains:
                chains[pdb_id] = [self.emaps[pdb_id].chains[0]]
            if pdb_id not in eta_moieties:
                eta_moieties[pdb_id] = [
                    item for item in self.emaps[pdb_id].eta_moieties.keys() if extract_chain(item) in chains
                ]
            eta_moieties[pdb_id] = moieties_on_chains(chains[pdb_id], eta_moieties[pdb_id])
        self._emap_parameters = kwargs
        self._included_chains = chains
        self._included_eta_moieties = eta_moieties
        self._include_residues = include_residues


    def _process_emap(self, pdb_id):
        '''Processes :class:`~pyemap.emap` object in order to generate protein graph. 
        
        :func:`_setup_process` should be executed first if using this function.

        Parameters
        ----------
        pdb_id: str
            PDB ID corresponding to :class:`~pyemap.emap`

        '''
        process(self._emaps[pdb_id],
                chains=self._included_chains[pdb_id],
                eta_moieties=self._included_eta_moieties[pdb_id],
                include_residues=self._include_residues,
                **self._emap_parameters)
        print("Finished:" + str(pdb_id))

    def process_emaps(self, chains={}, eta_moieties={}, include_residues=['Y', 'W'], **kwargs):
        ''' Processes :class:`~pyemap.emap` objects in order to generate protein graphs. 
        
        Should be executed once all of the :class:`~pyemap.emap` objects have been added to the group.

        Parameters
        -----------
        chains: dict of str: list of str, optional
            Chains to include for each PDB. The special keyword 'All' is also accepted.
        eta_moieties: dict of str: list of str, optional
            Dict containing list of ETA moieties(specified by their residue label) to include for each PDB
        include_residues: list of str, optional
            List of 1-letter standard AA codes to include in the graph
        **kwargs
            For a list of accepted kwargs, see the documentation for :func:`~pyemap.process_data.process`.

        Examples
        --------
        >>> my_pg = pyemap.graph_mining.PDBGroup()
        >>> # Add pdbs 1u3d,1u3c,6PU0,4I6G,2J4D ...
        >>> eta_moieties = {'1u3d': ['FAD510(A)-2'], '1u3c': ['FAD510(A)-2'], '6PU0': ['FAD501(A)-2'], '4I6G': ['FAD900(A)-2'], '2J4D': ['FAD1498(A)-2']}
        >>> chains = {'1u3d': ['A'], '1u3c': ['A'], '6PU0': ['A'], '4I6G': ['A'], '2J4D': ['A']}
        >>> my_pg.process_emaps(chains=chains,eta_moieties=eta_moieties)

        '''
        self._reset_process()
        remove_pdbs = []
        kwargs = set_defaults(kwargs)
        try:
            if chains.upper()=='ALL':
                chains = {}
                for pdb_id in self._emaps:
                    chains[pdb_id] = self._emaps[pdb_id].chains
        except:
            pass
        for pdb_id in self._emaps:
            if pdb_id not in chains:
                chains[pdb_id] = [self._emaps[pdb_id].chains[0]]
            if pdb_id not in eta_moieties:
                eta_moieties[pdb_id] = [
                    item for item in self._emaps[pdb_id].eta_moieties.keys() if extract_chain(item) in chains[pdb_id]
                ]
            try:
                eta_moieties[pdb_id] = moieties_on_chains(chains[pdb_id], eta_moieties[pdb_id])
                process(self._emaps[pdb_id],
                        chains=chains[pdb_id],
                        eta_moieties=eta_moieties[pdb_id],
                        include_residues=include_residues,
                        **kwargs)
                print("Finished:" + str(pdb_id))
            except Exception as e:
                print(e)
                remove_pdbs.append(pdb_id)
                warnings.warn("Could not generate graph for: " + pdb_id + ". It will not be included in the analysis.")
        for pdb_id in remove_pdbs:
            self.emaps.pop(pdb_id)
        if len(self.emaps) < 2:
            raise PyeMapMiningException("Not enough graphs could be generated for mining.")
        self._align_sequences()
        self._emap_parameters = kwargs
        self._included_chains = chains
        self._included_eta_moieties = eta_moieties
        self._include_residues = include_residues

    def _apply_num_labels(self):
        ''' Adds numerical node labels to the :class:`networkx.Graph` objects as node attributes. 

        '''
        for emap in self.emaps.values():
            protein_graph = emap.init_graph
            for node in protein_graph.nodes:
                protein_graph.nodes[node]['num_label'] = get_numerical_node_label(node, self._node_labels)
            for edge in protein_graph.edges:
                protein_graph.edges[edge]['num_label'] = get_edge_label(protein_graph, edge, self._edge_thresholds)                

    def _report_header(self):
        ''' Generates header for reports.

        Includes timestamp and relevant parameters.

        '''
        full_str = ""
        full_str += "Generated:\n" + str(datetime.datetime.now()) + "\n"
        full_str += "Graph Parameters:\n"
        if not self._emap_parameters:
            full_str += "Custom.\n"
        else:
            full_str += str(self._emap_parameters)
            full_str += "\n"
        full_str += "Included residues:\n"
        if not self._include_residues:
            full_str += "Custom.\n"
        else:
            full_str += str(self._include_residues)
            full_str += "\n"
        full_str += "Mining parameters:\n"
        if not self._gspan_parameters:
            full_str += "Custom.\n"
        else:
            full_str += str(self._gspan_parameters)
            full_str += "\n"
        full_str += "Chains:\n"
        if not self._included_chains:
            full_str += "Custom.\n"
        else:
            full_str += str(self._included_chains)
            full_str += "\n"
        full_str += "Included non protein moieties:\n"
        if not self._included_eta_moieties:
            full_str += "Custom.\n"
        else:
            full_str += str(self._included_eta_moieties)
            full_str += "\n"
        full_str += "Edge thresholds:\n" + str(self._edge_thresholds) + "\n"
        full_str += "Node labels:\n" + str(self._node_labels) + "\n"
        full_str += "Residue categories:\n" + str(self._residue_categories) + "\n"
        return full_str

    def mining_report(self, dest=None):
        ''' Generates general report of all subgraph patterns found in the analysis.

        Parameters
        -----------
        dest: str, optional
            Destination to write report to file

        Returns
        -------
        report: str
            General report of all subgraph patterns found in the analysis.

        '''
        full_str = "Overview of all subgraphs:\n"
        full_str += self._report_header()
        full_str += "\nSubgraphs found:\n\n"
        for fsg in self.subgraph_patterns:
            full_str += self.subgraph_patterns[fsg].general_report() + "\n"
        if dest:
            fi = open(dest, "w")
            fi.write(full_str)
            fi.close()
        return full_str

    def _set_node_labels(self):
        ''' Sets node labels for mining. 

        Specified AA residues will be given a unique category, all others will be labeled as 'X'. 
        Non-protein residues will be labeled as '#'.

        Parameters
        ----------
        nodes: list of str
            Residue labels which receive their own category. If None, all residues in `self._include_residues` 
            receive their own category.

        '''
        num_label = 2
        for res in self._include_residues:
            if res not in self._substitutions:
                self._node_labels[res] = num_label
                self._residue_categories[num_label] = res
                num_label += 1
        if len(self._node_labels) > 1:
            num_label = max(self._node_labels.values()) + 1
        self._residue_categories[num_label] = "X"
        self._node_labels["X"] = num_label
        self._residue_categories[num_label + 1] = "#"
        self._node_labels["#"] = num_label + 1

    def add_emap(self, emap_obj):
        ''' Adds an :class:`~pyemap.emap` object to the PDB group.

        Parameters
        ----------
        emap_obj: :class:`~pyemap.emap` object 
            Parsed PDB generated by :func:`~pyemap.parser.parse` or :func:`~pyemap.parser.fetch_and_parse`

        '''
        if emap_obj.pdb_id not in self._emaps:
            self.emaps[emap_obj.pdb_id] = emap_obj
            self._emaps[emap_obj.pdb_id] =  emap_obj
        else:
            print("An emap object with PDB ID:" + str(emap_obj.pdb_id) + " is already in the data set. Skipping...")

    def generate_graph_database(self, sub=[], edge_thresh=[]):
        ''' Generates graph database for mining.

        Parameters
        ----------
        sub: list of str, optional
            List of 1-character amino acid codes to be labeled as "X". All other included standard amino acids receive their own category.
        edge_thresholds: list of float, optional
            List of edge thresholds. Edges with weight below the first value will be given the label 2, edges
            between the 1st and second values will be labeled as 3, and so on. 

        Examples
        ---------
        >>> pg.generate_graph_database(['W','Y'],[12,15])

        '''
        self._clean_graph_database()
        try:
            for i in range(0,len(edge_thresh)-1):
                assert float(edge_thresh[i]) <= float(edge_thresh[i + 1])
            self._edge_thresholds = edge_thresh.copy()
        except Exception as e:
            raise PyeMapGraphDatabaseException("Invalid specification of edge thresholds.") from e
        try:
            for x in sub:
              assert x in char_to_res_name  
              assert x.upper() in self._include_residues
            self._substitutions = sub.copy()
            self._set_node_labels()
        except Exception as e:
            raise PyeMapGraphDatabaseException("Invalid specification of substitutions.") from e
        f = StringIO("")
        for i, key in enumerate(self.emaps):
            G = self.emaps[key].init_graph
            f.write("t # " + str(i) + "\n")
            for i, node in enumerate(G.nodes):
                f.write("v " + str(i) + " " + str(get_numerical_node_label(node, self._node_labels)) + "\n")
            for i, edge in enumerate(G.edges):
                f.write("e " + str(list(G.nodes()).index(edge[0])) + " " + str(list(G.nodes()).index(edge[1])) + " " +
                        str(get_edge_label(G, edge, self._edge_thresholds)+1) + "\n")
        f.write("t # -1")
        self._graph_database = f.getvalue()
        f.close()
        self._apply_num_labels()

    def run_gspan(self, min_support, min_num_vertices=4, max_num_vertices=float('inf'), **kwargs):
        ''' Mines for common subgraphs using gSpan algorithm. Results are stored as :class:`~pyemap.graph_mining.SubgraphPattern` objects 
        in the `subgraph_patterns` dictionary.

        References
        ----------
        Yan, Xifeng, and Jiawei Han. "gspan: Graph-based substructure pattern mining." 2002 IEEE International Conference on Data Mining, 
        2002. Proceedings.. IEEE, 2002.

        Parameters
        ----------
        min_support: int
            Minimum support number of subgraphs in the search space 
        min_num_vertices: int, optional
            Minimum number of nodes for subgraphs in the search space
        max_num_vertices: int, optional
            Maximum number of nodes for subgraphs in the search space
        **kwargs
            See https://github.com/betterenvi/gSpan for a list of accepted kwargs.

        '''
        self._clean_subgraphs()
        self._gspan_parameters["min_support"] = min_support
        self._gspan_parameters['min_num_vertices'] = min_num_vertices
        self._gspan_parameters['max_num_vertices'] = max_num_vertices
        self._gspan_parameters["graph_specification"] = []
        db = NamedTemporaryFile(mode='w', delete=False)
        print(self._graph_database, file=db)
        db.close()
        old_stdout = sys.stdout
        sys.stdout = mystdout = StringIO()
        gs = gSpan(database_file_name=db.name,
                   min_support=min_support,
                   min_num_vertices=min_num_vertices,
                   max_num_vertices=max_num_vertices,
                   where=True,
                   **kwargs)
        gs.run()
        sys.stdout = old_stdout
        self._gspan_results = mystdout.getvalue()
        self._generate_subgraph_patterns()

    def _generate_subgraph_patterns(self):
        ''' Parses gSpan output and generates :class:`~pyemap.graph_mining.SubgraphPattern` objects. Results are stored in
        `self.subgraph_patterns`.

        '''
        with StringIO(self._gspan_results) as buff:
            subgraphs = []
            lines = buff.readlines()
            line_idx = 0
            while line_idx < len(lines):
                line = lines[line_idx]
                if len(line.split()) == 3 and line.split()[0] == "t" and line.split()[1] == "#":
                    graph_number = int(line.split()[2])
                    line_idx += 1
                    G = nx.Graph()
                    line = lines[line_idx]
                    while "---" not in line:
                        line = lines[line_idx]
                        if len(line.split()) > 1 and line.split()[0] == "v":
                            node_idx = int(line.split()[1])
                            node_label = int(line.split()[2])
                            G.add_node(node_idx)
                            G.nodes[node_idx]['label'] = self._residue_categories[node_label]
                            G.nodes[node_idx]['num_label'] = node_label
                        if len(line.split()) > 1 and line.split()[0] == "e":
                            idx1 = int(line.split()[1])
                            idx2 = int(line.split()[2])
                            edge_label = int(line.split()[3])-1
                            G.add_edge(idx1, idx2, label=edge_label)
                            G.edges[(idx1, idx2)]['num_label'] = edge_label
                        if "where" in line:
                            pdb_list_by_index = line[7:-2].strip('][').split(', ')
                            pdb_list_by_index = list(np.array(pdb_list_by_index, dtype=int))
                            pdb_list = list(self.emaps.keys())
                            support = {}
                            for idx in pdb_list_by_index:
                                support[pdb_list[idx]] = self.emaps[pdb_list[idx]]
                            subgraphs.append(
                                SubgraphPattern(G, graph_number, support, self._node_labels, self._edge_thresholds))
                        line_idx += 1
                line_idx += 1
            buff.close()
        subgraphs.sort(reverse=True)
        for graph_number,sg in enumerate(subgraphs):
            sg._update_id(graph_number)
            self.subgraph_patterns[sg.id] = sg

    def find_subgraph(self, graph_specification):
        ''' Finds a specified subgraph by searching for monomorphisms in each protein graph.

        Currently, only linear chains are supported. If edge thresholds are used, all possible combinations of 
        edges will be searched for.

        Parameters
        ----------
        G: str or :class:`networkx.Graph`, optional
            String specifying linear chain to search for, or a graph object to search for

        Notes
        ------
        Special characters for graph_specification:
        * - wildcard character
        # - non-protein residue
        
        Only linear chains are supported when G is a string. Any valid :class:`networkx.Graph` can be used, however, care 
        should be taken to set node and edge labels appropriately, as the analysis relies on the 'num_label' attribute for both 
        nodes and edges.

        Examples
        ---------
        >>> my_pg.find_subgraph('WWW#')

        '''
        self._clean_subgraphs()
        self._gspan_parameters["support"] = None
        self._gspan_parameters["min_num_vertices"] = None
        self._gspan_parameters["max_num_vertices"] = None
        self._gspan_parameters["graph_specification"] = graph_specification
        try:
            node_combs, edge_combs = nodes_and_edges_from_string(graph_specification, self._edge_thresholds, list(self._residue_categories.values()))
        except:
            raise PyeMapMiningException("Could not parse graph from string.")
        try:
            frequent_subgraphs = []
            for node_list in node_combs:
                G = nx.Graph()
                for node_idx, node in enumerate(node_list):
                    G.add_node(node_idx)
                    G.nodes[node_idx]['label'] = node
                    G.nodes[node_idx]['num_label'] = self._node_labels[node]
                    if node_idx > 0:
                        G.add_edge(node_idx - 1, node_idx)
                for edge_comb in edge_combs:
                    for j, edge in enumerate(G.edges):
                        G.edges[edge]['num_label'] = edge_comb[j]
                        G.edges[edge]['label'] = edge_comb[j]
                    support = {}
                    for pdb_id in self.emaps:
                        GM = get_graph_matcher(self.emaps[pdb_id].init_graph, G)
                        if GM.subgraph_is_monomorphic():
                            support[pdb_id] = self.emaps[pdb_id]
                    if len(support) > 0:
                        frequent_subgraphs.append(SubgraphPattern(G,len(frequent_subgraphs),support,self._node_labels,self._edge_thresholds))
            frequent_subgraphs.sort(reverse=True)
            for graph_number,fs in enumerate(frequent_subgraphs):
                fs._update_id(graph_number)
                self.subgraph_patterns[fs.id] = fs
        except Exception as e:
            raise PyeMapMiningException("Could not generate graphs using the specified string.") from e