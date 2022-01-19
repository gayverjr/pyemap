import unittest
import pyemap
from pyemap.graph_mining import PDBGroup
import os
from math import isclose

class PDBGroupProcess(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pdb_ids = ['1IQR', '1NP7', '1U3C', '1U3D', '2J4D', '3FY4', '3ZXS', '4GU5', '4I6G', '4U63', '6FN2', '6KII', '6LZ3', '6PU0'] 
        cls.pg = PDBGroup('My Group') 
        for pdb in cls.pdb_ids: 
            cls.pg.add_emap(pyemap.fetch_and_parse(pdb)) 
    
    @classmethod
    def tearDownClass(cls):
        for pdb in cls.pdb_ids:
            os.remove(pdb+".pdb")
         
    def test_gspan(self):
        self.pg.process_emaps(edge_prune='DEGREE')
        self.pg.generate_graph_database(edge_thresh=[12]) 
        self.pg.run_gspan(14,min_num_vertices=4,max_num_vertices=5)
        assert len(self.pg.subgraph_patterns)==19
        self.pg.run_gspan(14,min_num_vertices=3,max_num_vertices=3)
        assert len(self.pg.subgraph_patterns)==15
        self.pg.run_gspan(14,4,4) 
        assert len(self.pg.subgraph_patterns)==12

    def test_graph_db(self):
        self.pg.process_emaps(edge_prune='DEGREE')
        self.pg.generate_graph_database()
        self.pg.find_subgraph('WWW*')
        assert len(self.pg.subgraph_patterns) == 3
        self.pg.process_emaps(edge_prune='DEGREE',include_residues=['Y','W','F'])
        self.pg.generate_graph_database()
        self.pg.find_subgraph('WWW*')
        assert len(self.pg.subgraph_patterns) == 4
        self.pg.process_emaps(edge_prune='DEGREE')
        self.pg.generate_graph_database(edge_thresh=[12])
        self.pg.find_subgraph('WWW*')
        assert len(self.pg.subgraph_patterns) == 23   
        self.pg.process_emaps(edge_prune='DEGREE',include_residues=['Y','W','F'])
        self.pg.generate_graph_database(sub=['Y','F'])
        self.pg.find_subgraph('WWW*')
        assert len(self.pg.subgraph_patterns) == 3

    def test_clustering(self):
        self.pg.process_emaps(edge_prune='DEGREE')
        self.pg.generate_graph_database()
        self.pg.find_subgraph('WWW#')
        sg = self.pg.subgraph_patterns['1_WWW#_14']
        sg.find_protein_subgraphs(clustering_option="structural")
        assert len(sg.protein_subgraphs)==137
        g1 = sg.groups[1]
        assert len(g1) == 12
        sg.set_clustering("sequence")
        g1 = sg.groups[1]
        assert len(g1) == 14

    def test_reports(self):
        self.pg.process_emaps(edge_prune='DEGREE')
        self.pg.generate_graph_database()
        self.pg.find_subgraph('WWW#')
        sg = self.pg.subgraph_patterns['1_WWW#_14']
        sg.find_protein_subgraphs()
        assert self.pg.mining_report() != None
        assert sg.general_report() != None
        assert sg.full_report() != None

    







