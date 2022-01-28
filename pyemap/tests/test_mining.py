## Copyright (c) 2017-2022, James Gayvert, Ruslan Tazhigulov, Ksenia Bravaya
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import unittest
import pyemap
from pyemap.graph_mining import PDBGroup
import os
import numpy as np

class PDBGroupProcess(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pdb_ids = ['1IQR', '1NP7', '1U3C', '1U3D', '2J4D', '3FY4', '3ZXS', '4GU5', '4I6G', '4U63', '6FN2', '6KII', '6LZ3', '6PU0'] 
        cls.pg = PDBGroup('My Group') 
        for pdb in cls.pdb_ids: 
            cls.pg.add_emap(pyemap.fetch_and_parse(pdb)) 
        cls.pg.add_emap(pyemap.fetch_and_parse('1IQR')) 
        assert len(cls.pg._emaps) == len(cls.pdb_ids)
    
    @classmethod
    def tearDownClass(cls):
        for pdb in cls.pdb_ids:
            os.remove(pdb+".pdb")

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
        eta_moieties = {}
        for pdb_id in self.pg.emaps:
            eta_moieties[pdb_id] = []
        eta_moieties['1U3D'] = ['FAD510(A)-2']
        self.pg.process_emaps(edge_prune='DEGREE',eta_moieties=eta_moieties)
        self.pg.generate_graph_database(edge_thresh=[12,14])
        self.pg.find_subgraph('WWW#')
        sg = self.pg.subgraph_patterns['1_WWW#_1']
        sg.find_protein_subgraphs()
        assert len(sg.groups)==1

    def test_no_subgraphs(self):
        self.pg.process_emaps(edge_prune='DEGREE')
        self.pg.generate_graph_database()
        self.pg.run_gspan(20)
        assert len(self.pg.subgraph_patterns) == 0
        self.pg.find_subgraph('WWWX')
        assert len(self.pg.subgraph_patterns) == 0

    def test_default_chains_and_moieties(self):
        self.pg.process_emaps(edge_prune='DEGREE')
        total_length_chains = np.sum([len(x) for x in self.pg._included_chains.values()])
        assert total_length_chains == 14
        total_length_moieties = np.sum([len(x) for x in self.pg._included_eta_moieties.values()])
        assert total_length_moieties == 40

        self.pg.process_emaps(edge_prune='DEGREE',chains='ALL')
        total_length_chains = np.sum([len(x) for x in self.pg._included_chains.values()])
        assert total_length_chains == 25
        total_length_moieties = np.sum([len(x) for x in self.pg._included_eta_moieties.values()])
        assert total_length_moieties == 75

        self.pg.process_emaps(edge_prune='DEGREE',eta_moieties={'2J4D':['MHF1502(B)-1']})
        total_length_chains = np.sum([len(x) for x in self.pg._included_chains.values()])
        assert total_length_chains == 14
        total_length_moieties = np.sum([len(x) for x in self.pg._included_eta_moieties.values()])
        assert total_length_moieties == 36
         
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

    def test_reports(self):
        self.pg.process_emaps(edge_prune='DEGREE')
        self.pg.generate_graph_database()
        self.pg.find_subgraph('WWW#')
        sg = self.pg.subgraph_patterns['1_WWW#_14']
        sg.find_protein_subgraphs()
        assert self.pg.mining_report() != None
        assert sg.general_report() != None
        assert sg.full_report() != None

    def test_visualize(self):
        self.pg.process_emaps(edge_prune='DEGREE')
        self.pg.generate_graph_database()
        self.pg.find_subgraph('WWW#')
        sg = self.pg.subgraph_patterns['1_WWW#_14']
        sg.find_protein_subgraphs()
        sg.subgraph_to_Image()
        sg.subgraph_to_file(dest='test.png')
        os.remove('test.png')
        sg.subgraph_to_Image(id=next(iter(sg.protein_subgraphs)))
        sg.subgraph_to_file(id=next(iter(sg.protein_subgraphs)),dest='test.png')
        os.remove('test.png')