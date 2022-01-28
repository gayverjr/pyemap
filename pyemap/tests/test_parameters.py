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
import os
from math import isclose

class SingleChainParams(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.my_emap = pyemap.fetch_and_parse("1u3d") 
    
    @classmethod
    def tearDownClass(cls):
        os.remove("1u3d.pdb")

    def test_eta_moities(self):
        assert len(self.my_emap.eta_moieties) >= 3

    def test_aromatic_amino_acid_options(self):
        pyemap.process(self.my_emap)
        ##default, tryptophan and tyrosine on
        resnames=[]
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert "TYR" in resnames
        assert "TRP" in resnames
        #phenylalanine on
        pyemap.process(self.my_emap,include_residues=["W","Y","F"])
        resnames=[]
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert "TYR" in resnames
        assert "TRP" in resnames
        assert "PHE" in resnames
        #histidine on
        resnames=[]
        pyemap.process(self.my_emap,include_residues=["W","Y","H"])
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert "TYR" in resnames
        assert "TRP" in resnames
        assert "HIS" in resnames
        #tyrosine off
        resnames=[]
        pyemap.process(self.my_emap,include_residues=["W"])
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert "TYR" not in resnames
        assert "TRP" in resnames
        #tryptophan off
        resnames=[]
        pyemap.process(self.my_emap,include_residues=["C","Y","W"])
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert "TYR" in resnames
        assert "TRP" in resnames
        assert "CYS" in resnames


    def test_edge_prune_options(self):
        pyemap.process(self.my_emap,edge_prune="DEGREE")
        G = self.my_emap.init_graph.copy()
        pyemap.process(self.my_emap,edge_prune="PERCENT")
        G2 = self.my_emap.init_graph.copy()
        pyemap.process(self.my_emap,edge_prune="DEGREE",max_degree=2)
        G3 = self.my_emap.init_graph.copy()
        assert G.size()!=G2.size() and G2.size()!=G.size() and G2.size()!=G3.size()
        assert max(G.degree, key=lambda x: x[1])[1] == 4
        assert max(G3.degree, key=lambda x: x[1])[1] == 2


    def test_distance_options(self):   
        pyemap.process(self.my_emap)
        com_weight = self.my_emap.init_graph["W385(A)"]["Y53(A)"]['weight']
        assert isclose(com_weight,6.18,abs_tol=1e-2)
        pyemap.process(self.my_emap,dist_def=1)
        closest_atom_weight = self.my_emap.init_graph["W385(A)"]["Y53(A)"]['weight']
        assert isclose(closest_atom_weight,3.97,abs_tol=1e-2)

    def test_choose_eta_moieties(self):        
        pyemap.process(self.my_emap)
        eta_moieties_resnames = []
        for res in self.my_emap.eta_moieties:
            eta_moieties_resnames.append(res)
        #default
        resnames=[]
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert all(elem in resnames for elem in eta_moieties_resnames)
        #deselected eta moiety
        pyemap.process(self.my_emap,eta_moieties=eta_moieties_resnames[:-1])
        resnames=[]
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert all(elem in resnames for elem in eta_moieties_resnames[:-1])
        assert eta_moieties_resnames[-1] not in resnames
    
    def test_custom_atom_range(self):
         #custom residues
        pyemap.process(self.my_emap,eta_moieties=[], custom = "(3960-3969),(3970-3980,3982,3984-3987)")
        resnames=[]
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert all(elem in resnames for elem in ["CUST-1","CUST-2"])
        try:
            self.my_emap.residue_to_file('CUST-1')
            assert False
        except KeyError:
            assert True
        try:
            self.my_emap.residue_to_Image('CUST-1')
            assert False
        except KeyError:
            assert True
        #example for bad custom residues, shoud raise exception
        try:
            pyemap.process(self.my_emap, custom = "(28-41)")
            assert False
        except Exception as e:
            assert True

    def test_graph_parameters(self):
        pyemap.process(self.my_emap,num_st_dev_edges=5,percent_edges=3)
        g1 = self.my_emap.init_graph
        # check st_dev_edges works
        pyemap.process(self.my_emap,num_st_dev_edges=1,percent_edges=3)
        g2 = self.my_emap.init_graph
        assert g1.edges() != g2.edges()
        #check distance cutoff works    
        pyemap.process(self.my_emap,num_st_dev_edges=5,percent_edges=3,distance_cutoff=30)
        g3 = self.my_emap.init_graph
        assert g1.edges() != g3.edges()
        #check percent edges works
        pyemap.process(self.my_emap,num_st_dev_edges=5,percent_edges=1,distance_cutoff=30)
        g4 = self.my_emap.init_graph
        assert g1.edges() != g4.edges()

    def test_penalty_function_parameters(self):  
        #alpha
        pyemap.process(self.my_emap,coef_alpha=0.5)
        alpha_weight = self.my_emap.init_graph["W385(A)"]["Y53(A)"]['weight']
        assert isclose(alpha_weight,6.48,abs_tol=1e-2)
        #Roffset
        pyemap.process(self.my_emap,r_offset=1)
        alpha_weight = self.my_emap.init_graph["W385(A)"]["Y53(A)"]['weight']
        assert isclose(alpha_weight,5.18,abs_tol=1e-2)
        #Roffset
        pyemap.process(self.my_emap,exp_beta=5)
        beta_weight = self.my_emap.init_graph["W385(A)"]["Y53(A)"]['weight']
        assert isclose(beta_weight,13.44,abs_tol=1e-2)


class MultiChainParams(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.my_emap = pyemap.fetch_and_parse("2oal")

    def test_chains(self):
        #default all chains on
        pyemap.process(self.my_emap,chains=["A","B"])
        labels=[]
        for residue in self.my_emap.residues.values():
            labels.append(residue.node_label)
        assert all(elem in labels for elem in ["W18(A)" and 'Y298(B)'])
        #B chain off
        pyemap.process(self.my_emap,chains=["A"])
        labels=[]
        for residue in self.my_emap.residues.values():
            labels.append(residue.node_label)
        assert "W18(A)" in labels
        assert 'Y298(B)' not in labels
