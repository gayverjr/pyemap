import os
import sys
import unittest
import warnings
import pyemap
from math import isclose

class SingleChainParams(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1u3d.pdb")) 

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
        assert "PHE" not in resnames
        assert "HIS" not in resnames
        #phenylalanine on
        pyemap.process(self.my_emap,phe=True)
        resnames=[]
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert "TYR" in resnames
        assert "TRP" in resnames
        assert "PHE" in resnames
        assert "HIS" not in resnames
        #histidine on
        resnames=[]
        pyemap.process(self.my_emap,his=True)
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert "TYR" in resnames
        assert "TRP" in resnames
        assert "PHE" not in resnames
        assert "HIS" in resnames
        #tyrosine off
        resnames=[]
        pyemap.process(self.my_emap,tyr=False)
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert "TYR" not in resnames
        assert "TRP" in resnames
        assert "PHE" not in resnames
        assert "HIS" not in resnames
        #tryptophan off
        resnames=[]
        pyemap.process(self.my_emap,trp=False)
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert "TYR" in resnames
        assert "TRP" not in resnames
        assert "PHE" not in resnames
        assert "HIS" not in resnames


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
        pyemap.process(self.my_emap,eta_moieties=[], custom = "(134-140),(109-117)")
        resnames=[]
        for residue in self.my_emap.residues.values():
            resnames.append(residue.resname)
        assert all(elem in resnames for elem in ["CUST-1","CUST-2"])
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
        cls.my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/2oal.pdb")) 

    def test_chains(self):
        #default all chains on
        pyemap.process(self.my_emap)
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
