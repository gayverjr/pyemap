from Bio.PDB import PDBIO, FastMMCIFParser, PDBParser
import pyemap
import os
import sys
from math import isclose

def test_parse():
    #fetch from pdb database
    #my_emap = pyemap.fetch_and_parse("1u3d")
    #os.remove("1u3d.pdb")
    #assert len(list(my_emap.structure.get_residues()))==577 
    #from file
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1u3d.pdb"))  
    assert len(list(my_emap.structure.get_residues()))==577 

def test_eta_moities():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1u3d.pdb"))  
    assert len(my_emap.eta_moieties) >= 3
    print(my_emap.eta_moieties)

def test_aromatic_amino_acid_options():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1u3d.pdb")) 
    pyemap.process(my_emap)
    ##default, tryptophan and tyrosine on
    resnames=[]
    for residue in my_emap.residues.values():
        resnames.append(residue.resname)
    assert "TYR" in resnames
    assert "TRP" in resnames
    assert "PHE" not in resnames
    assert "HIS" not in resnames
    #phenylalanine on
    pyemap.process(my_emap,phe=True)
    resnames=[]
    for residue in my_emap.residues.values():
        resnames.append(residue.resname)
    assert "TYR" in resnames
    assert "TRP" in resnames
    assert "PHE" in resnames
    assert "HIS" not in resnames
    #histidine on
    resnames=[]
    pyemap.process(my_emap,his=True)
    for residue in my_emap.residues.values():
        resnames.append(residue.resname)
    assert "TYR" in resnames
    assert "TRP" in resnames
    assert "PHE" not in resnames
    assert "HIS" in resnames
    #tyrosine off
    resnames=[]
    pyemap.process(my_emap,tyr=False)
    for residue in my_emap.residues.values():
        resnames.append(residue.resname)
    assert "TYR" not in resnames
    assert "TRP" in resnames
    assert "PHE" not in resnames
    assert "HIS" not in resnames
    #tryptophan off
    resnames=[]
    pyemap.process(my_emap,trp=False)
    for residue in my_emap.residues.values():
        resnames.append(residue.resname)
    assert "TYR" in resnames
    assert "TRP" not in resnames
    assert "PHE" not in resnames
    assert "HIS" not in resnames


def test_surface_options():
    rd_check = 'Y53(A)'
    asa_check = 'W61(A)' 
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1u3d.pdb")) 
    #residue depth
    pyemap.process(my_emap)
    G = my_emap.init_graph
    rd_ser = []
    for n, d in G.nodes(data=True):
        if d['shape'] == "box":
            rd_ser.append(n)
    #solvent accessibility
    pyemap.process(my_emap,sdef=1)
    G = my_emap.init_graph
    asa_ser = []
    for n, d in G.nodes(data=True):
        if d['shape'] == "box":
            asa_ser.append(n)
    assert rd_check in rd_ser
    assert asa_check not in rd_ser
    assert rd_check  not in asa_ser
    assert asa_check in asa_ser
    
def test_distance_options():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1u3d.pdb")) 
    pyemap.process(my_emap)
    com_weight = my_emap.init_graph["W385(A)"]["Y53(A)"]['weight']
    assert isclose(com_weight,6.18,abs_tol=1e-2)
    pyemap.process(my_emap,dist_def=1)
    closest_atom_weight = my_emap.init_graph["W385(A)"]["Y53(A)"]['weight']
    assert isclose(closest_atom_weight,3.97,abs_tol=1e-2)

def test_chains():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/2oal.pdb")) 
    #default all chains on
    pyemap.process(my_emap)
    labels=[]
    for residue in my_emap.residues.values():
        labels.append(residue.node_label)
    assert all(elem in labels for elem in ["W18(A)" and 'Y298(B)'])
    #B chain off
    pyemap.process(my_emap,chains=["A"])
    labels=[]
    for residue in my_emap.residues.values():
        labels.append(residue.node_label)
    assert "W18(A)" in labels
    assert 'Y298(B)' not in labels

def test_additional_residues():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1u3d.pdb")) 
    pyemap.process(my_emap)
    #default
    resnames=[]
    for residue in my_emap.residues.values():
        resnames.append(residue.resname)
    assert all(elem in resnames for elem in ["FAD510(A)-1","FAD510(A)-2","ANP511(A)"])
    #deselected eta moiety
    pyemap.process(my_emap,eta_moieties=["FAD510(A)-1","FAD510(A)-2"])
    resnames=[]
    for residue in my_emap.residues.values():
        resnames.append(residue.resname)
    assert all(elem in resnames for elem in ["FAD510(A)-1","FAD510(A)-2"])
    assert "ANP511(A)" not in resnames
    #custom residues
    pyemap.process(my_emap,eta_moieties=[], custom = "(3943-3952),(3953-3970)")
    resnames=[]
    for residue in my_emap.residues.values():
        resnames.append(residue.resname)
    assert all(elem in resnames for elem in ["CUST-1","CUST-2"])
    #example for bad custom residues, shoud raise exception
    try:
        pyemap.process(my_emap, custom = "(3943-3952),(3953-3970)")
        assert False
    except Exception as e:
        assert True

def test_graph_parameters():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1u3d.pdb")) 
    pyemap.process(my_emap,num_st_dev_edges=5,percent_edges=3)
    g1 = my_emap.init_graph
    # check st_dev_edges works
    pyemap.process(my_emap,num_st_dev_edges=1,percent_edges=3)
    g2 = my_emap.init_graph
    assert g1.edges() != g2.edges()
    #check distance cutoff works    
    pyemap.process(my_emap,num_st_dev_edges=5,percent_edges=3,distance_cutoff=30)
    g3 = my_emap.init_graph
    assert g1.edges() != g3.edges()
    #check percent edges works
    pyemap.process(my_emap,num_st_dev_edges=5,percent_edges=1,distance_cutoff=30)
    g4 = my_emap.init_graph
    assert g1.edges() != g4.edges()

def test_penalty_function_parameters():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1u3d.pdb")) 
    #alpha
    pyemap.process(my_emap,coef_alpha=0.5)
    alpha_weight = my_emap.init_graph["W385(A)"]["Y53(A)"]['weight']
    assert isclose(alpha_weight,6.48,abs_tol=1e-2)
    #Roffset
    pyemap.process(my_emap,r_offset=1)
    alpha_weight = my_emap.init_graph["W385(A)"]["Y53(A)"]['weight']
    assert isclose(alpha_weight,5.18,abs_tol=1e-2)
    #Roffset
    #pyemap.process(my_emap,exp_beta=1)
    #beta_weight = my_emap.init_graph["W385(A)"]["Y53(A)"]['weight']
    #assert isclose(beta_weight,13.44,abs_tol=1e-2)

def test_shortest_paths():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1u3d.pdb")) 
    pyemap.process(my_emap)
    pyemap.find_pathways(my_emap,"W385(A)")
    assert len(my_emap.paths)>0
    pyemap.find_pathways(my_emap,"W385(A)",target="Y53(A)")
    assert len(my_emap.paths)==10








