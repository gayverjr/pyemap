from Bio.PDB import PDBIO, FastMMCIFParser, PDBParser
import pyemap
import os
import sys

def test_fetch():
    my_emap = pyemap.fetch_and_parse("1u3d")
    os.remove("1u3d.pdb")
    assert len(list(my_emap.structure.get_residues()))==577 
    assert len(my_emap.eta_moieties) == 3

def test_parse():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/Test1/1u3d.pdb"))  
    assert len(list(my_emap.structure.get_residues()))==577 
    assert len(my_emap.eta_moieties) == 3

def test_aromatic_amino_acid_options():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/Test1/1u3d.pdb")) 
    pyemap.process(my_emap)
    default = len(my_emap.residues)
    assert default == 39  
    pyemap.process(my_emap,phe=True)
    phe_on = len(my_emap.residues)
    assert phe_on == 63 
    pyemap.process(my_emap,his=True)
    his_on = len(my_emap.residues)
    assert his_on == 53 
    pyemap.process(my_emap,tyr=False)
    trp_off = len(my_emap.residues)
    assert trp_off == 24 
    pyemap.process(my_emap,trp=False)
    tyr_off = len(my_emap.residues)
    assert tyr_off == 18


def test_surface_options():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/Test1/1u3d.pdb")) 
    correct_rd_ser = ['W45(A)', 'Y53(A)', 'Y111(A)', 'Y142(A)', 'W145(A)', 'W162(A)', 'Y170(A)', 'W217(A)', 'W275(A)', 'Y299(A)', 'Y302(A)', 
    'Y309(A)', 'Y330(A)', 'Y341(A)', 'W352(A)', 'W356(A)', 'W377(A)', 'W379(A)', 'Y383(A)', 'Y402(A)', 'Y424(A)', 'W436(A)', 'W447(A)', 'W452(A)', 'W492(A)', 'ANP511(A)']
    pyemap.process(my_emap)
    G = my_emap.init_graph
    rd_ser = []
    for n, d in G.nodes(data=True):
        if d['shape'] == "box":
            rd_ser.append(n)
    assert correct_rd_ser == rd_ser

    correct_asa_ser = ['W45(A)', 'W61(A)', 'W145(A)', 'W162(A)', 'Y170(A)', 'W213(A)', 'W217(A)', 'Y235(A)', 'W275(A)', 'Y309(A)', 'W324(A)', 
    'Y330(A)', 'W352(A)', 'W356(A)', 'W377(A)', 'W379(A)', 'Y383(A)', 'W385(A)', 'Y402(A)', 'Y424(A)', 'W436(A)', 'W447(A)', 'W452(A)', 'W492(A)']
    pyemap.process(my_emap,sdef=1)
    G = my_emap.init_graph
    asa_ser = []
    for n, d in G.nodes(data=True):
        if d['shape'] == "box":
            asa_ser.append(n)
    assert asa_ser == correct_asa_ser
    
def test_distance_options():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/Test1/1u3d.pdb")) 
    pyemap.process(my_emap)
    com_graph = my_emap.init_graph
    assert len(com_graph.edges()) == 52
    pyemap.process(my_emap,dist_def=1)
    closest_graph = my_emap.init_graph
    assert len(closest_graph.edges())== 52
    assert com_graph.edges() != closest_graph.edges()

def test_additional_residues():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/Test1/1u3d.pdb")) 
    pyemap.process(my_emap,eta_moieties=["FAD510(A)-1","FAD510(A)-2"])
    assert len(my_emap.residues) == 38
    pyemap.process(my_emap,eta_moieties=[], custom = "(3943-3952),(3953-3970)")
    assert len(my_emap.residues) == 38









