import os
import pyemap
import tempfile
from sys import platform

def test_save_functions():
    my_emap = pyemap.fetch_and_parse("4dja")
    #cluster
    fout = tempfile.NamedTemporaryFile(suffix=".png")
    #aromatic eta moiety
    my_emap.residue_to_Image("SF4603(A)")
    my_emap.residue_to_file("SF4603(A)",dest=fout.name)
    my_emap.residue_to_Image("FAD601(A)-1")
    my_emap.residue_to_file("FAD601(A)-1",dest=fout.name)
    if platform == "linux":
        pyemap.process(my_emap,sdef=1)
    elif platform == "darwin":
        pyemap.process(my_emap)
    #standard residue
    my_emap.residue_to_Image("Y443(A)")
    my_emap.residue_to_file("Y443(A)",dest=fout.name)
    #init graph
    my_emap.init_graph_to_Image()
    my_emap.init_graph_to_file(dest=fout.name)
    pyemap.find_paths(my_emap,"Y443(A)",target="Y437(A)")
    #paths graph
    my_emap.paths_graph_to_Image()
    my_emap.paths_graph_to_file(dest=fout.name)
    #check that report does something
    assert my_emap.report() != None
    os.remove("4dja.pdb") 

def test_ligands():
    my_emap = pyemap.fetch_and_parse("1A4A")
    #cluster
    fout = tempfile.NamedTemporaryFile(suffix=".png")
    #aromatic eta moiety
    my_emap.residue_to_Image("CU130(A)")
    my_emap.residue_to_file("CU130(A)",dest=fout.name)
    my_emap.residue_to_Image("CU130(A)")
    my_emap.residue_to_file("CU130(A)",dest=fout.name)
    if platform == "linux":
        pyemap.process(my_emap,sdef=1)
    elif platform == "darwin":
        pyemap.process(my_emap)
    pyemap.find_paths(my_emap,"CU130(A)",target="W48(A)")
    #check that report does something
    assert my_emap.report() != None
    os.remove("1A4A.pdb") 