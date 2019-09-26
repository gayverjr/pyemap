import os
import sys
import unittest
import warnings
import pyemap
import tempfile

def test_save_functions():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1dnp.pdb")) 
    #cluster
    fout = tempfile.NamedTemporaryFile(suffix=".png")
    my_emap.residue_to_Image("MHF473(B)")
    my_emap.residue_to_file("MHF473(B)",dest=fout.name)
    #aromatic eta moiety
    my_emap.residue_to_Image("FAD472(B)-1")
    my_emap.residue_to_file("FAD472(B)-1",dest=fout.name)
    #after file is processed
    pyemap.process(my_emap)
    #standard residue
    my_emap.residue_to_Image("W306(A)")
    my_emap.residue_to_file("FAD472(B)-1",dest=fout.name)
    #init graph
    my_emap.init_graph_to_Image()
    my_emap.init_graph_to_file(dest=fout.name)
    pyemap.find_paths(my_emap,"FAD472(A)-2")
    #paths graph
    my_emap.paths_graph_to_Image()
    my_emap.paths_graph_to_file(dest=fout.name)
    assert True
    
