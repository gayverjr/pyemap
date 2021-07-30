import os
from os import path
import sys
import unittest
import warnings
import pyemap
import tempfile
from sys import platform

def test_save_functions():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/4DJA.pdb")) 
    #cluster
    fout = tempfile.NamedTemporaryFile(suffix=".png")
    fout_svg = tempfile.NamedTemporaryFile(suffix=".svg")
    #aromatic eta moiety
    my_emap.residue_to_Image("SF4603(A)")
    my_emap.residue_to_file("SF4603(A)",dest=fout.name)
    my_emap.residue_to_Image("FAD601(A)-1")
    my_emap.residue_to_file("FAD601(A)-1",dest=fout.name)
    my_emap.residue_to_file("FAD601(A)-1",dest=fout_svg.name)
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
    
