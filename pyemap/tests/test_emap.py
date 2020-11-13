import os
from os import path
import sys
import unittest
import warnings
import pyemap
import tempfile

def test_save_functions():
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/4DJA.pdb")) 
    #cluster
    fout = tempfile.NamedTemporaryFile(suffix=".png")
    my_emap.residue_to_Image("SF4603(A)")
    my_emap.residue_to_file("SF4603(A)",dest=fout.name)
    ''' commented out for now until we can find a solution with newer RDKit
    #aromatic eta moiety
    my_emap.residue_to_Image("FAD601(A)-1")
    my_emap.residue_to_file("FAD601(A)-1",dest=fout.name)
    '''
    #after file is processed
    pyemap.process(my_emap)
    #standard residue
    ''' commented out for now until we can find a solution with newer RDKit
    my_emap.residue_to_Image("Y443(A)")
    my_emap.residue_to_file("Y443(A)",dest=fout.name)
    '''
    #init graph
    my_emap.init_graph_to_Image()
    my_emap.init_graph_to_file(dest=fout.name)
    pyemap.find_paths(my_emap,"Y443(A)")
    #paths graph
    my_emap.paths_graph_to_Image()
    my_emap.paths_graph_to_file(dest=fout.name)
    #check that report does something
    assert my_emap.report() != None
    
