import os
import sys
import unittest
import warnings
import pyemap

def test_parse():
    #from PDB file
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1u3d.pdb"))  
    assert len(list(my_emap.structure.get_residues()))==577 
    #from CIF file
    my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/1u3d.cif"))  
    assert len(list(my_emap.structure.get_residues()))==577 


def test_fetch():
    #fetch from pdb database
    my_emap = pyemap.fetch_and_parse("1u3d")
    os.remove("1u3d.pdb")
    assert len(list(my_emap.structure.get_residues()))==577
