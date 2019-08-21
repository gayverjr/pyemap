import pyemap
import os
import sys
from math import isclose
import unittest

class PathwaysTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.my_emap = pyemap.parse(os.path.join(sys.path[0],"pyemap/tests/test_pdbs/2oal.pdb")) 
        pyemap.process(cls.my_emap)

    def test_dijkstras_paths(self):
        pyemap.find_paths(self.my_emap,"Y167(B)")
        assert len(self.my_emap.paths)>0

    def test_yens_paths(self):
        pyemap.find_paths(self.my_emap,"Y167(B)",target="Y369(B)")
        assert len(self.my_emap.paths)==10
        pyemap.find_paths(self.my_emap,"Y167(B)",target="Y369(B)",max_paths=15)
        assert len(self.my_emap.paths)==15
        try:
            pyemap.find_paths(self.my_emap,"Y167(B)",target="Y167(A)")
            assert False
        except Exception as e:
            assert True










