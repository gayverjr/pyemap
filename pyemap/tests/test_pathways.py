import pyemap
import os
import unittest

class PathwaysTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.my_emap = pyemap.fetch_and_parse("2oal")
        pyemap.process(cls.my_emap,chains=["A","B"])
        cls.my_emap2 = pyemap.fetch_and_parse("1u3d")
        pyemap.process(cls.my_emap2)
        cls.my_emap3 = pyemap.fetch_and_parse("1u3d")
        pyemap.process(cls.my_emap3,sdef=None)
        os.remove("2oal.pdb") 

    def test_dijkstras_paths(self):
        pyemap.find_paths(self.my_emap,"W466(A)")
        assert len(self.my_emap.paths) == 5
        pyemap.find_paths(self.my_emap2,"W400(A)")
        assert len(self.my_emap2.paths) == 12
        try:
            pyemap.find_paths(self.my_emap3,"W400(A)")
            assert False
        except Exception:
            assert True

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











