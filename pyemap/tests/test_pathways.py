## Copyright (c) 2017-2022, James Gayvert, Ruslan Tazhigulov, Ksenia Bravaya
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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











