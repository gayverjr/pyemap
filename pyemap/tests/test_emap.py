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