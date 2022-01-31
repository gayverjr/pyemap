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

res_name_to_char = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    "ACE": ">",
    "NME": "<",
}

res_name_to_smiles = {
    'ALA': '[CH3]',
    'CYS': '*[CH2]',
    'ASP': '*C([CH2])=O',
    'GLU': '*C(=O)[CH2][CH2]',
    'PHE': '[CH2]c1[cH][cH][cH][cH][cH]1',
    'GLY': '*C(=O)[CH2]N',
    'HIS': '[CH2]C1:*:C:*:[CH]1',
    'ILE': '[CH3][CH][CH2][CH3]',
    'LYS': '*[CH2][CH2][CH2][CH2]',
    'LEU': '[CH2][CH]([CH3])[CH3]',
    'MET': '[CH2][CH2]S[CH3]',
    'ASN': '[CH2]C([NH2])=O',
    'PRO': '*C(=O)[CH]1*[CH2][CH2][CH2]1',
    'GLN': 'CCC(N)=O',
    'ARG': '*=C([NH2])[NH][CH2][CH2][CH2]',
    'SER': '[CH2][OH]',
    'THR': '[CH3][CH][OH]',
    'VAL': '[CH3][CH][CH3]',
    'TRP': '[CH2]c1[cH][nH]c2[cH][cH][cH][cH]c12',
    'TYR': '*c1[cH][cH]c([CH2])[cH][cH]1'
}

char_to_res_name = {v: k for k, v in res_name_to_char.items()}

aromatic_aa = ["HIS", "PHE", "TRP", "TYR"]
polar_aa = ["SER", "THR", "CYS", "PRO", "ASN", "GLN"]
pos_aa = ["LYS", "ARG"]
neg_aa = ["ASP", "GLU"]
nonpolar_aa = ["GLY", "ALA", "VAL", "LEU", "MET", "ILE"]

radii = {
    'H': 1.20,
    'N': 1.55,
    'NA': 2.27,
    'CU': 1.40,
    'CL': 1.75,
    'C': 1.70,
    'O': 1.52,
    'I': 1.98,
    'P': 1.80,
    'B': 1.85,
    'BR': 1.85,
    'S': 1.80,
    'SE': 1.90,
    'F': 1.47,
    'FE': 1.80,
    'K': 2.75,
    'MN': 1.73,
    'MG': 1.73,
    'ZN': 1.39,
    'HG': 1.80,
    'XE': 1.80,
    'AU': 1.80,
    'LI': 1.80,
    '.': 1.80
}

# pulled from: https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1992.pdf page 27
side_chain_atoms = {
    "ALA": ['CB'],
    "CYS": ['SG'],
    "ASP": ['CG', 'OD1', 'OD2'],
    "GLU": ['CG', 'CD', 'OE1', 'OE2'],
    "PHE": ['CG', 'CD1', 'CD2', 'CE1', 'CZ', 'CE2'],
    "GLY": ['CA'],
    "HIS": ['CG', 'ND1', 'CE1', 'NE2', 'CD2', 'AD1', 'AE1', 'AE2', 'AD2'],
    "ILE": ['CG1', 'CD1', 'CG2'],
    "LYS": ['CG', 'CD', 'CE', 'NZ'],
    "LEU": ['CG', 'CD1', 'CD2'],
    "MET": ['CG', 'SD', 'CE'],
    "ASN": ['CG', 'OD1', 'AD1', 'AD2', 'ND2'],
    "PRO": ['N', 'CA', 'CD', 'CG', 'CB'],
    "GLN": ['CG', 'CD', 'OE1', 'AE1', 'NE2', 'AE2'],
    "ARG": ['CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
    "SER": ['OG'],
    "THR": ['OG1', 'CG2'],
    "VAL": ['CG1', 'CG2'],
    "TRP": ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    "TYR": ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH']
}


SB_means = {
    'CC': 1.54,
    'CN': 1.49,
    'NC': 1.49,
    'CO': 1.43,
    'OC': 1.43,
    'CS': 1.81,
    'SC': 1.81,
    'CP': 1.87,
    'PC': 1.87,
    'NN': 1.43,
    'NO': 1.44,
    'ON': 1.44,
    'NP': 1.73,
    'PN': 1.73,
    'NS': 1.60,
    'SN': 1.60,
    'OS': 1.58,
    'SO': 1.58,
    'OO': 1.47,
    'OP': 1.62,
    'PO': 1.62,
    'SS': 2.05,
    'PP': 2.26
}

# http://hbcponline.com/faces/documents/09_01/09_01_0005.xhtml
SB_std_dev = {
    'CC': 0.015,
    'CN': 0.014,
    'NC': 0.014,
    'CO': 0.019,
    'OC': 0.019,
    'CS': 0.019,
    'SC': 0.019,
    'CP': 0.019,
    'PC': 0.019,
    'NN': 0.027,
    'NO': 0.009,
    'ON': 0.009,
    'NP': 0.017,
    'PN': 0.017,
    'SN': 0.012,
    'NS': 0.012,
    'OS': 0.015,
    'SO': 0.015,
    'OO': 0.012,
    'OP': 0.007,
    'PO': 0.007,
    'SS': 0.026,
    'PP': 0.025
}

import os
# clusters
clusters = os.listdir(os.path.abspath(os.path.dirname(__file__)) + '/data/clusters')
clusters = [cluster.replace('.svg', '') for cluster in clusters]

metal_ligands = {'CU1': 1, 'CU': 2, 'CU3': 3,  'FE': 3, 'FE2': 2, 'MN': 2, 'MN3': 3,   
'CO': 2, '3CO': 3, 'CR': 3,  'NI': 2, 'MO': 0, '4MO': 4, '6MO': 6}
