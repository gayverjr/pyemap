# PyeMap: A python package for automatic identification of electron and hole transfer pathways in proteins.
# Copyright(C) 2017-2020 Ruslan Tazhigulov, James Gayvert, Ksenia Bravaya (Boston University, USA)
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

TRP_sc = ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
TYR_sc = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH']
PHE_sc = ['CG', 'CD1', 'CD2', 'CE1', 'CZ', 'CE2']
HIS_sc = ['CG', 'ND1', 'CE1', 'NE2', 'CD2', 'AD1', 'AE1', 'AE2', 'AD2']

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
clusters = os.listdir(os.path.abspath(
    os.path.dirname(__file__)) + '/data/clusters')
clusters = [cluster.replace('.svg', '') for cluster in clusters]
