# This Python script is a part of
# eMAP: online mapping of electron transfer channels in biomolecules
# Copyright(C) 2017-2018 Ruslan Tazhigulov, James Gayvert, Melissa Wei, Ksenia Bravaya (Boston University, USA)
"""Parser used to parse PDB and MMCIF files.

Main module used for step 1. Takes a filename input from the view, and returns a Bio.PDB structure object,
a list of chains, and a list of Bio.PDB residue objects.

Usage
-----
Called by the views module to parse the .pdb/.mmcif file and obtain the Structure object, chains, and
custom Residue objects.

See Also
________
custom_residues: used to identify custom residues
Bio.PDB: module used to handle all protein data
views: view functions which call parse
data: contains useful dictionaries

"""
from Bio.PDB import PDBIO, FastMMCIFParser, PDBParser
from .custom_residues import *
from .data import *
import os


def upload(filename):
    """Parses file and returns structure, chain list, and custom residue list.

    Uses the Bio.PDB module to parse the file and obtain a Structure object.
    The chains present in the structure are identified, and then the custom_residues
    module is used to identify any custom residues present in the structure. If
    a .mmcif file is processed, a .pdb copy is also written to file and stored
    in the appropriate tempFiles directory on the server.

    Parameters
    ----------
    filename: string
        file specified by user

    Returns
    -------
    (structure,chain_list,custom_residue_list): tuple
        A Bio.PDB structure object, a list of string chain IDs,
        a list of Bio.PDB Residue objects identified as custom residues.

    Raises
    ______
    Exception
        Bad PDB/mmcif file.

    Notes
    ------
    For .pdb/.mmcif files with multiple models, only the first model is processed by eMap.

    Bio.PDB.FastMMCIFParser does not store atom serial numbers in the Residue objects, so
    a .pdb copy must be written to file and parsed in order to use atom serial numbers later.

    """
    try:
        parser = PDBParser()
        pdb_file = filename + "/file" + filename + ".pdb"
        structure = parser.get_structure("protein", pdb_file)
    except Exception as e:
        parser = FastMMCIFParser()
        cif_file = filename + "/file" + filename + ".cif"
        structure = parser.get_structure("protein", cif_file)
        io = PDBIO()
        io.set_structure(structure)
        io.save(filename + "/file" + filename + ".pdb")
        parser = PDBParser()
        pdb_file = filename + "/file" + filename + ".pdb"
        structure = parser.get_structure("protein", pdb_file)
    chain_list = []
    num_models = 0
    for model in structure.get_models():
        num_models += 1
    if num_models < 1:
        raise Exception(
            "Error Code 1: Unable to parse PDB. Please try another file.")
    for chain in structure[0].get_chains():
        chain_list.append(chain.id)
    non_standard_residue_list = []
    for res in structure[0].get_residues():
        if res.resname not in data.res_name_to_char:
            res.get_full_id()
            arom_res = res.copy()
            non_standard_residue_list.append(arom_res)
    custom_residue_list = process_custom_residues(
        non_standard_residue_list, filename, len(chain_list),"")
    return structure, chain_list, custom_residue_list

def fetch(filename):
    print("Fetching file " + filename + " from rcsb database...")
    cmd = (
           'wget --no-check-certificate --read-timeout=1 -t 1 -nc -P {0} https://files.rcsb.org/download/'
           + filename)
    os.system(cmd)
    print("Success!")

