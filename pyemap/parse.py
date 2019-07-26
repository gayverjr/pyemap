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
from .emap import *

def fetch_and_parse(filename,dest=os.getcwd(),quiet=False,pdb=False):
    if not quiet:
        print("Fetching file " + filename + " from RSCB Database...")
    if not pdb:
        cmd = (
            'wget --no-check-certificate --quiet --read-timeout=1 -t 1 -nc -P {0} https://files.rcsb.org/download/'
            + filename+".cif").format(dest)
        os.system(cmd)
    if os.path.exists(dest+"/"+filename+".cif"):
        if not quiet:
            print("Success!")
        return parse(dest+"/"+filename+".cif",quiet,pdb)
    else:
        if not quiet and not pdb:
            print("Couldn't find .cif, trying .pdb ...")
        cmd = (
        'wget --no-check-certificate --quiet --read-timeout=1 -t 1 -nc -P {0} https://files.rcsb.org/download/'
        + filename+".pdb").format(dest)
        if os.path.exists(dest+"/"+filename+".pdb"):
            if not quiet:
                print("Success!")
            return parse(dest+"/"+filename+".pdb",quiet,pdb)
        else:
            raise Exception("Couldn't find entry matching PDB ID:" + str(filename))

def parse(filename,quiet=False,pdb=False):
    if pdb:
        parser = PDBParser()
        structure = parser.get_structure("protein", filename)
    else:
        parser = FastMMCIFParser()
        structure = parser.get_structure("protein", filename)
        io = PDBIO()
        fn=filename[:-4]+".pdb"
        io.set_structure(structure)
        io.save(fn)
        parser = PDBParser()
        structure = parser.get_structure("protein", fn)
    chain_list = []
    num_models = 0
    for model in structure.get_models():
        num_models += 1
    if num_models < 1:
        raise Exception(
            "Unable to parse structure. Please try another file.")
    for chain in structure[0].get_chains():
        chain_list.append(chain.id)
    non_standard_residue_list = []
    for res in structure[0].get_residues():
        if res.resname not in res_name_to_char:
            res.get_full_id()
            arom_res = res.copy()
            non_standard_residue_list.append(arom_res)
    custom_residue_list,smiles_str_list = process_custom_residues(
        non_standard_residue_list,len(chain_list))
    if not quiet:
        print("Identified " + str(len(custom_residue_list)) + " nonprotein ET active moieties." )
    my_emap = emap(filename[:-4],structure,chain_list,custom_residue_list,smiles_str_list)
    return my_emap

