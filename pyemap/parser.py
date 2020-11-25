# PyeMap: A python package for automatic identification of electron and hole transfer pathways in proteins.
# Copyright(C) 2017-2020 Ruslan Tazhigulov, James Gayvert, Ksenia Bravaya (Boston University, USA)
"""Parser used to parse PDB and MMCIF files. 

Constructs an emap object containing parsed Bio.PDB.Structure and a list of customized Bio.PDB.Residue objects
corresponding to automatically identified electron transfer moieties.
"""
from Bio.PDB import PDBIO, FastMMCIFParser, PDBParser
from .custom_residues import process_custom_residues
from .data import res_name_to_char
from .emap import emap
import os

def fetch_and_parse(pdb_id, dest="", quiet=False):
    '''Fetches pdb from database and parses the file.

    Parameters
    ----------
    pdb_id: str
        RCSB PDB ID
    dest: str, optional
        Full path to where file should be saved
    quiet: bool, optional
        Supresses output when set to true

    Returns
    -------
    emap: :class:`~pyemap.emap`
        emap object ready for processing.
    '''
    if not dest:
        dest = os.getcwd()
    if not quiet:
        print("Fetching file " + pdb_id + " from RSCB Database...")
    cmd = ('wget -nc --no-check-certificate --quiet --read-timeout=1 -t 1 -P {0} https://files.rcsb.org/download/' +
    pdb_id + ".pdb").format(dest)
    import subprocess
    subprocess.check_output(cmd, shell=True)
    if not quiet:
        print("Success!")
    return parse(dest + "/" + pdb_id + ".pdb", quiet)

def parse(filename, quiet=False):
    '''Parses pdb file and returns emap object.
    
    Parameters
    ----------
    filename: str
        Full path to file which needs to be parsed
    quiet: bool, optional
        Supresses output when set to true
    
    Returns
    -------
    my_emap: :class:`~pyemap.emap`
        emap object reading for parsing
    '''
    try:
        parser = PDBParser()
        structure = parser.get_structure("protein", filename)
    except Exception as e:
        parser = FastMMCIFParser()
        structure = parser.get_structure("protein", filename)
        io = PDBIO()
        fn = filename[:-4] + ".pdb"
        io.set_structure(structure)
        io.save(fn)
        parser = PDBParser()
        structure = parser.get_structure("protein", fn)
    chain_list = []
    num_models = 0
    for model in structure.get_models():
        num_models += 1
    if num_models < 1:
        raise RuntimeError("Unable to parse file.")
    for chain in structure[0].get_chains():
        chain_list.append(chain.id)
    non_standard_residue_list = []
    for res in structure[0].get_residues():
        if res.resname not in res_name_to_char:
            res.get_full_id()
            arom_res = res.copy()
            non_standard_residue_list.append(arom_res)
    custom_residue_list = process_custom_residues(non_standard_residue_list)
    if not quiet:
        print("Identified " + str(len(custom_residue_list)) + " non-protein ET active moieties.")
    my_emap = emap(filename, structure, custom_residue_list, chain_list)
    return my_emap
