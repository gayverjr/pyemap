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

"""Parser used to parse PDB and MMCIF files.

Constructs an emap object containing parsed Bio.PDB.Structure and a list of customized Bio.PDB.Residue objects
corresponding to automatically identified electron transfer moieties.
"""
from Bio.PDB import PDBIO, FastMMCIFParser, PDBParser
from .custom_residues import process_custom_residues
from .data import res_name_to_char
from .emap import emap
from .pyemap_exceptions import *
import os
from pathlib import Path
import requests


def _download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        R = requests.get(url, allow_redirects=True)
        if R.status_code != 200:
            raise ConnectionError('could not download {}\nerror code: {}'.format(url, R.status_code))
        Path(outfnm).write_bytes(R.content)
        return outfnm
    except Exception as e:
        print(e)
        raise PyeMapParseException("Could not fetch PDB {}.".format(pdbcode)) from e


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
        print("Fetching PDB " + pdb_id + " from RSCB Database...")
    outfnm = _download_pdb(pdb_id, dest)
    if not quiet:
        print("Success!")
    return parse(outfnm, quiet)


def parse(filename, quiet=True):
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
    if not quiet:
        print("Parsing file: " + str(filename))
    try:
        os.listdir()
        parser = PDBParser()
        structure = parser.get_structure("protein", filename)
    except Exception as e:
        try:
            parser = FastMMCIFParser()
            structure = parser.get_structure("protein", filename)
            io = PDBIO()
            fn = filename[:-4] + ".pdb"
            io.set_structure(structure)
            io.save(fn)
            parser = PDBParser()
            structure = parser.get_structure("protein", fn)
        except Exception:
            raise PyeMapParseException("Error: could not parse file {}.".format(filename)) from e
    chain_list = []
    sequences = {}
    non_standard_residue_list = []
    num_models = 0
    for model in structure.get_models():
        num_models += 1
    if num_models < 1:
        raise PyeMapParseException("Error: structure " + structure.header['idcode'] + " does not contain any models.")
    if structure.header['idcode'] == "":
        idcode = "CUST"
    else:
        idcode = structure.header['idcode']
    for chain in structure[0]:
        chain_list.append(chain.id)
        seq = []
        seq_idx = 1
        for residue in chain:
            if residue.resname in res_name_to_char:
                seq.append(res_name_to_char[residue.resname])
                residue.sequence_index = seq_idx
            else:
                residue.get_full_id()
                non_standard_residue_list.append(residue.copy())
                residue.sequence_index = 'X'
            seq_idx += 1
        seq_str = ">" + idcode + ":" + chain.id + "\n" + ''.join(seq)
        sequences[chain.id] = seq_str
    custom_residue_list = process_custom_residues(non_standard_residue_list)
    if not quiet:
        print("Identified " + str(len(custom_residue_list)) + " non-protein ET active moieties.")
    my_emap = emap(filename, idcode, custom_residue_list, chain_list, sequences)
    my_emap._structure = structure
    return my_emap
