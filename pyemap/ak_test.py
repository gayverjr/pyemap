import propka.atom
import propka.protonate
import logging
from pathlib import Path
from io import StringIO
import pytest
import propka.run
from propka.parameters import Parameters
from propka.molecular_container import MolecularContainer
from propka.input import read_parameter_file, read_molecule_file,read_pdb
from propka.output import open_file_for_writing, write_pdb_for_conformation,write_pdb_for_protein 
from propka.lib import loadOptions
import propka.protonate
import propka.run

         #   propka.run.single("protein.pdb",
          
          #      optargs=["--mutation=N25R/N181D", "-v", "--pH=7.2"])
def protonate(pdbid, chain):
	mymolecule=propka.run.single("{}.pdb".format(pdbid),  optargs=["-c {}".format(chain), "--protonate-all"],stream=None, write_pka=False)
	write_pdb_for_conformation(mymolecule.conformations['1A'], '{}_protonated.txt'.format(pdbid))
	write_pdb_for_conformation(mymolecule.conformations['1B'], '{}_protonated.txt'.format(pdbid))

#print(mymolecule.conformation_names)
#print(mymolecule.conformations)
#propka.protonate.Protonate(mymolecule)
#print(mymolecule)
#write_pdb_for_conformation(conformation, filename)

#write_pdb_for_conformation(mymolecule, chains= mymolecule.options.chains, 'output.txt')


#propka.input.read_molecule_file("1u3d.pdb", mymolecule, stream=None)
#mymolecule= propka.molecular_container.MolecularContainer("1u3d",options=None)
#newmolecule= propka.hydrogens.protonate_30_style(mymolecule)


#molecule = MolecularContainer()

#toprotonate= read_pdb('1u3d.pdb', parameters, molecule)





#propka.run.single(filename: str, optargs: tuple = (), stream=None, write_pka: bool = True)[source]

#H.protonate(molecule)

#H.molecule

#propka.output.write_pdb_for_conformation(H,'output.pdb')

#print(H)









