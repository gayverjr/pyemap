import pyemap
from pyemap.common_paths import PDBGroup
import time

'''
f = open("flavoproteins.txt")
lines = f.readlines()
f.close()
pdb_ids = []
for line in lines:
    l = line.strip()
    pdb_ids += l.split(',')
'''

pdb_list = "6v59,5nx0,1eb7,2vhd,1rz5,2c1u,2c1v,1iqc,1rz6,1nml,1zzh,1c1u,1c1v"
pdb_ids = pdb_list.split(',')

chains = {}

pg = PDBGroup("flavoproteins")
for i in range(0,len(pdb_ids)):
    emap_obj = pyemap.parse(pdb_ids[i]+".pdb")
    pg.add_emap(emap_obj)

#for pdb_id,emap in pg.emaps.items():
#    chains[pdb_id]=emap.chains

pg.process_emaps(dist_def=1,include_residues=["TRP","HIS","TYR"])
pg.generate_graph_database()
pg.find_subgraph('#W#')
sg = next(iter(pg.subgraph_patterns.items()))[1]
sg.find_protein_subgraphs()
print(sg.support)
from pandas import DataFrame
print(DataFrame(sg.L).to_string())
import numpy.linalg as LA
import numpy as np
L = sg.L
eigv, eigvc = LA.eig(L)
eigv = np.real(eigv)
eigvc = np.real(eigvc)
idx = eigv.argsort()
eigv = eigv[idx]
eigvc = eigvc[:, idx]
# second lowest eigenvector
eigvc2 = eigvc[:, 1]
print(eigvc2)
print(pg.general_report())
print(sg.full_report())
