import pyemap
from pyemap.common_paths import PDBGroup
import time


f = open("flavoproteins.txt")
lines = f.readlines()
f.close()
pdb_ids = []
for line in lines:
    l = line.strip()
    pdb_ids += l.split(',')

pg = PDBGroup("flavoproteins")
for i in range(0,len(pdb_ids)):
    emap_obj = pyemap.fetch_and_parse(pdb_ids[i])
    pg.add_emap(emap_obj)

pg.process_emaps(sdef=None,dist_def=0)
#pg.generate_graph_database()
pg.generate_graph_database()
pg.find_subgraph('WW*#')
print(pg.subgraph_patterns)

# get graph w/ largest support
sg = next(iter(pg.subgraph_patterns.items()))[1]
sg.find_protein_subgraphs()
print("RMSD between: 1U3D(1)-1 and 1U3C(1)-1")
print(sg.subgraph_rmsd("1U3D(1)-1","1U3C(1)-1"))
print("Report")
print(sg.full_report())

