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

pg.process_emaps(dist_def=0,include_residues=["TYR","TRP"])
pg.generate_graph_database()
pg.run_gspan(18,4)

print(pg.general_report())
sg = next(iter(pg.frequent_subgraphs.items()))[1]
sg.find_protein_subgraphs()
print(sg.full_report())

