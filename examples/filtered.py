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

pg.process_emaps(sdef=None,dist_def=1)
pg.generate_graph_database()
pg.run_gspan(18)
print("Full set of patterns:")
print(list(pg.subgraph_patterns.keys()))
matched_patterns = pg.apply_filter("W**X")
print("Matched patterns:")
print(list(matched_patterns.keys()))
