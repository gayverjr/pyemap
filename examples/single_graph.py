import pyemap
from pyemap.common_paths import PDBGroup
from pyemap.common_paths.protein_group import FrequentSubgraph
import time


f = open("cryptochromes.txt")
lines = f.readlines()
f.close()
pdb_ids = []
for line in lines:
    l = line.strip()
    pdb_ids += l.split(',')

pg = PDBGroup("cryptochromes",temp_dir="/Users/JG/Documents/Software/pyemap/test_dir/tmp_dir")
for i in range(0,len(pdb_ids)):
    emap_obj = pyemap.fetch_and_parse(pdb_ids[i])
    pg.add_emap(emap_obj)

pg.process_emaps(dist_def=1,sdef=None)
pg.generate_graph_database()
pg.find_subgraph('WWWX')

print(pg.subgraph_report('1_WWWX_18'))

