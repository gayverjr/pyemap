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
pg.find_subgraph('WWWX')


print(pg.subgraph_report('0_WWWX_18'))

fs = pg.subgraph_patterns['0_WWWX_18']
for key,val in fs.protein_subgraphs.items():
    for key2,val2 in fs.protein_subgraphs.items():
        if not key==key2:
            print(str(key) + " group:" + str(val.graph['group_val']))
            print(str(key2) + " group:" + str(val2.graph['group_val']))
            print(pg.subgraph_rmsd(fs.id,key,key2))
            print()



