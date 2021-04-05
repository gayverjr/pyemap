import pyemap
from pyemap.common_paths import PDBGroup
from pyemap.common_paths.protein_group import FrequentSubgraph
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

pg.process_emaps(dist_def=1,sdef=None)
pg.generate_graph_database()
pg.find_subgraph('WWWX')

print(pg.subgraph_report('1_WWWX_18'))

fs = pg.frequent_subgraphs['1_WWWX_18']
for key,val in fs.specific_subgraphs.items():
    for key2,val2 in fs.specific_subgraphs.items():
        if not key==key2:
            print(str(key) + " group:" + str(val.graph['group_val']))
            print(str(key2) + " group:" + str(val2.graph['group_val']))
            print(pg.subgraphs_rmsd(fs.id,key,key2))
            print()



