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
pg.generate_graph_database(edge_thresholds=[8.5,12.0])
pg.run_gspan(int(0.5*len(pdb_ids)),lower_bound=4)

print(pg.frequent_subgraphs.keys())


my_sg = pg.frequent_subgraphs['7_WWWX_18']
id = my_sg.id
my_sg.clustering()


my_emap = pg.emaps["1U3D"]

for key in my_sg.eigenvector_sorted:
    print(str(key)+ " group:")
    graphs = my_sg.eigenvector_sorted[key]
    print(str(len(graphs))+ " members")
    for graph in graphs:
        print(my_sg._report_for_graph(graph))
        #img = my_emap._graph_to_Image(graph)
        #img.show()


