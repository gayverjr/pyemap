import pyemap
from pyemap.common_paths import protein_group
from pyemap.common_paths.protein_group import Subgraph




protein_ids = ["1u3d","1u3c","6PU0","4I6G","2J4D"]

pg = protein_group("cryptochromes")#temp_dir="/Users/JG/Documents/Software/pyemap/test_dir/tmp_dir")

for i in range(0,len(protein_ids)):
    emap_obj = pyemap.fetch_and_parse(protein_ids[i])
    pg.add_emap(emap_obj)

for emap_id in pg.emaps.keys():
    cur_emap = pg.emaps[emap_id]
    pyemap.process(cur_emap,dist_def=1)

pg.generate_graph_db(edge_thresholds=[8.0,12.0])
pg.run_gspan(support=5,lower_bound=6)

#pg.subgraphs is a dictionary
# grab subgraph with ID 0_WWWWWY_5

report= open("full_report.txt", "w")
for graph_id in pg.subgraphs:
    report.write("\n")
    report.write("subgraph " + pg.subgraphs[graph_id].id)
    report.write("\n")
    my_sg = pg.subgraphs[graph_id]



 #   print("generic example")
  #  print(my_sg.generate_generic_report())
   # print("specific example")
    report.write(my_sg.specific_subgraph_example())

