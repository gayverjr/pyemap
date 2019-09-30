import pyemap
my_emap = pyemap.fetch_and_parse("1u3d")
# view residue
my_emap.residue_to_Image("FAD510(A)-2").show()
# process file
pyemap.process(my_emap)
my_emap.init_graph_to_Image().show()
# find paths
pyemap.find_paths(my_emap,"FAD510(A)-2")
my_emap.paths_graph_to_Image().show()
# print report
print(my_emap.report())
