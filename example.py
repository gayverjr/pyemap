import pyemap

my_emap = pyemap.fetch_and_parse("1u3d")
pyemap.process(my_emap)
#my_emap.save_init_graph(dest="1u3d.png")
#pyemap.process(my_emap,phe=True,distance_cutoff=15)
#my_emap.save_init_graph(dest="phe.png")
pyemap.find_paths(my_emap,"FAD510(A)-2")
#my_emap.save_paths_graph(dest="source_only.png")
pyemap.find_paths(my_emap,"FAD510(A)-2",target="W324(A)")
#my_emap.save_paths_graph(dest="target.png")
#my_emap.show_paths_graph()
#print(my_emap.paths["1a"])
#print(my_emap.paths["1a"].selection_strs)
weight = my_emap.init_graph["FAD510(A)-2"]["W400(A)"]['weight']
print(weight)
#my_emap.save_residue("FAD510(A)-1")
#my_emap.show_init_graph()

'''
pyemap.find_pathways(my_emap,"Y330(A)",target="Y309(A)")
my_emap.save_paths_graph(dest="graph1.png")
pyemap.find_pathways(my_emap,"Y330(A)",target="W447(A)")
my_emap.save_paths_graph(dest="graph2.png")
'''
