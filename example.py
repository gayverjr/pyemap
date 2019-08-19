import pyemap

my_emap = pyemap.fetch_and_parse("1u3d")
pyemap.process(my_emap)
pyemap.find_pathways(my_emap,"James")
#my_emap.save_residue("FAD510(A)-1")
#my_emap.show_init_graph()

'''
pyemap.find_pathways(my_emap,"Y330(A)",target="Y309(A)")
my_emap.save_paths_graph(dest="graph1.png")
pyemap.find_pathways(my_emap,"Y330(A)",target="W447(A)")
my_emap.save_paths_graph(dest="graph2.png")
'''
