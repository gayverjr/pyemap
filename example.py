import pyemap

my_emap = pyemap.fetch_and_parse("1u3d")
pyemap.process(my_emap)
pyemap.find_paths(my_emap,"FAD510(A)-2",target = "W324(A)")
my_emap.report()
#my_emap.save_init_graph("custom.png")



'''
pyemap.find_pathways(my_emap,"Y330(A)",target="Y309(A)")
my_emap.save_paths_graph(dest="graph1.png")
pyemap.find_pathways(my_emap,"Y330(A)",target="W447(A)")
my_emap.save_paths_graph(dest="graph2.png")
'''
