import pyemap

my_emap = pyemap.fetch_and_parse("1u3d")
my_emap.save_residue("FAD510(A)-1")
pyemap.process(my_emap,my_emap.chains,my_emap.eta_moieties)
pyemap.find_pathways(my_emap,"Y330(A)",target="Y309(A)")
my_emap.save_paths_graph(dest="graph1.png")
pyemap.find_pathways(my_emap,"Y330(A)",target="W447(A)")
my_emap.save_paths_graph(dest="graph2.png")
#for pt in my_emap.shortest_paths:
#    print(pt)
