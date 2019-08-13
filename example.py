import pyemap

my_emap = pyemap.fetch_and_parse("1u3d")
#print(len(list(my_emap.structure.get_residues())))
#my_emap.save_residue("FAD510(A)-1")
#custom_atm_string="3943-3952"

#pyemap.process(my_emap)
print(my_emap.eta_moieties)

'''
pyemap.find_pathways(my_emap,"Y330(A)",target="Y309(A)")
my_emap.save_paths_graph(dest="graph1.png")
pyemap.find_pathways(my_emap,"Y330(A)",target="W447(A)")
my_emap.save_paths_graph(dest="graph2.png")
'''
