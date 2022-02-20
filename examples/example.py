# PyeMap: A python package for automatic identification of electron and hole transfer pathways in proteins.
# Copyright(C) 2017-2020 Ruslan Tazhigulov, James Gayvert, Ksenia Bravaya (Boston University, USA)
import pyemap
my_emap = pyemap.fetch_and_parse("2ij2")
# view residue
#my_emap.residue_to_Image("FAD510(A)-2").show()
# process file

pyemap.process(my_emap,eta_moieties=[],chains=['A'])
my_emap.init_graph_to_Image().show()
pyemap.find_paths(my_emap,'W96(A)','Y334(A)')
my_emap.paths_graph_to_Image().show()
print(my_emap.report())