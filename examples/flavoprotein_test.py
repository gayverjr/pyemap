import pyemap
from pyemap.common_paths import PDBGroup
import time


def compare_graph_strings(str1, str2):
    if not len(str1) == len(str2):
        return False
    for i in range(0, len(str1)):
        if not str1[i] == str2[i] and not str1[i] == "*" and not str2[i] == "*":
            return False
    return True

def apply_filter(subgraph_patterns,filter_str):
    matched_subgraphs = {}
    for key in subgraph_patterns.keys():
        idx1 = key.index("_")
        idx2 = key.index("_", idx1 + 1)
        graph_str = key[idx1 + 1:idx2]
        if compare_graph_strings(graph_str, filter_str) or compare_graph_strings(graph_str[::-1], filter_str):
            matched_subgraphs[key] = subgraph_patterns[key]
    return matched_subgraphs

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

pg.process_emaps(dist_def=1,sdef=None,percent_edges=20,include_residues=["TYR","TRP"])
pg.generate_graph_database()
pg.run_gspan(int(0.5*len(pdb_ids)),lower_bound=4)
print(len(pg.subgraph_patterns))
filtered_graphs = apply_filter(pg.subgraph_patterns,"XWWW")

sg = filtered_graphs['111_WWWX_17']
sg.find_protein_subgraphs()
print(sg.full_report())

