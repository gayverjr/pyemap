#!/usr/bin/env python3 
import pyemap 
from pyemap.graph_mining import PDBGroup 
# Step 1: parse PDBs 
pdb_ids = ['1X0P', '1DNP', '1EFP', '1G28', '1IQR', '1IQU', '1NP7', '1O96', '1O97', '1QNF', '1U3C', '1U3D', '2IYG', '2J4D', '2WB2', '2Z6C', '3FY4', '3ZXS', '4EER', '4GU5', '4I6G', '4U63', '6FN2', '6KII', '6LZ3', '6PU0', '6RKF'] 
pg = PDBGroup('My Group') 
for pdb in pdb_ids: 
    pg.add_emap(pyemap.fetch_and_parse(pdb)) 
# Step 2: Generate graphs 
process_kwargs = {'sdef': None, 'dist_def': 'COM', 'rsa_thresh': 0.05, 'rd_thresh': 3.03, 'distance_cutoff': 20.0, 'percent_edges': 1.0, 'max_degree': 4.0, 'edge_prune': 'DEGREE', 'num_st_dev_edges': 1.0, 'coef_alpha': 1.0, 'exp_beta': 2.3, 'r_offset': 0.0} 
chains = {'1X0P': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'], '1DNP': ['A', 'B'], '1EFP': ['A', 'B', 'C', 'D'], '1G28': ['A', 'B', 'C', 'D'], '1IQR': ['A'], '1IQU': ['A'], '1NP7': ['A', 'B'], '1O96': ['A', 'B', 'C', 'D', 'E', 'F', 'Q', 'Z'], '1O97': ['C', 'D'], '1QNF': ['A'], '1U3C': ['A'], '1U3D': ['A'], '2IYG': ['A', 'B'], '2J4D': ['A', 'B'], '2WB2': ['A', 'C', 'D'], '2Z6C': ['A', 'B'], '3FY4': ['A', 'B', 'C'], '3ZXS': ['A', 'B', 'C'], '4EER': ['A'], '4GU5': ['A', 'B'], '4I6G': ['A', 'B'], '4U63': ['A'], '6FN2': ['A'], '6KII': ['A'], '6LZ3': ['A', 'B', 'C', 'D'], '6PU0': ['A'], '6RKF': ['A', 'B', 'C', 'D', 'E', 'F']} 
include_residues = ['W', 'Y'] 
included_eta_moieties = {'1X0P': ['FAD9150(A)', 'FAD9151(B)', 'FAD9152(C)', 'FAD9153(D)', 'FAD9154(E)', 'FAD9155(F)', 'FAD9156(G)', 'FAD9157(H)', 'FAD9158(I)', 'FAD9159(J)'], '1DNP': ['FAD472(A)-1', 'FAD472(A)-2', 'MHF473(A)', 'FAD472(B)-1', 'FAD472(B)-2', 'MHF473(B)'], '1EFP': ['FAD399(A)-1', 'FAD399(A)-2', 'AMP400(A)', 'FAD399(C)-1', 'FAD399(C)-2', 'AMP400(C)'], '1G28': ['FMN1033(A)', 'FMN2033(B)', 'FMN3033(C)', 'FMN4033(D)'], '1IQR': ['FAD421(A)-1', 'FAD421(A)-2'], '1IQU': ['TDR500(A)', 'FAD421(A)-1', 'FAD421(A)-2'], '1NP7': ['FAD500(A)-1', 'FAD500(A)-2', 'FAD501(B)-1', 'FAD501(B)-2'], '1O96': ['AMP1263(A)', 'FAD1319(B)-1', 'FAD1319(B)-2', 'AMP1262(C)', 'FAD1319(D)-1', 'FAD1319(D)-2', 'AMP1262(E)', 'FAD1319(F)-1', 'FAD1319(F)-2', 'AMP1262(Q)', 'FAD1319(Z)-1', 'FAD1319(Z)-2'], '1O97': ['AMP1262(C)', 'FAD1319(D)-1', 'FAD1319(D)-2'], '1QNF': ['FAD485(A)-1', 'FAD485(A)-2', 'HDF486(A)'], '1U3C': ['FAD510(A)-1', 'FAD510(A)-2'], '1U3D': ['FAD510(A)-1', 'FAD510(A)-2', 'ANP511(A)'], '2IYG': ['FMN1122(A)', 'FMN1122(B)'], '2J4D': ['FAD1498(A)-1', 'FAD1498(A)-2', 'MHF1499(A)-1', 'MHF1499(A)-2', 'FAD1501(B)-1', 'FAD1501(B)-2', 'MHF1502(B)-1', 'MHF1502(B)-2'], '2WB2': ['FAD1510(A)-1', 'FAD1510(A)-2', 'DA1(C)', 'DC2(C)', 'DA3(C)', 'DG4(C)', 'DC5(C)', 'DG6(C)', 'DG7(C)', 'Z9(C)', 'DG10(C)', 'DC11(C)', 'DA12(C)', 'DG13(C)', 'DG14(C)', 'DT15(C)', 'DT1(D)', 'DA2(D)', 'DC3(D)', 'DC4(D)', 'DT5(D)', 'DG6(D)', 'DC7(D)', 'DG8(D)', 'DA9(D)', 'DC10(D)', 'DC11(D)', 'DG12(D)', 'DC13(D)', 'DT14(D)', 'DG15(D)'], '2Z6C': ['FMN500(A)', 'FMN500(B)'], '3FY4': ['IMD901(A)', 'IMD902(A)', 'IMD905(A)', 'FAD900(A)-1', 'FAD900(A)-2', 'IMD901(B)', 'IMD903(B)', 'IMD904(B)', 'IMD905(B)', 'IMD906(B)', 'FAD900(B)-1', 'FAD900(B)-2', 'IMD901(C)', 'IMD902(C)', 'FAD900(C)-1', 'FAD900(C)-2'], '3ZXS': ['FAD1509(A)-1', 'FAD1509(A)-2', 'DLZ1511(A)', 'FAD1509(B)-1', 'FAD1509(B)-2', 'DLZ1511(B)', 'FAD1509(C)-1', 'FAD1509(C)-2', 'DLZ1511(C)', 'SF41510(A)', 'SF41510(B)', 'SF41510(C)'], '4EER': ['FMN1001(A)'], '4GU5': ['FAD602(A)-1', 'FAD602(A)-2', 'FAD602(B)-1', 'FAD602(B)-2'], '4I6G': ['FAD900(A)-1', 'FAD900(A)-2', 'FAD900(B)-1', 'FAD900(B)-2'], '4U63': ['MHF1001(A)', 'FAD1002(A)-1', 'FAD1002(A)-2'], '6FN2': ['FAD601(A)-1', 'FAD601(A)-2', 'HDF602(A)'], '6KII': ['FAD501(A)-1', 'FAD501(A)-2', 'MHF502(A)-1', 'MHF502(A)-2'], '6LZ3': ['FAD701(A)-1', 'FAD701(A)-2', 'FAD701(B)-1', 'FAD701(B)-2', 'FAD701(C)-1', 'FAD701(C)-2', 'FAD701(D)-1', 'FAD701(D)-2'], '6PU0': ['FAD501(A)-1', 'FAD501(A)-2'], '6RKF': ['FAD401(A)-1', 'FAD401(A)-2', 'FAD401(B)-1', 'FAD401(B)-2', 'FAD401(C)-1', 'FAD401(C)-2', 'FAD401(D)-1', 'FAD401(D)-2', 'FAD401(E)-1', 'FAD401(E)-2', 'FAD401(F)-1', 'FAD401(F)-2']} 
pg.process_emaps(chains=chains, eta_moieties=included_eta_moieties, include_residues=include_residues, **process_kwargs) 
# Step 3: Generate graph database 
substitutions = [] 
edge_thresholds = [12] 
pg.generate_graph_database(sub=substitutions,edge_thresh=edge_thresholds) 
# Step 4: Mine for subgraphs 
pg.find_subgraph('WWW#') 
# Select subgraph pattern 
sg = pg.subgraph_patterns['1_WWW#_18'] 
# Final step: Find protein subgraphs 
sg.find_protein_subgraphs('structural')
# To switch clustering, call 'sg.set_clustering('sequence')
# Visualize pattern
sg.subgraph_to_Image().show()
# Visualize protein subgraph
sg.subgraph_to_Image(id='6KII_23').show()
# Print reports
print(pg.mining_report())
print(sg.full_report())
