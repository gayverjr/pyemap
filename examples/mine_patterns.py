#!/usr/bin/env python3 
import pyemap 
from pyemap.graph_mining import PDBGroup 
# Step 1: parse PDBs 
pdb_ids = ['1C52', '1CC5', '1COT', '1CTJ', '1E29', '1FCD', '1GKS', '1GU2', '1H1O', '1H32', '1I8O', '1IQC', '1JMX', '1KB0', '1KV9', '1M70', '1M7O', '1PBY', '1PPJ', '1QKS', '1QL3', '1RZ5', '1WVE', '1YCC', '2BH4', '2C1V', '2C85', '2FW5', '2FYN', '2VHD', '2XIL', '2XTS', '351C', '3A9F', '3AYF', '3HQ6', '3MK7', '3O5C', '3VRD', '3WFB', '4B2N', '4D6T', '4FA1', '5Z25', '6QKN', '6V59'] 
pg = PDBGroup('My Group') 
for pdb in pdb_ids: 
    pg.add_emap(pyemap.fetch_and_parse(pdb)) 
# Step 2: Generate graphs 
process_kwargs = {'sdef': None, 'dist_def': 'CATM', 'rsa_thresh': 0.05, 'rd_thresh': 3.03, 'distance_cutoff': 20.0, 'percent_edges': 1.0, 'max_degree': 4.0, 'edge_prune': 'DEGREE', 'num_st_dev_edges': 1.0, 'coef_alpha': 1.0, 'exp_beta': 2.3, 'r_offset': 0.0} 
chains = {'1C52': ['A'], '1CC5': ['A'], '1COT': ['A'], '1CTJ': ['A'], '1E29': ['A'], '1FCD': ['A', 'C', 'B', 'D'], '1GKS': ['A'], '1GU2': ['A', 'B'], '1H1O': ['A', 'B'], '1H32': ['A', 'B'], '1I8O': ['A'], '1IQC': ['A', 'B', 'C', 'D'], '1JMX': ['A', 'B', 'G'], '1KB0': ['A'], '1KV9': ['A'], '1M70': ['A', 'B', 'C', 'D'], '1M7O': ['A', 'B'], '1PBY': ['A', 'B', 'C'], '1PPJ': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W'], '1QKS': ['A', 'B'], '1QL3': ['A', 'B', 'C', 'D'], '1RZ5': ['A'], '1WVE': ['A', 'C', 'B', 'D'], '1YCC': ['A'], '2BH4': ['X'], '2C1V': ['A', 'B'], '2C85': ['A'], '2FW5': ['A'], '2FYN': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R'], '2VHD': ['A', 'B'], '2XIL': ['A'], '2XTS': ['A', 'B', 'C', 'D'], '351C': ['A'], '3A9F': ['A', 'B'], '3AYF': ['A'], '3HQ6': ['A', 'B'], '3MK7': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'U', 'X', 'Y', 'Z'], '3O5C': ['A', 'B', 'C', 'D'], '3VRD': ['A', 'B'], '3WFB': ['L', 'H', 'B', 'C'], '4B2N': ['A', 'B'], '4D6T': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W'], '4FA1': ['A', 'B', 'C', 'D', 'E', 'F'], '5Z25': ['A'], '6QKN': ['A', 'B'], '6V59': ['A', 'B']} 
include_residues = ['H', 'C'] 
included_eta_moieties = {'1C52': ['HEM200(A)'], '1CC5': ['HEM1(A)'], '1COT': ['HEC130(A)'], '1CTJ': ['HEM91(A)'], '1E29': ['HEC201(A)'], '1FCD': ['FAD699(A)-1', 'FAD699(A)-2', 'HEC901(C)', 'HEC902(C)', 'FAD699(B)-1', 'FAD699(B)-2', 'HEC901(D)', 'HEC902(D)'], '1GKS': ['HEM0(A)'], '1GU2': ['HEC125(A)', 'HEC125(B)'], '1H1O': ['HEM1184(A)', 'HEM1185(A)', 'HEM1385(B)', 'HEM1386(B)'], '1H32': ['HEC1263(A)', 'HEC1264(A)', 'HEC1139(B)'], '1I8O': ['HEC115(A)'], '1IQC': ['HEC401(A)', 'HEC402(A)', 'HEC401(B)', 'HEC402(B)', 'HEC401(C)', 'HEC402(C)', 'HEC401(D)', 'HEC402(D)'], '1JMX': ['HEC1001(A)', 'HEC1002(A)', 'TRQ43(G)', 'NI2001(A)'], '1KB0': ['TRO512(A)', 'HEC802(A)', 'PQQ1800(A)'], '1KV9': ['PQQ801(A)-1', 'PQQ801(A)-2', 'HEC901(A)'], '1M70': ['HEC199(A)', 'HEC200(A)', 'HEC199(B)', 'HEC200(B)', 'HEC199(C)', 'HEC200(C)', 'HEC199(D)', 'HEC200(D)'], '1M7O': [], '1PBY': ['HEC991(A)', 'HEC992(A)', 'TRW43(C)'], '1PPJ': ['HEM501(C)', 'HEM502(C)', 'SMA2001(C)', 'ANY2002(C)', 'HEC501(D)', 'HEM501(P)', 'HEM502(P)', 'SMA3001(P)', 'ANY3002(P)', 'HEC501(Q)', 'FES501(E)', 'FES501(R)'], '1QKS': ['HEC601(A)', 'DHE602(A)', 'HEC601(B)', 'DHE602(B)'], '1QL3': ['HEC101(A)', 'HEC101(B)', 'HEC101(C)', 'HEC101(D)'], '1RZ5': ['HEC401(A)', 'HEC402(A)'], '1WVE': ['FAD599(A)-1', 'FAD599(A)-2', 'HEM699(C)', 'FAD599(B)-1', 'FAD599(B)-2', 'HEM699(D)'], '1YCC': ['HEC118(A)'], '2BH4': ['HEC1123(X)'], '2C1V': ['HEC401(A)', 'HEC402(A)', 'HEC401(B)', 'HEC402(B)'], '2C85': ['MN1567(A)'], '2FW5': ['HEC803(A)', 'HEC805(A)'], '2FYN': ['HEM501(A)', 'HEM502(A)', 'SMA503(A)', 'HEM301(B)', 'HEM501(D)', 'HEM502(D)', 'SMA503(D)', 'HEM301(E)', 'HEM501(G)', 'HEM502(G)', 'SMA503(G)', 'HEM301(H)', 'HEM501(J)', 'HEM502(J)', 'SMA503(J)', 'HEM301(K)', 'HEM501(M)', 'HEM502(M)', 'SMA503(M)', 'HEM301(N)', 'HEM501(P)', 'HEM502(P)', 'SMA503(P)', 'HEM301(Q)', 'FES200(C)', 'FES200(F)', 'FES200(I)', 'FES200(L)', 'FES200(O)', 'FES200(R)'], '2VHD': ['HEC401(A)', 'HEC402(A)', 'HEC401(B)', 'HEC402(B)'], '2XIL': ['HEM1301(A)'], '2XTS': ['MTE500(A)', 'HEC500(B)', 'MTE500(C)', 'HEC500(D)', 'CO1431(A)', 'CO1431(C)'], '351C': ['HEM0(A)'], '3A9F': ['HEC207(A)', 'HEC207(B)'], '3AYF': ['HEM801(A)', 'HEM802(A)'], '3HQ6': ['HEM400(A)', 'HEM401(A)', 'HEM400(B)', 'HEM401(B)'], '3MK7': ['HEM501(A)', 'HEM502(A)', 'HEC211(B)', 'HEC321(C)', 'HEC322(C)', 'HEM501(D)', 'HEM502(D)', 'HEC211(E)', 'HEC321(F)', 'HEC322(F)', 'HEM501(G)', 'HEM502(G)', 'HEC211(H)', 'HEC321(I)', 'HEC322(I)', 'HEM501(K)', 'HEM502(K)', 'HEC211(L)', 'HEC321(M)', 'HEC322(M)', 'CU503(A)', 'CU503(D)', 'CU503(G)', 'CU503(K)'], '3O5C': ['HEM401(A)', 'HEM501(A)', 'HEM402(B)', 'HEM502(B)', 'HEM403(C)', 'HEM503(C)', 'HEM404(D)', 'HEM504(D)'], '3VRD': ['HEC201(A)', 'HEC202(A)', 'FAD501(B)-1', 'FAD501(B)-2'], '3WFB': ['HEM801(B)', 'HEM802(B)', 'HEC201(C)', 'FE803(B)'], '4B2N': ['HEC700(A)', 'HEC701(A)', 'HEC700(B)', 'HEC701(B)'], '4D6T': ['HEM501(C)', 'HEM502(C)', '4X9503(C)', 'HEC501(D)', 'HEM501(P)', 'HEM502(P)', '4X9503(P)', 'HEC501(Q)', 'FES501(R)'], '4FA1': ['HEC402(A)', 'HEC403(A)', 'HEC402(B)', 'HEC403(B)', 'TRQ57(C)', 'TRQ57(E)'], '5Z25': ['HEC101(A)'], '6QKN': ['HEC401(A)', 'HEC402(A)', 'HEC401(B)', 'HEC402(B)'], '6V59': ['HEC601(A)', 'HEC602(A)', 'HEC601(B)', 'HEC602(B)']} 
pg.process_emaps(chains=chains, eta_moieties=included_eta_moieties, include_residues=include_residues, **process_kwargs) 
# Step 3: Generate graph database 
substitutions = [] 
edge_thresholds = [12] 
pg.generate_graph_database(sub=substitutions,edge_thresh=edge_thresholds) 
# Step 4: Mine for subgraphs 
pg.run_gspan(42,4,6) 
# Select subgraph pattern 
sg = pg.subgraph_patterns['1_#CHC_43'] 
# Final step: Find protein subgraphs 
sg.find_protein_subgraphs('structural')
# To switch clustering, call 'sg.set_clustering('sequence')
# Visualize pattern
sg.subgraph_to_Image().show()
# Visualize protein subgraph
sg.subgraph_to_Image(id='2FYN_6').show()
# Print reports
print(pg.mining_report())
print(sg.full_report())
