Algorithms
=========================================================

Protein graph mining with PyeMap involves 4 steps. 

#. Generation of the protein graphs for each PDB in the analysis. This is done identically to the single protein anaylsis, and won't be discussed here.
#. Classification of the nodes and edges of each protein graph in order to generate a graph database.
#. Mining the graph database for shared patterns. There are two different types of mining algorithms available in PyeMap.
#. Graph matching to identify protein subgraphs, and clustering them based on similarity.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   classification
   mining
   clustering
