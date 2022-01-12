Classification
=================

Introduction
-------------
The efficiency and descriptive power of graph mining is enhanced when the algorithms are 
able to distinguish between different types of nodes and edges. Graph mining in PyeMap relies 
on each node and edge in the graph database being assigned a numerical label which corresponds to its category. In the context of protein graph mining with PyeMap, there is already a natural classification of nodes based on their amino acid residue type, but this categorization is by no means definitive for all purposes. As such, PyeMap allows users some control over how nodes/edges are classified in order to broaden 
or narrow the mining search space. 

Nodes
------

By default, each standard amino acid residue receives its own category (and thus its own numerical label), and all non-standard residues 
included in the analysis are labeled as 'NP' for non-protein (processed internally as '#'), which has its own numerical label. Currently, the customization options are limited, but we allow users to assign residue types their own category, or to labeled as 'X' for unknown residue 
type. This allows for some subsitutition within identified subgraph patterns. This can be done by passing the :py:attr:`node_categories` argument 
to the :py:func:`PDBGroup.generate_graph_database` function. :py:attr:`node_categories` should be formatted as a list of 1-character amino acid codes, and residues which are not included in this list will be classified as the unknown residue type 'X' for mining purposes.

**Example**

Set TRP and TYR to be their own category, all other amino acids will be labeled as 'X'.

.. code-block:: python

    node_categories = ['W','Y']
    pg.generate_graph_database(node_categories=node_categories)


Edges
------

By default, all edges are assigned the same numerical label of 1. One can classify edges based on their weights by passing the :py:attr:`edge_thresholds` argument to the :py:func:`PDBGroup.generate_graph_database` function. :py:attr:`edge_thresholds` should be formatted 
as a list of floats in ascending order, where each value indicates a cutoff threshold for an edge category. 

**Example**

4 types of edges: 

* weight < 8
* 8 <= weight <10
* 10 <= weight <12
* weight >= 12

.. code-block:: python

    edge_thresholds = [8,10,12]
    pg.generate_graph_database(edge_thresholds=edge_thresholds)




