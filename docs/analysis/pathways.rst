==============================
Find Shortest Pathways
==============================

Introduction
-------------
Once you have a processed emap object which contains the graph theory model of the protein structure, you can search for pathways.
There are two modes of search: source only and specified target. Pathways are stored as :class:`~pyemap.ShortestPath` objects, which are
organized into :class:`~pyemap.Branch` objects. The data is stored in the :class:`~pyemap.emap` object which was passed in. For a report of 
pathways found by pyemap, use :func:`~pyemap.emap.report`.

Source only
------------
When only a source node is selected, a NetworkX implementation of `Dijkstra’s algorithm`_ is used to
calculate the shortest path from the source to each surface-exposed
residue. In the output, the pathways are organized into "branches"
based on the first surface-exposed residue reached during the course
of the pathway. The source node will be colored yellow in the graph visualization.

Specified target
-----------------
If a target is specified, a NetworkX procedure based on `Yen’s algorithm`_ is used to 
calculate the shortest paths from source to target. The target does not
need to be a surface-exposed residue. The target node will be colored blue in the graph visualization.

.. _Yen’s algorithm: https://en.wikipedia.org/wiki/Yen%27s_algorithm
.. _Dijkstra’s algorithm: https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm

Examples
--------
**Source only:**

   >>> import pyemap
   >>> my_emap = pyemap.fetch_and_parse("1u3d")
   >>> pyemap.process(my_emap)
   >>> pyemap.find_paths(my_emap,"FAD510(A)-2")
   >>> print(my_emap.report())
   Branch: W356(A)
   1a: ['FAD510(A)-2', 'W356(A)'] 9.48
   1b: ['FAD510(A)-2', 'W356(A)', 'Y432(A)', 'W436(A)'] 26.44
   1c: ['FAD510(A)-2', 'W356(A)', 'W213(A)', 'W62(A)', 'W217(A)'] 34.72
   Branch: ANP511(A)
   2a: ['FAD510(A)-2', 'FAD510(A)-1', 'ANP511(A)'] 14.15
   ...

   >>> my_emap.paths_graph_to_Image().show()

.. image:: images/source_only.png

**Specified Target:**

   >>> pyemap.find_paths(my_emap,"FAD510(A)-2", target = "W324(A)", max_paths=10)
   >>> my_emap.report()
   Branch: W324(A)
   1a: ['FAD510(A)-2', 'W400(A)', 'W377(A)', 'W324(A)'] 24.15
   1b: ['FAD510(A)-2', 'W385(A)', 'Y53(A)', 'W377(A)', 'W324(A)'] 35.25
   1c: ['FAD510(A)-2', 'W400(A)', 'W334(A)', 'W379(A)', 'W324(A)'] 36.37
   1d: ['FAD510(A)-2', 'W400(A)', 'W377(A)', 'W492(A)', 'W324(A)'] 49.20
   1e: ['FAD510(A)-2', 'W385(A)', 'Y53(A)', 'Y309(A)', 'W377(A)', 'W324(A)'] 50.67
   ...

   >>> my_emap.paths_graph_to_Image().show()

.. image:: images/target.png


Find Paths
------------
.. automodule:: pyemap.pathway_analysis
   :members: find_paths

