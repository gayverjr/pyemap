Shortest paths
=========================================================

Introduction
-------------
Once you have a processed :ref:`emap <emap>` object, you can search for electron/hole transfer pathways. 
There are two modes: `source only`, or `specified target`. 

Source only
------------
When only a source is specified, `Dijkstra’s algorithm`_ is used to find the
shortest path between the specified source and every surface-exposed
residue in the graph. The pathways are classified into branches based
on the first surface-exposed residue reached along the pathway. Within
a branch, the pathways are ranked according to their score.

For an example of how the branches are structured, refer to the example below. In
this example, the flavin group on FAD (FAD510(A)-2) is specified as the source, and W356(A), W436(A), and ANP511(A) have been 
identified as surface-exposed residues. For W436(A), the shortest path involves first
going through W356(A), which itself is a surface-exposed residue. Thus
this path is assigned to branch W356(A), and given the ID 1b, as it is the
second shortest path in this branch. The shortest path to ANP511(A) is given
the ID 2a, and belongs to a different branch.

**Example 1: Source only**

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

.. image:: ../analysis/images/source_only.png


Specified target
-----------------
When a source and a target are specified, a NetworkX_ procedure based
on `Yen’s algorithm`_ is used to find the 5 shortest paths from source to
target. Yen’s algorithm exploits the idea that shortest paths are likely
to share common steps, and is able to compute the k shortest paths
between nodes in a graph with non-negative weights. The procedure
first finds the shortest path, then finds the next 4 shortest deviations.
The pathways are ranked according to their score. See the example below.

**Example 2: Specified Target:**

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

.. image:: ../analysis/images/target.png

.. _Yen’s algorithm: https://en.wikipedia.org/wiki/Yen%27s_algorithm
.. _Dijkstra’s algorithm: https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm
.. _NetworkX: https://networkx.github.io/

Source
-------
.. autosummary::
   :toctree: autosummary

   pyemap.shortest_paths.dijkstras_shortest_paths
   pyemap.shortest_paths.yens_shortest_paths