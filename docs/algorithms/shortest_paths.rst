Shortest Paths
=========================================================

Introduction
-------------
Once you have a processed :ref:`emap <emap>` object, you can search for electron/hole transfer pathways. 
There are two modes: "source only", or "specified target". 

Source only
------------
When only a source is specified, Dijkstra’s algorithm is used to find the
shortest path between the specified source and every surface-exposed
residue in the graph. The pathways are classified into branches based
on the first surface-exposed residue reached along the pathway. Within
a branch, the pathways are ranked according to their score.

For an example of how the branches are structured, refer to the figure below. In
this example, W377, Y53, Y309, W492, and Y383 are surface-exposed
residues, and each of the pathways displayed is the shortest path to
that residue from W400. For W324, the shortest path involves first
going through W377, which itself is a surface-exposed residue. Thus
this path is assigned to branch W377, and given the ID 1b, as it is the
second shortest path in this branch. The shortest path to Y383 is given
the ID 2a, and belongs to a different branch.


Specified target
-----------------
When a source and a target are specified, a NetworkX procedure based
on Yen’s algorithm is used to find the 5 shortest paths from source to
target. Yen’s algorithm exploits the idea that shortest paths are likely
to share common steps, and is able to compute the k shortest paths
between nodes in a graph with non-negative weights. The procedure
first finds the shortest path, then finds the next 4 shortest deviations.
The pathways are ranked according to their score.

Source
-------
.. autosummary::
   :toctree: autosummary

   pyemap.shortest_paths.dijkstras_shortest_paths
   pyemap.shortest_paths.yens_shortest_paths