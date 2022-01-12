Graph construction
=========================================================

Introduction
-------------
The first step to constructing the graph theory model is constructing a pairwise distance matrix for the selected amino acid residues.
The distance is calculated either between centers of mass of the side chains, or between their
closest atoms. For standard protein residues, only side chain atoms
are considered in the calculation. All atoms of automatically identified
non-protein ET active moieties and user-specified custom fragments
are considered in the distance calculations. From the distance matrix,
an undirected weighted graph is constructed using NetworkX_, with
the calculated distances as weights. 

One of two algorithms are then used to prune the edges of the graph, which is specified by the `edge_prune` keyword argument.

**Percent-based algorithm**

The following thresholds are then imposed on the graph:

* Only the shortest `percent_edges`% of edges or 2 edges, whichever is greater, are kept for each vertex
* Of those edges, all edges > `distance_cutoff` are removed
* Of the remaining edges, only those with length :math:`l \leq \overline{l}_{vertex} + n\sigma_{vertex}` are kept, where :math:`σ_{vertex}` is the standard deviation in length in the remaining set of edges for given vertex, and n = `num_std_dev_edges` 
* All disconnected vertices are removed

Specify `edge_prune` as 1 to use this algorithm.

**Degree-based algorithm**

First, all edges > `distance_cutoff` are removed. Then, for each node, the edges are sorted by weight, and all edges except the smallest 
`max_degree` are marked as candidates for removal. We then sort this group of edges by weight, and iterate over all edges starting from the largest edge. 
At each iteration, we remove the edge from the graph if the degree of both associated is larger than `max_degree`.

Specify `edge_prune` as 0 to use this algorithm.

**Penalty Functions**

After the edges are pruned, the weights are recast as modified distance dependent penalty functions:

.. math::
   P'=-log_{10}(\epsilon)

where: 

.. math::
   \epsilon = \alpha \exp(-\beta(R-R_{offset}))

α, β, and :math:`R_{offset}` are hopping parameters, similar to the through-space
tunneling penalty function in the Pathways model [Beratan1992]_. All subsequent
calculations are performed using the modified penalty functions as
edge weights. When using default hopping parameters (α = 1.0,
β = 2.3, Roffset = 0.0), the edge weights will be equal to the distances
(multiplied by a prefactor of :math:`2.3*log_{10}(e)` ≈ 1).

Distance thresholds and penalty function parameters can be modified at the process step. 

Visualization and further analysis
-----------------------------------
The graph can be interacted with and written to file using the :class:`~pyemap.emap` object. The graph is visualized using PyGraphviz_ and 
Graphviz_. The graph is stored as a :class:`networkx.Graph` object in the `init_graph` and `paths_graph` attributes of the :class:`~pyemap.emap` object.

.. _PyGraphviz: https://pygraphviz.github.io/
.. _Graphviz: http://www.graphviz.org/
.. _NetworkX: https://networkx.github.io/

	>>> G = my_emap.init_graph
	>>> print(G.edges[('W17(A)', 'W45(A)')]['distance'])
	>>> 12.783579099370808

Source
-------

.. autosummary::
   :toctree: autosummary

   pyemap.process_data.create_graph
   pyemap.process_data.pathways_model
   pyemap.process_data.closest_atom_dmatrix
   pyemap.process_data.com_dmatrix