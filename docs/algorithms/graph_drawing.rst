Graph construction
=========================================================

Introduction
-------------
The first step to constructing the graph theory model is constructing a pairwise distance matrix is built for the selected  The distance is calculated
either between centers of mass of the side chains, or between their
closest atoms. For standard protein residues, only side chain atoms
are considered in the calculation. All atoms of automatically identified
non-protein ET active moieties and user-specified custom fragments
are considered in the distance calculations. From the distance matrix,
an undirected weighted graph is constructed using NetworkX_, with
the calculated distances as weights. The following thresholds are then
imposed on the graph:

* Only the shortest 1% of edges or 2 edges, whichever is greater, are kept for each vertex
* Of those edges, all edges > 20 Å are removed
* Of the remaining edges, only those with length :math:`l \leq \overline{l}_{vertex} + n\sigma_{vertex}` are kept, where :math:`σ_{vertex}` is the standard deviation in length in the remaining set of edges for given vertex, and n = 1 by default
* All disconnected vertices are removed

For the remaining edges, the weights are recast as modified distance dependent penalty functions:

.. math::
   P'=-log_{10}(\epsilon)

where: 

.. math::
   \epsilon = \alpha \exp(-\beta(R-R_{offset}))

α, β, and :math:`R_{offset}` are hopping parameters, similar to the through-space
tunneling penalty function in the Pathways model. All subsequent
calculations are performed using the modified penalty functions as
edge weights. When using default hopping parameters (α = 1.0,
β = 2.3, Roffset = 0.0), the edge weights will be equal to the distances
(multiplied by a prefactor of :math:`2.3*log_{10}(e)` ≈ 1).

Distance thresholds and penalty function parameters can be modified at the process step. 

Visualization
--------------
The graph can be interacted with and written to file using the :ref:`emap <emap>` object. The graph is visualized using PyGraphviz_ and 
Graphviz_.

.. _PyGraphviz: https://pygraphviz.github.io/
.. _Graphviz: http://www.graphviz.org/
.. _NetworkX: https://networkx.github.io/


Source
-------

.. autosummary::
   :toctree: autosummary

   pyemap.process_data.create_graph
   pyemap.process_data.pathways_model
   pyemap.process_data.closest_atom_dmatrix
   pyemap.process_data.com_dmatrix