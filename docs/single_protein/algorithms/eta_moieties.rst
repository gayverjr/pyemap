Non-protein electron transfer active moieties
=============================================

Introduction
------------
Frequently, protein crystal structures contain residues which are not
amino acids, and do not belong to the polypeptide chain(s). Unless they
are solvent molecules or salt ions not belonging to any co-factors, they
can play significant role in electron/hole transfer. PyeMap automatically
identifies those non-protein electron/hole transfer (ET) active moieties,
and gives users the option to include them in the analysis.
In the current implementation, non-protein ET active moieties identified by PyeMap are
non-amino acid aromatic sites, extended conjugated systems, and a pre-defined list of metal clusters. For a given non-standard
co-factor (e.g., flavin adenine dinucleotide), there can be multiple non-protein ET
active moieties identified by PyeMap, and they will appear as separate nodes on the graph
if selected for analysis.

Identification
--------------

**Aromatic moieties and extended conjugated chains**

After initial parsing, non-protein residues are analyzed for detection of
ET active moieties. For each non-standard residue, a chemical graph
is constructed using the NetworkX library, consisting of the O, C, N,
P and S atoms in the residue. To isolate the conjugated systems, an
edge is only drawn between two atoms j and k if:

.. math::
   r_{\text{jk}} \leq \ \ \overline{x} - 2\sigma_{\overline{x}}

where x is the mean single-bond distance between those two elements,
and :math:`\sigma_\bar{x}` is the standard deviation. If there are any conjugated systems, the resulting chemical graph will
be a forest of connected component subgraphs. Each subgraph that
contains a cycle, or consists of 10 or more atoms will be considered a
non-protein ET active moiety, and can be selected for the analysis.

**Clusters**

The PyeMap repository contains a list of 66 inorganic clusters which are automatically identified by their 3 character residue names. All atoms
in the residue are collected as part of the customized residue object, and a pre-rendered image is used for visualization of chemical structure.
Otherwise, they can be used and interacted with just like any other residue. The list of clusters and pre-rendered images were obtained from the 
Protein Data Bank in Europe (PDBe_).

Visualization
-------------
Chemical structures of residues(not including user-specified residues) can be visualized using the :func:`~pyemap.emap.residue_to_Image()`, :func:`~pyemap.emap.init_graph_to_Image()` functions.
SMARTS strings and `NGL Viewer`_ selection strings are also accessible through the :class:`~pyemap.emap` object. Note that SMARTS strings are
not available for clusters and user-specified residues.


.. _NGL Viewer: http://nglviewer.org/ngl/api/
.. _PDBe: https://www.ebi.ac.uk/pdbe/

Source
------

.. toctree::
   :maxdepth: 1

   custom_residues
   chemical_structures
