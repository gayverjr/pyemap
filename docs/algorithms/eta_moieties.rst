Non-protein electron transfer active moieties
=============================================

Introduction
------------
Frequently, protein crystal structures contain residues which are not
amino acids, and do not belong to the polypeptide chain(s). Unless they
are solvent molecules or salt ions not belonging to any co-factors, they
can play significant role in electron/hole transfer. pyemap automatically
identifies those non-protein electron/hole transfer (ET) active moieties,
and gives users the option to include them in the analysis.
In the current implementation, non-protein ET active moieties identified by pyemap are 
non-amino acid aromatic sites or extended conjugated systems. For a given non-standard 
co-factor (e.g., flavin adenine dinucleotide), there can be multiple non-protein ET 
active moieties identified by pyemap, and they will appear as separate nodes on the graph
if selected for analysis.

Identification
--------------
After initial parsing, non-protein residues are analyzed for detection of
ET active moieties. For each non-standard residue, a chemical graph
is constructed using the NetworkX library, consisting of the O, C, N,
P and S atoms in the residue. To isolate the conjugated systems, an
edge is only drawn between two atoms j and k if:

.. math::
   r_{\text{jk}} \leq \ \ \overline{x} - 3\sigma_{\overline{x}}

where x is the mean single-bond distance between those two elements,
and σx is the standard deviation. The data was obtained from the
online CRC Handbook of Chemistry and Physics. If there are any conjugated systems, the resulting chemical graph will
be a forest of connected component subgraphs. Each subgraph that
contains a cycle, or consists of 10 or more atoms will be considered a
non-protein ET active moiety, and can be selected for the analysis.

Visualization
-------------
Chemical structures of non-protein electron transfer active moieties are available in pyemap as SMILES strings, and 
`NGL Viewer`_ selection strings. For a given residue, both can be accessed through the :ref:`emap <emap>` object. 


.. _NGL Viewer: http://nglviewer.org/ngl/api/

Source
------

.. toctree::
   :maxdepth: 1

   smiles
   custom_residues