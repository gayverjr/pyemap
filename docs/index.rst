.. pyemap documentation main file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Overview of PyeMap
=========================================================
.. image:: https://github.com/gayverjr/pyemap/workflows/ubuntu/badge.svg
   :target: https://github.com/gayverjr/pyemap/actions
.. image:: https://codecov.io/gh/gayverjr/pyemap/branch/main/graph/badge.svg
  :target: https://codecov.io/gh/gayverjr/pyemap/branch/main

- Website:  https://emap.bu.edu
- News:     https://twitter.com/eMap_protein

PyeMap is a python package aimed at automatic identification of electron and hole transfer
pathways in proteins. The analysis is based on a coarse-grained version of Beratan and
Onuchic’s Pathway model, and only accounts for through-space hopping between
aromatic residues side chains [Beratan1992]_. Side chains of aromatic residues and non-protein electron
transfer active moieties are modeled as vertices in a weighted graph, where the edge
weights are modified distance dependent penalty functions. Shortest path algorithms are
used to compute the shortest pathways from a specified electron or hole donor to the
surface of the protein, or to a user-specified acceptor.


PyeMap analysis is organized into 3 steps. The first step is parsing a PDB or CIF file provided
by the user or fetched from the RCSB_ database. The next step is constructing the
graph theory model of the protein crystal structure. Finally, the shortest paths between a
specified electron/hole source to the surface or to a specified electron/hole acceptor are calculated.

.. _RCSB: http://www.rcsb.org/

PyeMap serves as the backend for the web application eMap_, 
and can also be used as a fully functional Python package.

.. _eMap: http://emap.bu.edu/

Current Features
----------------
* Identification of most probable electron/hole transfer pathways from a specified donor to the protein surface or a specified electron/hole acceptor
* Accepts valid .pdb or .cif structures provided by the user or fetched from RCSB_ database
* Automatic detection of non-protein aromatic moieties such as porphyrins, nucleobases, and other aromatic cofactors
* Automatic detection of 60+ inorganic clusters such as iron-sulfur clusters and others
* User specified custom fragments
* Visualization of chemical structures and graphs
* Automatic identification of surface exposed residues using residue depth or solvent accessibility criteria
* Control over various parameters which determine connectivity of graph theory model
* Tested on structures as large as 5350 residues (51599 atoms)



In Development
----------------
* Tools for screening families of proteins for common electron/hole transfer pathways
* Generalization to DNA, protein-DNA complexes etc.
* Improving the physical model of electron transfer by incorporating information on geometry-dependent electronic couplings and site sensitive energetics


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   install
   tutorial/tutorial
   reference
   bibliography
   cite
   credits


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
