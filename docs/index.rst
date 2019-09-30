.. pyemap documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Overview of PyeMap
=========================================================
.. image:: https://travis-ci.org/gayverjr/pyemap.svg?branch=master
   :target: https://travis-ci.org/gayverjr/pyemap
.. image:: https://codecov.io/gh/gayverjr/pyemap/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/gayverjr/pyemap/branch/master

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


PyeMap analysis happens in 3 steps. The first step is parsing a PDB or CIF file provided
by the user or fetched from the RCSB_ database. The next step is constructing the
graph theory model of the protein crystal structure. Finally, the shortest paths between a
specified electron/hole source to the surface or to a specified electron/hole acceptor are calculated.

.. _RCSB: http://www.rcsb.org/

PyeMap is intended to be used as the backend for the web application
eMap_, and as a standalone python package.

.. _eMap: http://emap.bu.edu/

Current Features
----------------
* Accepts valid .pdb or .cif structures provided by the user or fetched from RCSB_ database
* Tested on structures as large as 5350 residues (51599 atoms)
* Automatic detection of non-protein aromatic moieties such as porphyrins, nucleobases, and other aromatic cofactors
* Automatic detection of 60+ inorganic clusters such as iron-sulfur clusters and others
* Automatic identification of surface exposed residues using residue depth or solvent accessibility criteria
* Visualization of chemical structures and graphs
* User specified custom fragments
* Control over various parameters which determine connectivity of graph theory model
* Identification of most probable electron/hole transfer pathways from a specified donor to the protein surface or a specified electron/hole acceptor


Planned Features
----------------
* Tools for screening families of proteins for common electron/hole transfer pathways
* Generalization to DNA, protein-DNA, and other relevant biomolecules
* Improving the physical model of electron transfer by incorporating information on geometry-dependent electronic couplings and chemical-dependent energetics

Audience
--------
The aim of this software is to efficiently identify possible electron hopping channels to
be investigated further in quantitative and experimental studies. As such, our audience
includes computational and experimental chemists, biologists, and physicists interested
in gaining insight into potentially relevant electron/hole transfer pathways in proteins.

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
