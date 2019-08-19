.. pyemap documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Overview of pyemap
=========================================================
pyemap is a python package aimed at automatic identification of electron and hole transfer
pathways in protein. The analysis is based on a coarse-grained version of Beratan and 
Onuchic’s Pathway model, and only accounts for the through-space hopping between
aromatic residues side chains. Side chains of aromatic residues and non-protein electron 
transfer active moieties are modeled as vertices in a weighted graph, where the edge 
weights are modified distancedependent penalty functions. Shortest path algorithms are 
used to compute the shortest pathways from a specified electron or hole donor to the 
surface of the protein, or to a user-specified acceptor. 


pyemap analysis happens in 3 steps. The first step is parsing a PDB or CIF file provided
by the user or fetched from the RCSB_ database. The next step is to construct the graph 
theory model of the protein crystal structure based on user specifications. 
Finally, the shortest paths between a specified electron/hole source to the surface or 
to a specified electron/hole acceptor are calculated. 

.. _RCSB: http://www.rcsb.org/

pyemap is intended to be used as the backend for the web application 
eMap_, and as a standalone python package.

.. _eMap: http://emap.bu.edu/

Audience
--------
The aim of this software is to efficiently identify possible electron hopping channels to 
be investigated further in quantitative and experimental studies. As such, our audience
includes computational and experimental chemists, biologists, and physicists interested
in gaining insight into potentially relevant electron/hole transfer pathways in proteins.



.. toctree::
   :maxdepth: 1
   :caption: Contents:

   reference
   tutorial
 


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
