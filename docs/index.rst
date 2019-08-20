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


Citation
--------

Primary citation:|br| 
Tazhigulov, R. N., Gayvert, J. R., Wei, M., & Bravaya, K. B. (2019). 
eMap: A Web Application for Identifying and Visualizing Electron or Hole Hopping Pathways
in Proteins. *The Journal of Physical Chemistry B*. https://doi.org/10.1021/acs.jpcb.9b04816

In addition to citing eMap, please cite the third party software we depend on:

Biopython_: |br|
Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff 
F, Wilczynski B and de Hoon MJL (2009) *Bioinformatics*, 25, 1422-1423. |br|
Hamelryck T and Manderick B (2003) *Bioinformatics*, 22, 2308-2310.


MSMS_ (Residue Depth): |br|
Sanner, M. F.; Olson, A. J.; Spehner, J. C. *Biopolymers*, 1996, 38, 305-320.

DSSP_ (Solvent Accessibility): |br|
Touw, W. G.; Baakman, C.; Black, J.; te Beek, T. A.; Krieger, E.; Joosten, R. P.; Vriend, G. Nucleic Acids Res., 2015, 43, D364-D368. |br|
Wolfgang, K.; Christian, S. *Biopolymers*, 1983, 22, 2577-2637. 


.. _MSMS: http://mgltools.scripps.edu/packages/MSMS/

.. _DSSP: https://swift.cmbi.umcn.nl/gv/dssp/index.html

.. _Biopython: https://biopython.org/wiki/Documentation

.. |br| raw:: html

  <br/>

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   
   tutorial
   reference
   algorithms
 


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
