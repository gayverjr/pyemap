Surface exposed residues
==============================

Introduction
-------------

Electron (or hole) transfer often proceeds from/to surface residues
to/from an acceptor/donor inside the protein. Therefore, identification of surface-exposed
residues is a key step for prediction of relevant electron/hole transfer pathways. 
Users can select one of two parameters to classify residues as buried/exposed: 
residue depth and relative solvent accessibility. In the graph images, buried residues will appear as ovals, 
while exposed residues will appear as rectangles.

Residue depth
-------------

Residue depth is a measure of solvent exposure that describes the extent to which a 
residue is buried within the protein structure. The parameter was first introduced by 
Chakravarty [Chakravarty1999]_ and coworkers, and is computed in PyeMap using the freely available program MSMS. MSMS computes a solvent-excluded surface
by rolling a probe sphere along the surface of the protein, which is represented as 
atomic spheres. The boundary of the volume reachable by the probe is taken to be the 
solvent-excluded surface The residue depth for each residue is calculated as the 
average distance of its respective atoms from the solvent-excluded surface [Sanner1996]_. In PyeMap, 
the threshold for classifying residues as buried/exposed is:

.. math::
   \mathbf{RD \leq}~\mathbf{3.03Å}

which is the threshold proposed by Tan and coworkers [Tan2009]_. Residues
3.03 Å and shallower will be classified as exposed in the final graph;
those deeper will be classified as buried.

Relative Solvent Accessibility
-------------------------------
Accessible surface area is a measure of solvent exposure, first introduced by 
Lee and Richards, which describes the surface area of a biomolecule that is accessible 
to solvent molecules [Lee1971]_. To calculate the accessible surface of each atom, a water sphere is 
rolled along the surface of the protein, making the maximum permitted van der Waals 
contacts without penetrating neighboring atoms [Shrake1973]_. The total accessible surface area for a
residue is the sum of the solvent accessible surface areas of its respective atoms.
In order to develop a threshold to classify residues as buried or exposed,
calculated ASA values need to be normalized based on corresponding reference values for a 
given residue. This requires precomputed or predefined maximal accessible surface area 
(MaxASA) for all residues. MaxASA is the maximal possible solvent accessible surface area
for a given residue. MaxASA values are obtained from theoretical calculations of Gly-X-Gly
tripeptides in water, where X is the residue of interest. From ASA and MaxASA, the relative
solvent accessibility (RSA) can be calculated by the formula:

.. math::
   RSA=\frac{ASA}{Max ASA}

Several scales for MaxASA have been published. PyeMap uses the most
recent theoretical scale from Tien and coworkers [Tien2013]_.
Relative solvent accessibility is calculated using the DSSP program developed by Kabsch and Sander [Sander1983]_.
In PyeMap, the RSA threshold chosen for exposed residues is:

.. math::
   \mathbf{RSA \geq 0.05}

as recommended by Tien and coworkers. Residues with RSA greater
than equal to 0.05 will be classified as exposed, those with lower RSA
values will be classified as buried.

Source
------

.. autosummary::
   :toctree: autosummary

   pyemap.process_data.calculate_residue_depth
   pyemap.process_data.calculate_asa

