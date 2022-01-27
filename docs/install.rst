Installation
=========================================================
PyeMap officially supports Python versions 3.6 and later, and has been tested for Linux and OSX platforms.

**Conda (recommended)**

In order to set up the pre-requisites for PyeMap, you need to install the conda_ package manager. The easiest ways to get conda are
through the Anaconda_ or Miniconda_ distributions.

.. _conda: https://docs.conda.io/en/latest/

.. _Anaconda: https://www.anaconda.com/

.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html

Once you have a working copy of conda, create and activate a new virtual environment::

    $ conda create -n pyemap_env
    $ conda activate pyemap_env

Now add conda channels in order to download dependencies::

    $ conda config --add channels conda-forge --add channels gayverjr
    $ conda update --all

And finally, install PyeMap::

    $ conda install pyemap

Our conda recipe does not include MSMS_ or DSSP_, which are used to classify residues as surface 
exposed or buried. Please install these packages separately if you need that functionality.

**Pip**

Pip installation will only install python dependencies, and requires a working Graphviz_ installation.
This is sufficient to run PyeMap analysis, but some features will be missing::

    $ pip install pyemap

For full functionality, install the following packages:

    - RDKit_: visualization of chemical stuctures
    - MSMS_: residue depth criterion for surface exposed residues (not available on MacOS Catalina)
    - DSSP_: solvent accessibility criterion for surface exposed residues
    - MUSCLE_: Multiple sequence alignment

All of these packages can be downloaded free of charge from their respective owners, and build recipes are available on the
Anaconda_ cloud for some platforms.

.. _here: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.591.840&rep=rep1&type=pdf
.. _MSMS: http://mgltools.scripps.edu/packages/MSMS
.. _DSSP: https://github.com/cmbi/xssp/releases
.. _Graphviz: https://graphviz.gitlab.io/
.. _RDKit: https://www.rdkit.org/docs/Install.html
.. _MUSCLE: http://www.drive5.com/muscle/
.. _RCSB: https://www.rcsb.org/


