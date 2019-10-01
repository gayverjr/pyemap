Installation
=========================================================
PyeMap officially supports Python versions 3.5 and later, and has been tested for Linux and OSX platforms.

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

    $ conda config --add channels conda-forge --add channels salilab --add channels bioconda --add channels gayverjr
    $ conda update --all

And finally, install PyeMap::

    $ conda install pyemap

**Pip**

Pip installation will only install python dependencies, and requires Graphviz_ in order to work.
This is sufficient to run PyeMap analysis and view graph images, but some features will be missing::

    $ pip install pyemap

For full functionality, install the following packages:

    - RDKit_: visualization of chemical stuctures
    - MSMS_: residue depth criterion for surface exposed residues
    - DSSP_: solvent accessibility criterion for surface exposed residues
    - wget_: fetching PDBs from RCSB_ database

All of these packages can be downloaded free of charge from their respective owners, and build recipes are available on the
Anaconda_ cloud for some platforms.

**Graphviz**

PyeMap uses the Graphviz_ software to visualize the graphs. For graphs with <200 vertices, we use the `neato` program,
which works by minimizing a global energy function. Within the neato program we have found that an experimental mode called `ipsep`
(which you can read more about here_) provides the cleanest looking graphs. The versions of Graphviz distributed through conda and homebrew
do not come with ipsep enabled. To get around this, we suggest building graphviz from source, and adding a compiler argument which
enables ipsep. Here's how to do it:

.. _here: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.591.840&rep=rep1&type=pdf
.. _MSMS: http://mgltools.scripps.edu/packages/MSMS
.. _DSSP: https://github.com/cmbi/xssp/releases
.. _Graphviz: https://graphviz.gitlab.io/
.. _RDKit: https://www.rdkit.org/docs/Install.html
.. _wget: https://www.gnu.org/software/wget/
.. _RCSB: https://www.rcsb.org/

First, remove any prior installations of Graphviz from your conda environment::

   $ conda activate pyemap_env
   $ conda remove graphviz --force-remove

Then, download the latest Graphviz and compile from source::

   $ ./configure --with-ipsepcola=yes
   $ make
   $ make install

And then add the executables to a directory on your systems’ path.
