Installation
=========================================================

.. rubric:: Conda(recommended)

In order to set up the pre-requisites for PyeMap, you need to install the conda_ package manager. The easiest ways to get conda are
through the Anaconda_ and Miniconda_ distributions.

.. _conda: https://docs.conda.io/en/latest/

.. _Anaconda: https://www.anaconda.com/

.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html

Once you have a working copy of conda, create and activate a new virtual environment::

    $ conda create -n pyemap_env python=3.7
    $ conda activate pyemap_env

Now add conda channels in order to download dependencies::

    $ conda config --add channels conda-forge --add channels salilab --add channels bioconda --add channels gayverjr
    $ conda update --all

And finally, install pyemap::

    $ conda install pyemap 
   
.. rubric:: Pip

Pip installation will only install python dependencies, which is sufficient to run PyeMap analysis, but will be missing some features such as surface exposure and visualization::

    $ pip install --extra-index-url https://testpypi.python.org/pypi pyemap

For full functionality, you can download and install  MSMS_, DSSP_, and Graphviz_ separately.    

.. rubric:: Graphviz

PyeMap uses the Graphviz_ software to visualize the graphs. For graphs with <200 vertices, we use the `neato` program, 
which works by minimizing a global energy function. Within the neato program we have found that an experimental mode called `ipsep` 
(which you can read more about here_) provides the nicest looking graphs. The versions of graphviz available on conda unfortunately
do not come with ipsep enabled. To get around this, we suggest building graphviz from source, and adding a compiler argument which 
enables ipsep. Here's how to do it:

.. _here: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.591.840&rep=rep1&type=pdf
.. _MSMS: http://mgltools.scripps.edu/packages/MSMS
.. _DSSP: https://github.com/cmbi/xssp/releases
.. _Graphviz: https://graphviz.gitlab.io/

First, remove the graphviz that came with pygraphviz from your conda environment::

   $ conda activate pyemap_env
   $ conda remove graphviz --force-remove

Then, download the latest graphviz and compile from source::

   $ ./configure --with-ipsepcola=yes
   $ make
   $ make install

We emphasize that the safest way to do this is inside a newly prepared virtual environment. It is not recommended to try this within your base conda environment.





