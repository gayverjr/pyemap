Installation
=========================================================
pyemap requires python version 3.5 or later. Ensure that you have working installations of git_ and pip_ on your machine.

While not a requirement, for nicer looking graphs, please see our note below regarding graphviz.

.. _git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

.. _pip: https://pypi.org/project/pip/


Setting up conda environment
-----------------------------
In order to set up the pre-requisites for pyemap, you need to install the conda_ package manager. The easiest ways to get conda are
through the Anaconda_ and Miniconda_ distributions.

.. _conda: https://docs.conda.io/en/latest/

.. _Anaconda: https://www.anaconda.com/

.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html

Once you have a working copy of conda, create and activate a new virtual environment::

    $ conda create -n pyemap_env
    $ conda activate pyemap_env

And install the necessary dependencies::

    $ conda install -c conda-forge pygraphviz rdkit 
    $ conda install numpy scipy biopython networkx pillow
   
Installing pyemap 
----------------------

First obtain the source code from Github::

   $ git clone https://github.com/gayverjr/pyemap.git

And then install the pyemap package with pip::

   $ pip install -e pyemap

To test your installation, try running the example in pyemap/examples::

   $ cd pyemap/examples
   $ python example.py

Graphviz
---------
pyemap uses the graphviz software to visualize the graphs. For graphs with <200 vertices, we use the `neato` program, 
which works by minimizing a global energy function. Within the neato program we have found that an experimental mode called `ipsep` 
(which you can read more about here) provides the nicest looking graphs. The versions of graphviz available on conda unfortunately
do not come with ipsep enabled. To get around this, we suggest building graphviz from source, and adding a compiler argument which 
enables ipsep. Here's how to do it:

First, remove the graphviz that came with pygraphviz from your conda environment::

   $ conda activate pyemap_env
   $ conda remove graphviz --force-remove

Then, download the latest graphviz and compile from source::

   $ ./configure --with-ipsepcola=yes
   $ make
   $ make install


For contributors
------------------
pyemap serves as the backend for the web application eMap, so any contributions need to be tested for compatibility with the web version. 
eMap uses flask as its framework. To install the necessary dependencies, do:

   $ pip install flask
   $ pip install flask-restful
   $ pip install flask-table

