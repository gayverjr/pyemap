<div align="center">
  <img src="https://github.com/gayverjr/pyemap/blob/master/docs/logo/pyemap_logo.png">
</div>

[![Build Status](https://travis-ci.org/gayverjr/pyemap.svg?branch=master)](https://travis-ci.org/gayverjr/pyemap) [![codecov](https://codecov.io/gh/gayverjr/pyemap/branch/master/graph/badge.svg)](https://codecov.io/gh/gayverjr/pyemap/branch/master) [![Documentation Status](https://readthedocs.org/projects/pyemap/badge/?version=latest)](https://pyemap.readthedocs.io/en/latest/?badge=latest) [![Website emap.bu.edu](https://img.shields.io/website-up-down-green-red/https/emap.bu.edu.svg)](https://emap.bu.edu/) [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/gayverjr/pyemap/blob/master/LICENSE) [![DOI:10.1021/acs.jpcb.9b04816](https://zenodo.org/badge/DOI/10.1021/acs.jpcb.9b04816.svg)](https://doi.org/10.1021/acs.jpcb.9b04816)

PyeMap is a Python package aimed at automatic identification of electron and hole transfer pathways in proteins. It is intended to serve as the backend for the web application [eMap](https://emap.bu.edu), and as a standalone package.

- **Documentation:** https://readthedocs.org/projects/pyemap/
- **Website:** https://emap.bu.edu
- **News:** https://twitter.com/eMap_protein
# Installation
PyeMap requires Python versions 3.5 and later, and a working copy of [wget](https://www.gnu.org/software/wget/) in the path. It has been tested on OSX and Linux platforms. This is an abbreviated version of the instructions provided in our [documentation](https://pyemap.readthedocs.io/en/latest/install.html).

### Conda(recommended OSX and Linux):
For OSX and Linux users, the conda recipe will install all dependencies necessary for full functionality.
```
# create new virtual environment
$ conda create -n pyemap_env
$ conda activate pyemap_env
# include channels for dependencies, only needs to be done once
$ conda config --add channels conda-forge --add channels salilab --add channels bioconda --add channels gayverjr
$ conda update --all
# install pyemap
$ conda install pyemap
```

### Pip(all platforms)
Pip installation will only install python dependencies, which is sufficient to run PyeMap analysis, but will be missing some features.
```
pip install --extra-index-url https://testpypi.python.org/pypi pyemap
```
For full functionality you will need to install [PyGraphviz](https://pygraphviz.github.io/), [RDKit](https://www.rdkit.org/docs/Install.html), [MSMS](http://mgltools.scripps.edu/packages/MSMS), [MKDSSP](https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html), and [Graphviz](https://graphviz.gitlab.io/), all of which can be downloaded free of charge from their owners, and are available on most platforms. Windows users should refer to our [documentation](https://pyemap.readthedocs.io/en/latest/install.html) for more information.


# Bugs and feature requests
Please report any bugs and make feature requests [here](https://github.com/gayverjr/pyemap/issues). For issues exclusive to the web version [eMap](https:emap.bu.edu), please send an email to <emap.bu@gmail.com>. We also greatly encourage users to contribute by making pull requests on [GitHub](https://github.com/gayverjr/pyemap)!

# License
Released under the 3-Clause BSD license (see LICENSE).

Copyright (C) 2019

James Gayvert <jrg444@gmail.com>

Ruslan Tazhigulov <ruslan.tazhigulov@gmail.com>

Ksenia Bravaya <bravaya@bu.edu>
