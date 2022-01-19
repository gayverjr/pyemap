<div align="center">
  <img src="https://github.com/gayverjr/pyemap/blob/main/docs/logo/pyemap_logo.png">
</div>

[![GitHub release](https://img.shields.io/github/release/gayverjr/pyemap.svg)](https://github.com/gayverjr/pyemap)[![Build Status](https://github.com/gayverjr/pyemap/workflows/ubuntu/badge.svg)](https://github.com/gayverjr/pyemap/actions) [![codecov](https://codecov.io/gh/gayverjr/pyemap/branch/main/graph/badge.svg)](https://codecov.io/gh/gayverjr/pyemap/branch/main) [![Documentation Status](https://readthedocs.org/projects/pyemap/badge/?version=latest)](https://pyemap.readthedocs.io/en/latest/?badge=latest) [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/gayverjr/pyemap/blob/main/LICENSE)

[![Website emap.bu.edu](https://img.shields.io/website-up-down-green-red/https/emap.bu.edu.svg)](https://emap.bu.edu/) [![DOI:10.1021/acs.jpcb.9b04816](https://zenodo.org/badge/DOI/10.1021/acs.jpcb.9b04816.svg)](https://doi.org/10.1021/acs.jpcb.9b04816) [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/gayverjr/pyemap.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/gayverjr/pyemap/context:python)

PyeMap is a Python package aimed at automatic identification of electron and hole transfer pathways in proteins.
It serves as the backend for the web application [eMap](https://emap.bu.edu/), and can also be used as a fully functional Python package.

- **Documentation:** https://readthedocs.org/projects/pyemap/
- **Website:** https://emap.bu.edu
- **News:** https://twitter.com/eMap_protein

# Installation
PyeMap officially supports Python versions 3.6, 3.7, and 3.8, and has been tested for Linux and OSX platforms. Below is an abbreviated version of the instructions provided in the [documentation](https://pyemap.readthedocs.io/en/latest/install.html).
### Conda (recommended):
```
# create new virtual environment
$ conda create -n pyemap_env python=3.7
$ conda activate pyemap_env
# include channels for dependencies, only needs to be done once
$ conda config --add channels conda-forge --add channels gayverjr
$ conda update --all
# install pyemap
$ conda install pyemap
```

Our conda recipe does not include [MSMS](http://mgltools.scripps.edu/packages/MSMS) (not available on Mac OS Catalina), or [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html), 
which are used for classifying residues as buried or surface exposed. Please install these packages separately if you need that functionality.

### Pip
Pip installation will only install python dependencies, and requires [Graphviz](https://graphviz.gitlab.io/) in order to work. This is sufficient to run PyeMap analysis and view graph images, but some features will be missing.
```
pip install pyemap
```
For full functionality, install [RDKit](https://www.rdkit.org/docs/Install.html), [MSMS](http://mgltools.scripps.edu/packages/MSMS) (not available on Mac OS Catalina), [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html), and [wget](https://www.gnu.org/software/wget/), all of which can be downloaded free of charge from their owners, and are available on most platforms.

# Getting started
Please see our [tutorial](https://pyemap.readthedocs.io/en/latest/tutorial/tutorial.html) on how to use PyeMap. This tutorial is also available as a Jupyter notebook in the examples directory.

# Bugs and feature requests
Please report any bugs and make feature requests [here](https://github.com/gayverjr/pyemap/issues). For issues exclusive to the web version [eMap](https:emap.bu.edu), please send an email to <emap.bu@gmail.com>. We also greatly encourage users to contribute by making pull requests on [GitHub](https://github.com/gayverjr/pyemap)!

# License
Released under the 3-Clause BSD license (see LICENSE).

Copyright (C) 2019

James Gayvert <jrg444@gmail.com>

Ruslan Tazhigulov <ruslan.tazhigulov@gmail.com>

Ksenia Bravaya <bravaya@bu.edu>
