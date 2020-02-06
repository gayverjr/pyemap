<div align="center">
  <img src="https://github.com/gayverjr/pyemap/blob/master/docs/logo/pyemap_logo.png">
</div>

[![GitHub release](https://img.shields.io/github/release/gayverjr/pyemap.svg)](https://github.com/gayverjr/pyemap)[![Build Status](https://travis-ci.org/gayverjr/pyemap.svg?branch=master)](https://travis-ci.org/gayverjr/pyemap) [![codecov](https://codecov.io/gh/gayverjr/pyemap/branch/master/graph/badge.svg)](https://codecov.io/gh/gayverjr/pyemap/branch/master) [![Documentation Status](https://readthedocs.org/projects/pyemap/badge/?version=latest)](https://pyemap.readthedocs.io/en/latest/?badge=latest) [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/gayverjr/pyemap/blob/master/LICENSE)

[![Website emap.bu.edu](https://img.shields.io/website-up-down-green-red/https/emap.bu.edu.svg)](https://emap.bu.edu/) [![DOI:10.1021/acs.jpcb.9b04816](https://zenodo.org/badge/DOI/10.1021/acs.jpcb.9b04816.svg)](https://doi.org/10.1021/acs.jpcb.9b04816) [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/gayverjr/pyemap.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/gayverjr/pyemap/context:python)

PyeMap is a Python package aimed at automatic identification of electron and hole transfer pathways in proteins.
It serves as the backend for the web application [eMap](https://emap.bu.edu/), and can also be used as a fully functional Python package.

- **Documentation:** https://readthedocs.org/projects/pyemap/
- **Website:** https://emap.bu.edu
- **News:** https://twitter.com/eMap_protein

# Installation
PyeMap officially supports Python versions 3.5 and later, and has been tested for Linux and OSX platforms. Below is an abbreviated version of the instructions provided in the [documentation](https://pyemap.readthedocs.io/en/latest/install.html).
### Conda (recommended):
The conda recipe will install all dependencies necessary for full functionality.
```
# create new virtual environment
$ conda create -n pyemap_env python=3.7
$ conda activate pyemap_env
# include channels for dependencies, only needs to be done once
$ conda config --add channels conda-forge --add channels salilab --add channels bioconda --add channels gayverjr
$ conda update --all
# install pyemap
$ conda install pyemap
```

### Pip
Pip installation will only install python dependencies, and requires [Graphviz](https://graphviz.gitlab.io/) in order to work. This is sufficient to run PyeMap analysis and view graph images, but some features will be missing.
```
pip install pyemap
```
For full functionality, install [RDKit](https://www.rdkit.org/docs/Install.html), [MSMS](http://mgltools.scripps.edu/packages/MSMS), [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html), and [wget](https://www.gnu.org/software/wget/), all of which can be downloaded free of charge from their owners, and are available on most platforms.

# Getting started
Please see our [tutorial](https://pyemap.readthedocs.io/en/latest/tutorial/tutorial.html) on how to use PyeMap.

# Bugs and feature requests
Please report any bugs and make feature requests [here](https://github.com/gayverjr/pyemap/issues). For issues exclusive to the web version [eMap](https:emap.bu.edu), please send an email to <emap.bu@gmail.com>. We also greatly encourage users to contribute by making pull requests on [GitHub](https://github.com/gayverjr/pyemap)!

# License
Released under the 3-Clause BSD license (see LICENSE).

Copyright (C) 2019

James Gayvert <jrg444@gmail.com>

Ruslan Tazhigulov <ruslan.tazhigulov@gmail.com>

Ksenia Bravaya <bravaya@bu.edu>
