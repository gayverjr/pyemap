<div align="center">
  <img src="https://github.com/gayverjr/pyemap/blob/master/docs/logo/pyemap_logo.png">
</div>

[![Build Status](https://travis-ci.org/gayverjr/pyemap.svg?branch=master)](https://travis-ci.org/gayverjr/pyemap) [![codecov](https://codecov.io/gh/gayverjr/pyemap/branch/master/graph/badge.svg)](https://codecov.io/gh/gayverjr/pyemap/branch/master) [![Documentation Status](https://readthedocs.org/projects/pyemap/badge/?version=latest)](https://pyemap.readthedocs.io/en/latest/?badge=latest) [![Website emap.bu.edu](https://img.shields.io/website-up-down-green-red/https/emap.bu.edu.svg)](https://emap.bu.edu/) [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/gayverjr/pyemap/blob/master/LICENSE) [![DOI:10.1021/acs.jpcb.9b04816](https://zenodo.org/badge/DOI/10.1021/acs.jpcb.9b04816.svg)](https://doi.org/10.1021/acs.jpcb.9b04816)

PyeMap is a Python package aimed at automatic identification of electron and hole transfer pathways in proteins. It is intended to serve as the backend for the web application [eMap](https://emap.bu.edu), and as a standalone package.

- **Documentation:** https://readthedocs.org/projects/pyemap/
- **Website:** https://emap.bu.edu
- **News:** https://twitter.com/eMap_protein
# Installation
For more detailed instructions, specifically on getting nicer looking graphs from graphviz please see the [documentation](https://readthedocs.org/projects/pyemap/).
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
Pip installation will only install python dependencies, which is sufficient to run PyeMap analysis, but will be missing some features such as surface exposure and visualization.
```
pip install --extra-index-url https://testpypi.python.org/pypi pyemap
```
# Bugs and feature requests
Please report any bugs and make feature requests [here](https://github.com/gayverjr/pyemap/issues). For issues exlcusive to the web version [eMap](https:emap.bu.edu), please send an email to <emap.bu@gmail.com>. We also greatly encourage users to contribute by making pull requests on [GitHub](https://github.com/gayverjr/pyemap)!

# License
Released under the 3-Clause BSD license (see LICENSE).

Copyright (C) 2019

James Gayvert <jrg444@gmail.com>

Ruslan Tazhigulov <ruslan.tazhigulov@gmail.com>

Ksenia Bravaya <bravaya@bu.edu>
