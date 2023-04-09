<div align="center">
  <img src="https://github.com/gayverjr/pyemap/blob/main/docs/logo/pyemap_logo.png">
</div>

[![GitHub release](https://img.shields.io/github/release/gayverjr/pyemap.svg)](https://github.com/gayverjr/pyemap) 
[![PyPI](https://badge.fury.io/py/pyemap.svg)](https://pypi.org/project/pyemap/)
[![Build Status](https://github.com/gayverjr/pyemap/workflows/ubuntu/badge.svg)](https://github.com/gayverjr/pyemap/actions) [![codecov](https://codecov.io/gh/gayverjr/pyemap/branch/main/graph/badge.svg)](https://codecov.io/gh/gayverjr/pyemap/branch/main) [![Documentation Status](https://readthedocs.org/projects/pyemap/badge/?version=latest)](https://pyemap.readthedocs.io/en/latest/?badge=latest) [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/gayverjr/pyemap/blob/main/LICENSE)

[![Website emap.bu.edu](https://img.shields.io/website-up-down-green-red/https/emap.bu.edu.svg)](https://emap.bu.edu/) [![DOI:10.1021/acs.jpcb.9b04816](https://zenodo.org/badge/DOI/10.1021/acs.jpcb.9b04816.svg)](https://doi.org/10.1021/acs.jpcb.9b04816)[![CodeQl](https://github.com/gayverjr/pyemap/actions/workflows/codeql.yml/badge.svg)](https://github.com/gayverjr/pyemap/actions)[![PyPI version](https://badge.fury.io/py/pyemap.svg)](https://badge.fury.io/py/pyemap)




PyeMap is a Python package aimed at automatic identification of electron and hole transfer pathways in proteins.
It serves as the backend for the web application [eMap](https://emap.bu.edu/), and can also be used as a fully functional Python package.

- **Documentation:** https://readthedocs.org/projects/pyemap/
- **Website:** https://emap.bu.edu

# Installation

### Pip
Pip installation will only install python dependencies, This is sufficient to run PyeMap analysis and view graph images, but some features will be missing.
```
pip install pyemap
```
For full functionality, install [RDKit](https://www.rdkit.org/docs/Install.html), [Graphviz](https://graphviz.gitlab.io/),  [PyGraphviz](https://pygraphviz.github.io/), [MUSCLE](https://www.drive5.com/muscle/), [MSMS](http://mgltools.scripps.edu/packages/MSMS) (not available on Mac OS Catalina and later), and [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html), all of which can be downloaded free of charge from their owners, and are available on most platforms. 

# Getting started
Please see our tutorials for [single](https://pyemap.readthedocs.io/en/latest/tutorial/single_protein.html) and [multiple](https://pyemap.readthedocs.io/en/latest/tutorial/mining.html) protein analysis to help get started using PyeMap. We also provide sample Jupyter notebooks in the examples directory.

# Bugs and feature requests
Please report any bugs and make feature requests [here](https://github.com/gayverjr/pyemap/issues). For issues exclusive to the web version [eMap](https:emap.bu.edu), please send an email to <emap.bu@gmail.com>. We also greatly encourage users to contribute by making pull requests on [GitHub](https://github.com/gayverjr/pyemap)!

# License
Released under the 3-Clause BSD license (see LICENSE).

Copyright (C) 2018-2022

James Gayvert <jrg444@gmail.com>

Ruslan Tazhigulov <ruslan.tazhigulov@gmail.com>

Ksenia Bravaya <bravaya@bu.edu>
