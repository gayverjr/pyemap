# PyeMap

[![Build Status](https://travis-ci.org/gayverjr/pyemap.svg?branch=master)](https://travis-ci.org/gayverjr/pyemap) [![codecov](https://codecov.io/gh/gayverjr/pyemap/branch/master/graph/badge.svg)](https://codecov.io/gh/gayverjr/pyemap/branch/master)

PyeMap is a Python package aimed at automatic identification of electron and hole transfer pathways in proteins. It is intended to serve as the backend for the web application [eMap](https://emap.bu.edu), and as a standalone package.

- **Documentation:** https://readthedocs.org/projects/pyemap/ 
- **Website:** https://emap.bu.edu
# Installation
For more detailed instructions, specifically on getting nicer looking graphs please see the [documentation](https://readthedocs.org/projects/pyemap/).
#### Conda(recommended):
The conda recipe will install all dependencies necessary for full functionality.
```
# create new virtual environment
$ conda create -n pyemap_env
$ conda activate pyemap_env
#include channels for dependencies, only needs to be done once
$ conda config --add channels conda-forge --add channels salilab --add channels bioconda 
$ conda update --all
# install pyemap
$ conda install -c gayverjr pyemap
```

### Pip
Pip installation will only install python dependencies, which is sufficient to run PyeMap analysis, but will be missing some features such as surface exposure and visualization. 
```
pip install --extra-index-url https://testpypi.python.org/pypi pyemap
```
For full functionality, you can download and install [MSMS](http://mgltools.scripps.edu/packages/MSMS), [DSSP](https://github.com/cmbi/xssp/releases), and [Graphviz](https://graphviz.gitlab.io/) separately. 

# Bugs and feature requests
Please report any bugs and make feature requests [here](https://github.com/gayverjr/pyemap/issues). For issues exlcusive to the web version [eMap](https:emap.bu.edu), please send an email to <emap.bu@gmail.com>. We also greatly encourage users to contribute by making pull requests on [Github](https://github.com/gayverjr/pyemap)!
 
# License
Released under the 3-Clause BSD license (see LICENSE).

Copyright (C) 2019 PyeMap developers

James Gayvert <jrg444@gmail.com>
Ruslan Tazhigulov <ruslan.tazhigulov@gmail.com>
Ksenia Bravaya <bravaya@bu.edu>




