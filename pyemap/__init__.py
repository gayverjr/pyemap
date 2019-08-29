"""
pyemap
Implementation of eMap analysis in the form of a python package.
"""

# Add imports here
from .parser import *
from .process_data import process
from .pathway_analysis import find_paths
from .emap import emap
from .shortest_paths import Branch,ShortestPath
# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

