"""
Unit and regression test for the pyemap package.
"""

# Import package, test suite, and other packages as needed
import pyemap
import pytest
import sys

def test_pyemap_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "pyemap" in sys.modules
