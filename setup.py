"""PyeMap: A python package for automatic identification of electron and hole transfer pathways in proteins.

PyeMap (pronounced Pie-Map) is a python package for automatic identification and visualization of electron and hole transfer pathways in proteins. 
For single PDB structures, PyeMap identifies and visualizes possible efficient pathways from a specified source residue to the surface 
or to a specified target residue. For a group of proteins, PyeMap identifies common pathways/motifs which are shared among the 
PDB structures by mining for common patterns in the set of protein graphs. These shared pathways are then clustered into 
groups based on their similarity.

All PyeMap wheels distributed on PyPI are 3-BSD licensed.

GitHub: https://github.com/gayverjr/pyemap

Documentation: https://pyemap.readthedocs.io/en/latest/
"""
import sys
from setuptools import setup, find_packages

DOCLINES = (__doc__ or '').split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []


setup(
    # Self-descriptive entries which should always be present
    name='pyemap',
    author='James Gayvert',
    author_email='jrg444@gmail.com',
    description= DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    version="2.0.0",
    license='BSD-3-Clause',

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # url='http://www.my_package.com',  # Website
    install_requires=['numpy',
    				  'networkx',
    				  'biopython',
    			      'scipy',
    				  'pillow',
                      'reportlab',
                      'pygraphviz',
                      'gspan_mining',
                      'requests',
                      'svglib'
                      ]
    				   # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    #python_requires=">=3.7"         # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
