Usage
========================
Graph mining with PyeMap occurs in 4 steps:

#. Generate protein graphs

#. Classify nodes and edges

#. Find subgraph patterns

#. Find and cluster protein subgraphs


**Step 1: Generate Protein Graphs**

The first step is to fetch and/or parse a list of PDBs using :func:`~pyemap.fetch_and_parse` or :func:`~pyemap.parse`, and add them 
to the :class:`~pyemap.graph_mining.PDBGroup` object. Once all of the eMap objects have been collected, the second step is to 
call :func:`~pyemap.graph_mining.PDBGroup.process_emaps`, which uses
the same infrastructure as :ref:`Single Protein Analysis <spa>` to generate the protein graphs.

.. code-block:: python

    pg = pyemap.graph_mining.PDBGroup("My Group")
    #pdb_ids = ['1u3d','1iqr'...]
    for pdb in pdb_ids: 
        pg.add_emap(pyemap.fetch_and_parse(pdb)) 

**Step 2: Classify nodes and edges**

The next step is to classify the nodes and edges using :func:`~pyemap.graph_mining.PDBGroup.generate_graph_database`. In 
many cases, this function can be called with no arguments, but in some cases it can be useful to allow for node substitutions, or to specify 
edge thresholds. See the :ref:`classification <classify>` section for more details.

.. code-block:: python

   pg.generate_graph_database()

**Step 3: Find Subgraph Patterns**

The next step is to find subgraph patterns which are shared among the protein graphs. One can either mine for 
all patterns up to a given support threshold:

.. code-block:: python

    pg.run_gspan(10)

Or search for a particular pattern or set of patterns:

.. code-block:: python

    pg.find_subgraph('WWW#')

See the :ref:`mining <mining_algo>` section for more details.

**Step 4: Find and cluster protein subgraphs**

The results of the mining calculation are stored in the :attr:`subgraph_patterns` dictionary as 
:class:`~pyemap.graph_mining.SubgraphPattern` objects. To find protein subgraphs, 
the :func:`~pyemap.graph_mining.SubgraphPattern.find_protein_subgraphs` function 
must be called for the subgraph pattern of interest. The identified protein subgraphs are stored in 
the :attr:`protein_subgraphs` attribute of the :class:`~pyemap.graph_mining.SubgraphPattern` object, 
and the clustering is described by the :attr:`groups` attribute. One can switch between different types of 
clustering using the :func:`~pyemap.graph_mining.SubgraphPattern.set_clustering` function.

.. code-block:: python

    sg = pg.subgraph_patterns['1_WWW#_18']
    sg.find_protein_subgraphs()
    sg.set_clustering("sequence")
    # print results, including clustering
    print(sg.full_report())

