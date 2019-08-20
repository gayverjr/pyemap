parse
==============================
PDB files can be either uploaded, or fetched from the RCSB_ database. At this stage, 
non-protein electron transfer active moieties are automatically identified by pyemap, and
saved to an :ref:`emap <emap>` object which is returned to the user.

.. _RCSB: http://www.rcsb.org/

.. autosummary::
   :toctree: autosummary

   pyemap.parse
   pyemap.fetch_and_parse