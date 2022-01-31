Branch
==============================
Branches are stored in the `branches` dictionary of the :class:`~pyemap.emap` object. The keys are
the branch numbers.

    >>> print(my_emap.branches[1])
    >>> Branch: W356(A)
    >>> 1a: ['FAD510(A)-2', 'W356(A)'] 9.48
 	>>> 1b: ['FAD510(A)-2', 'W356(A)', 'W213(A)'] 15.31
	>>> 1c: ['FAD510(A)-2', 'W356(A)', 'W213(A)', 'W61(A)'] 23.95
	>>> 1d: ['FAD510(A)-2', 'W356(A)', 'Y432(A)', 'W436(A)'] 26.44
	>>> 1e: ['FAD510(A)-2', 'W356(A)', 'W213(A)', 'W62(A)', 'W217(A)'] 34.72

.. autoclass:: pyemap.Branch
   :members:

   .. automethod:: __str__