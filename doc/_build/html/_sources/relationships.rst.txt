Pedigrees, Relationships, and Genealogies 
=========================================
In both pedigree datasets and forward-time simulations, pydigree stores genealogy data.
This data is broadly useful, and pydigree provides functionality for navigating and interpreting them.

Measures of relatedness
-----------------------

Kinship
^^^^^^^
Malecot's coefficient of coancestry (commonly called the kinship coefficient, :math:`\theta`) is a common measure of relatedness.
There are two main ways calculating this. 
The first involves enumerating all paths between two individuals in the genealogy. 
The other involves the recursively calculating kinships through the genealogy, since each relationship is a function of the parents' relationship.

When using pedigrees or other data where the number of pairwise relationships is relatively small, the recursive method will be more efficient since intermediate values can be cached. 
Caching also allows for fast computation of relationships for all individuals at once.
This is the method implemented in the Population class. 
Conversely, when the genealogy is large or deep, such as in population simulations, recursive calculation is prohibitively expensive.
The large number of pairs usually means that exhaustive calculation of pairwise relationships is not available (or desired).
The path-based kinship calculations may be more efficient for selected pairs of individuals.
These are implemented in `pydigree.paths`.  



Fraternity
^^^^^^^^^^
Fraternity coefficient (:math:`\Delta`) describes the probablity that both alleles in a pair of individuals are shared IBD. 
This is uniformly less likely than sharing only one pair of alleles IBD, and is smaller than the kinship coefficent. 
For two individuals, x and y, coefficient is defined as 

.. math::
    \Delta_{x,y} = \phi(x_m,y_m) \phi(x_f,y_f) + \phi(x_m,y_f) \phi(x_f,y_m)

where :math:`\phi(q,w)` is the kinship coefficent between individuals q and w, and :math:`a_b` is the father (b=f) or mother (b=m) of individual a.  

As with kinship coefficients, when working with pedigree data, the recursive implementations implemented in the pedigree objects are likely more efficient. 

Inbreeding
^^^^^^^^^^
Inbreeding coefficients (:math:`F`) represent the probability that two alleles within the same individual are IBD. 
This is a direct function of the relationship of the parents, and is calculated based on their kinship coefficent. 

Relationship Covariance Matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Statistical models incorporating relationships between individuals require matrices of covariance measures, which are calculated from relationships. 
Pydigree provides matrices for additive and dominance covariance matrices. 
The additive relationship matrix (:math:`\mathbf{A}`) is a function of kinships and is defined as

.. math::
    \mathbf{A}_{i,j} = 
    \begin{Bmatrix}
    2 \theta_{i,j} \mathrm{   if } i \neq j
    \\ 
    1 + F_i  \mathrm{   if } i = j
    \end{Bmatrix}

The dominance covariance matrix (:math:`\mathbf{D}`) is a function of fraternity:

.. math::
    \mathbf{D}_{i,j} = 
    \begin{Bmatrix}
    \Delta_{i,j} \mathrm{   if } i \neq j
    \\ 
    1 \mathrm{   if } i = j
    \end{Bmatrix}


Paths through Genealogies
-------------------------
Valid inheritance paths through pedigrees have one major difference to general graph traversal: you can only change direction once. 
When one individual is not a descendant of the other, paths move 'up' (backwards in time, to older generations) from one individual, then 'down' (forwards in time, to newer generations). 
Otherwise, they only move in one direction.

When general paths through genealogies are needed (for example, calculating allele transmission probabilities), pydigree provides useful functions in :py:mod:`pydigree.paths`.
`paths.common_ancestors` finds all shared ancestors between two individuals. 
`paths.paths` finds all paths between two individuals. 
`paths.paths_through_ancestors` will limit the returned paths to ones flowing through a specified common ancestor of two individuals.

Two special pathing functions are also provided for special paths. :py:meth:`pydigree.Individual.patriline` and :py:meth:`pydigree.Individual.matriline` recursively search through maternal and paternal lineage to find the founder of the lineage. 
This is useful for models incorporating mtDNA or Y chromosome inheritance. 
 
