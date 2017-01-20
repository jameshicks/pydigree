Basic datastructures
====================

Individuals
-----------
The :class:`Individual` class is the basic unit of pydigree. 
:class:`Individual` s represent sets of phenotypes and genotypes, with associated functions for working with them. 
When genealogy is present, individuals store references to their parents and children, and can use those references to compute more complex relationships (siblings, half-siblings, descendents, and ancestors).
Otherwise, individuals do not do much, but are acted on by other objects and functions.


Collections of Individuals
--------------------------
Groups of individuals are populations.
Implementation of common population functions are found in the the mixin class :class:`IndividualContainer`, and include selection of classes of individuals, batch functions on them, and allele frequency calculation. 

The most basic individual collection is :class:`Population`.
This collection is useful for datasets of unrelated individuals or populations from forward-time simulations. 
This is the only collection that supports general forward time simulation. 

Pedigrees and Pedigree Collections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Sets of individuals with known relationships have their own specialized containers. 
:class:`Pedigree` contains a single kindred and has efficient methods for quantifying the relationships between individuals and performing gene-dropping simulations. 
:class:`PedigreeCollection` holds multiple :class:`Pedigree` objects and preforms batch functions on them.  

Genotypes
---------
Genotype storage in pydigree comes in two parts.
:class:`ChromosomeTemplate` stores the basic information for each variant in a dataset, such as frequency, position, and variant name.
They are collected into :class:`ChromosomeSet` objects. 
:class:`ChromosomeSet` instance references are retained by the individual collections.
Each individual will have a complement of allele containers corresponding to each template in the population. 
These support basic slicing and copying, as well as finding missing values.

The :class:`Alleles` object stores a sequence of haploid variants in an efficent manner. 
When working with large datasets :class:`Alleles` may not be the optimal container.
For example, sequencing data can contain millions of variants.
Since most variation in the human genome is rare, storing an allele for every position would require memory for each wild-type allele.
Pydigree provides the :class:`SparseAlleles` container, which only stores minor alleles in a tree datastructure. 
This greatly reduces memory usage, but incurs a lookup penalty: element access is O(log n) for :py:class:`SparseAlleles`, compared to O(1) for the standard :class:`Alleles`. 
In datasets where genotypes are mostly dense (*i.e.* most variants are common), :class:`SparseAlleles` may show *higher* memory usage due to the bookkeeping for the tree.

