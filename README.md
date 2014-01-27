pydigree
========

A python package for operations on pedigree and genotype data, including simulation of genotype-phenotype associations under quantitative genetic models

Submodules
-----
In addition to basic pedigree data manipulation, pydigree also includes submodules for more complicated tasks:
* __simulation__: Provides classes for simulating genetic data
* __mixedmodel__: Provides classes for using mixed models with family data

Scripts
-----
Pydigree includes a few useful scripts for dealing with pedigree data including:
* __simulate\_pedigree\_data.py__: Simulates data from template pedigrees
* __bitsize.py__: Calculates bit sizes for each pedigree
* __kinship.py__: Caluclates inbreeding coefficients and pairwise kinship coefficients for pedigrees
* __genedrop.py__: Performs gene dropping simulations to approximate the actual probability of an IBD configuration
* __reml.py__: Calculates variance components for normally distributed continuous traits.

### References 
#### Pedigrees
* Charles II of Spain: http://en.wikipedia.org/wiki/Charles_II_of_Spain
