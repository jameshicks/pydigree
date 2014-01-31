pydigree
========

A python package for operations on pedigree and genotype data, including simulation of genotype-phenotype associations under quantitative genetic models

Submodules
-----
In addition to basic pedigree data manipulation, pydigree also includes submodules for more complicated tasks:
* __simulation__: Provides classes for simulating genetic data
* __mixedmodel__: Provides classes for using mixed models with family data
* __io__: Provides functions for importing/exporting data from common data formats, including:
 * `plink`: Functions for working with plink format PED/MAP data
 * `genomesimla`: Includes a function for reading genomeSIMLA format chromosome templates

Classes
-----
* `Invidiual`: Models an individual with pedigree and phenotype data
* `Population`: Models groups of Individuals with a common genetic background
 * `Pedigree`: A special case of Population for related individuals. Implements kinship/inbreeding functions
* `PedigreeCollection`: A container class handling multiple pedigrees
* `Chromosome`: Models a chromosome with information on allele frequency and marker position
* `MixedModel`: A class for fitting mixed-effect models with related individuals
* `Architecture`: A class describing the genetic architecture for a trait to be used in simulation
* `Simulation`: A base class from which other simulation objects inherit
 * `NaiveGeneDroppingSimulation`: Simulates genetic data for pedigrees by random gene dropping
 * `ConstrainedMendelianSimulation`: Simulates genetic data for pedigrees from a prespecified inheritance structure


Scripts
-----
Pydigree includes a few useful scripts for dealing with pedigree data including:
* `simulate_pedigree_data.py`: Simulates data from template pedigrees
* `bitsize.py`: Calculates bit sizes for each pedigree
* `kinship.py`: Caluclates inbreeding coefficients and pairwise kinship coefficients for pedigrees
* `genedrop.py`: Performs gene dropping simulations to approximate the actual probability of an IBD configuration
* `reml.py`: Calculates variance components for normally distributed continuous traits.

### References 
#### Pedigrees
* Charles II of Spain: http://en.wikipedia.org/wiki/Charles_II_of_Spain
