pydigree
========

A python package for operations on pedigree and genotype data, including simulation of genotype-phenotype associations under quantitative genetic models

Submodules
-----
In addition to basic pedigree data manipulation, pydigree also includes submodules for more complicated tasks:
* __simulation__: Provides classes for simulating genetic data
* __mixedmodel__: Provides classes for using mixed models with family data

Classes
-----
* Invidiual: Models an individual with pedigree and phenotype data
* Population: Models groups of Individuals with a common genetic background
* Pedigree: A special case of Population for related individuals. Implements kinship/inbreeding functions
* PedigreeCollection: A container class handling multiple pedigrees
* Chromosome: Models a chromosome with information on allele frequency and marker position
* MixedModel: A class for fitting mixed-effect models with related individuals
* Architecture: A class describing the genetic architecture for a trait to be used in simulation
* Simulation: A base class from which other simulation objects inherit
** NaiveGeneDroppingSimulation: Simulates genetic data for pedigrees by random gene dropping
** ConstrainedMendelianSimulation: Simulates genetic data for pedigrees from a prespecified inheritance structure


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
