pydigree
========

A python3 package for operations on pedigree and genotype data, including simulation of genotype-phenotype associations under quantitative genetic models

Requires: Python 3.4+, numpy, scipy, pandas, cython

Submodules
-----
In addition to basic pedigree data manipulation, pydigree also includes submodules for more complicated tasks:
* __simulation__: Provides classes for simulating genetic data
* __stats__: Classes and functions for statistical genetics
 * __mixedmodel__: Provides classes for using mixed models with family data
* __io__: Provides functions for importing/exporting data from common data formats, including:
 * `plink`: Functions for working with plink format PED/MAP data
 * `vcf`: Functions for working with the VCF genotype format 
 * `genomesimla`: Includes a function for reading genomeSIMLA format chromosome templates
* __sgs__: Functions for shared genomic segment data

Classes
-----
* `Invidiual`: Models an individual with pedigree and phenotype data
* `Population`: Models groups of Individuals with a common genetic background
 * `Pedigree`: A special case of Population for related individuals. Implements kinship/inbreeding functions
* `PedigreeCollection`: A container class handling multiple pedigrees
* `ChromosomeTemplate`: Models a chromosome with information on allele frequency and marker position
* `ChromosomeSet`: The set of `ChromosomeTemplate`s for a population
* `Alleles`: Stores a haploid set of alleles
 * `SparseAlleles`: Stores a haploid set of alleles as differences from a reference
 * `LabelledAlleles`: An efficient container for storing references to a founder chromosome
* `MixedModel`: A class for fitting mixed-effect models with related individuals
 * `MLEResult`: A class containing the maximum likelihood estimates of parameters and 
 values pertaining to the likelihood function at the MLE
* `Architecture`: A class describing the genetic architecture for a trait to be used in simulation
* `GeneDroppingSimulation`: A base class from which other gene-drop simulation objects inherit
 * `NaiveGeneDroppingSimulation`: Simulates genetic data for pedigrees by random gene dropping
 * `ConstrainedMendelianSimulation`: Simulates genetic data for pedigrees from a prespecified inheritance structure
* `SGSAnalysis`: A class containing the result of a shared genomic segment (SGS) analysis
 * `SGS`: A class containing the segments shared between a pair of individuals
 * `Segment`: A class describing the location of a shared segment between a pair of individuals

Exceptions
-----
* `IterationError`: Raised when at iterative algorithm exceeds the maximum allowed number of iterations
* `NotMeaningfulError`: Raised when a comparison does not make sense (e.g. is one genotype greater than the other)
* `SimulationError`: Raised when an error occurs in a simulation
* `FileFormatError`: Raised when an input file can't be parsed successfully 

Scripts
-----
Pydigree includes a few useful scripts for dealing with pedigree data including:
* `simulate_pedigree_data.py`: Simulates data from template pedigrees
* `bitsize.py`: Calculates bit sizes for each pedigree
* `kinship.py`: Caluclates inbreeding coefficients and pairwise kinship coefficients for pedigrees
* `genedrop.py`: Performs gene dropping simulations to approximate the actual probability of an IBD configuration
* `polygenic.py`: Calculates variance components for normally distributed continuous traits.

### Citation

J.E. Hicks (2017) Pydigree: a python module for manipulation and simulation and of genetic datasets. biorxiv preprint doi:[10.1101/213413](https://doi.org/10.1101/213413)

### References 
#### Pedigrees
* Charles II of Spain: http://en.wikipedia.org/wiki/Charles_II_of_Spain
