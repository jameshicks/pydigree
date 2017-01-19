File I/O
=========

Pedigree Files
--------------
Pydigree supports pedigree files specified in the following delimited format:

:: 
    
    PEDID ID FATHER MOTHER SEX

Reading these is supported through `pydigree.io.read_ped`.
Pedigree files are not required to be sorted, and IDs are not required to be numeric.
These files may contain extra data following sex, which can be handled with a user-specified function.



PLINK format
------------

PLINK format is the most common file format used for GWAS data. 
Pydigree provides input and output for PLINK formatted data in the 
`pydigree.io.plink` module. 
Genotypes must be in order per chromosome, or an exception will be thrown. 
Binary and transposed plink files are not currently supported.


VCF format
----------

Variant Call Format files can contain large amounts of genotype data.
Since the bulk of (human) genome variation is rare, naive input of VCF files
can result in huge amounts of memory usage. 
Pydigree reads VCF files into the SparseAlleles container to avoid this.
This functionality is found in `pydigree.io.vcf`. 
VCF files can contain large amounts of metadata, and parsing them is more intensive than other formats.

Pydigree does not store any metadata when reading VCF files, so writing them is not supported.
