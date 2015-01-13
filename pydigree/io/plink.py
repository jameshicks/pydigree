from itertools import izip, chain, imap

import numpy as np

from pydigree.common import grouper, cumsum
from pydigree.chromosome import Chromosome
from pydigree.io.base import read_ped
from pydigree.io.smartopen import smartopen as open

def translate_allele(g):
    return {'A': 1, 'C': 2,
            'G': 3, 'T': 4, '0': 0}[g]


def create_pop_handler_func(mapfile):
    def pop_handler(pop):
        pop.chromosomes = read_map(mapfile)
    return pop_handler


def plink_data_handler(ind, data):
    ind._init_genotypes(blankchroms=False)
    strand_a = data[0::2]
    strand_b = data[1::2]

    chromosomes = ind.chromosomes 
    sizes = [x.nmark() for x in chromosomes]
    starts = cumsum(sizes)
    

    start = 0
    for i, size in enumerate(sizes):
        stop = start + size
        chroma = strand_a[start:stop]
        chromb = strand_b[start:stop]
        ind.genotypes[i] = chroma, chromb
        start += size

def read_map(mapfile):
    """ Reads a PLINK map file into a list of chromosome objects """
    last_chr, last_pos = None, 0
    chroms = []
    chromosome = None
    with open(mapfile) as f:
        for i,line in enumerate(f):
            line = line.strip().split()
            chr, label, cm, pos = line
            cm, pos = float(cm), int(pos)
            if pos < 0:
                continue
            if chr != last_chr:
                # If this happens, we've moved on to a new chromosome,
                # or we've just started. If we haven't just started, We'll 
                # close up the old one
                if i > 0:
                    chroms.append(chromosome)
                # Make the next chromosome
                chromosome = Chromosome(label=chr)
            elif pos < last_pos:
                raise ValueError('Map file not sorted')
            chromosome.add_genotype(None, cm, label=label, bp=pos)
            last_chr, last_pos = chr, pos
    chroms.append(chromosome)
    return chroms


def read_plink(pedfile, mapfile, **kwargs):
    pop_handler = create_pop_handler_func(mapfile)
    return read_ped(pedfile, population_handler=pop_handler, data_handler=plink_data_handler, **kwargs)


def write_plink(pedigrees, filename_prefix, predicate=None, mapfile=False,
                compression=None):

    pedfile = filename_prefix + '.ped'
    if compression in {'gzip', 'gz'}:
        pedfile += '.gz'
    elif compression in {'bzip2', 'bz2'}:
        pedfile += '.bz2'
    write_ped(pedigrees, pedfile, predicate=predicate)
    if mapfile:
        write_map(pedigrees, filename_prefix + '.map')


def write_ped(pedigrees, pedfile,  delim=' ', predicate=None):
    """
    write_ped writes data in a plink-format PED file, and optionally a
    plink-format map file.

    Arguments
    ------

    pedigrees: An object of class PedigreeCollection containing what you
        want to output
    pedfile: a string giving the name out the file to output to.
    mapfile: the name of a mapfile to output, if you want to output one.
        an object that evaluates as False or None will skip the mapfile
    genotypes: Should genotypes be output True/False
    delim: Field seperator
    predicate: Which inputs to include in the output file. If not specified
        all are output. If the string is 'affected', only affected
        individuals are output. If the string is 'phenotyped', all individuals
        with phenotype information are output. Any other value of predicate
        must be a function to perform on the individual that evaluates to
        True/False for whether the individual should be output.

    Returns: Nothing
    """
    if not predicate:
        predicate = lambda x: True
    elif predicate == 'affected':
        predicate = lambda x: x.phenotypes['affected'] == 1
    elif predicate == 'phenotyped':
        predicate = lambda x: x.phenotypes['affected'] in set([0, 1])
    elif not callable(predicate):
        raise ValueError('Not a valid predicate!')

    afflab = {1: '2', 0: '1', None: '-9'}

    with open(pedfile, 'w') as f:
        for pedigree in pedigrees:
            for ind in pedigree:
                if not predicate(ind):
                    continue
                # Get the phenotype code
                aff = afflab[ind.phenotypes['affected']]
                # Prepare the 6-column identifier
                outline = [pedigree.label, ind.id,
                           ind.father.id if ind.father is not None else '0',
                           ind.mother.id if ind.mother is not None else '0',
                           1 if ind.sex == 0 else 2,
                           aff]
                # Get the genotypes in the format we need them
                g = []
                for chroma, chromb in ind.genotypes:
                    g.extend(chain.from_iterable(izip(chroma, chromb)))
                outline.extend(g)
                
                # Make strings
                outline = imap(str, outline)
                # Write it out
                f.write(delim.join(outline))
                f.write('\n')


def write_map(pedigrees, mapfile):
    with open(mapfile, 'w') as f:
        for ci, chromosome in enumerate(pedigrees.chromosomes):
            for mi, marker in enumerate(chromosome._iinfo()):
                label, cm, mb, frequency = marker
                if not mb:
                    mb = int(cm * 10e6)
                if not label:
                    label = 'SNP%s-%s' % (chromosome, mi)
                f.write('\t'.join(str(x) for x
                                  in [chromosome.label, label, cm, mb]) + '\n')
