from itertools import izip, chain, imap

from pydigree.common import grouper
from pydigree.io import read_ped


def translate_allele(g):
    return {'A': 1, 'C': 2, 'G': 3, 'T': 4}[g]


def create_pop_handler_func(mapfile):
    def pop_handler(pop):
        pop.chromosomes = read_map(mapfile)
    return pop_handler


def plink_data_handler(ind, data):
    genotypes = grouper(data, 2)
    chromosomes = ind.population.chromosomes
    try:
        for c_idx, chromosome in enumerate(chromosomes):
            for m_idx, marker in enumerate(chromosome):
                geno = genotypes.next()
                locus = c_idx, m_idx
                ind.set_genotype(locus, geno)
    except StopIteration:
        raise ValueError('Ran out of genotypes!')
    except KeyError:
        raise ValueError('Invalid genotype: %s' % )


def read_map(mapfile):
    """ Reads a PLINK map file into a list of chromosome objects """
    last_chr, last_pos = None, 0
    chroms = []
    chromosome = pydigree.Chromosome()
    with open(mapfile) as f:
        for line in f:
            line = line.strip().split()
            chr, label, cm, pos = line
            if int(pos) < 0:
                continue
            if chr != last_chr:
                # If this happens, we've moved on to a new chromosome
                # We'll close up the old one and start working on a new
                # object
                chroms.append(chromosome)
                chromosome = pydigree.Chromosome()
            elif pos < last_pos:
                raise ValueError('Map file not sorted')
            chromosome.add_genotype(None, cm, label=label, bp=pos)
            last_chr, last_pos = chr, pos
    return chroms


def read_plink(pedfile, mapfile):
    pop_handler = create_pop_handler_func(mapfile)
    return read_ped(pedfile, population_handler=pop_handler, data_handler=plink_data_handler


def write_ped(pedigrees, pedfile, mapfile=None, genotypes=True, delim=' ',
              predicate=None):
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
                if genotypes:
                    g = []
                    for chroma, chromb in ind.genotypes:
                        g.extend(chain.from_iterable(izip(chroma, chromb)))
                    outline.extend(g)
                # Make strings
                outline = imap(str, outline)
                # Write it out
                f.write(delim.join(outline))
                f.write('\n')

    if not mapfile:
        return

    with open(mapfile, 'w') as f:
        for ci, chromosome in enumerate(pedigrees.chromosomes):
            for mi, marker in enumerate(chromosome._iinfo):
                label, cm, mb, frequency = marker
                if not mb:
                    mb = int(cm * 10e6)
                if not label:
                    label = 'SNP%s-%s' % (chromosome, mi)
                f.write('\t'.join(str(x) for x
                                  in [chromosome, label, cm, mb]) + '\n')
