from itertools import izip, chain, imap

from pydigree.common import grouper, cumsum
from pydigree.genotypes import ChromosomeTemplate
from pydigree.io.base import read_ped, genotypes_from_sequential_alleles
from pydigree.exceptions import FileFormatError
from pydigree.io.smartopen import smartopen as open


def create_pop_handler_func(mapfile):
    chromosomes = read_map(mapfile)

    def pop_handler(pop):
        for chrom in chromosomes:
            pop.add_chromosome(chrom)

    return pop_handler


def plink_data_handler(ind, data):
    ''' A function to handle the data payload from a plink line '''
    ind.genotypes = genotypes_from_sequential_alleles(ind.chromosomes, data, missing_code='0')


def read_map(mapfile):
    """
    Reads a PLINK map file into a list of ChromosomeTemplate objects
    
    Arguments:
    mapfile: The file to be read

    Returns: a list of ChromosomeTemplate objects
    """
    last_chr, last_pos = None, 0
    chroms = []
    chromosome = None
    with open(mapfile) as f:
        for i, line in enumerate(f):
            line = line.strip().split()
            chr, label, cm, pos = line
            cm, pos = float(cm), int(pos)
            if pos < 0:
                raise FileFormatError('Invalid position: {}'.format(pos))
            if chr != last_chr:
                # If this happens, we've moved on to a new chromosome,
                # or we've just started. If we haven't just started, We'll
                # close up the old one
                if i > 0:
                    chroms.append(chromosome)
                # Make the next chromosome
                chromosome = ChromosomeTemplate(label=chr)
            elif pos < last_pos:
                raise FileFormatError('Map file not sorted')
            chromosome.add_genotype(None, cm, label=label, bp=pos)
            last_chr, last_pos = chr, pos
    chroms.append(chromosome)
    return chroms


def read_plink(pedfile=None, mapfile=None, prefix=None, **kwargs):
    '''
    Read a plink file by specifying pedfile and mapfile directly,
    or by using a prefix. Pass additional arguments to 
    pydigree.io.base.read_ped with kwargs

    Arguments:
    pedfile: a plink PED file to be read
    mapfile: a plink MAP file to be read
    prefix: sets mapfile to 'prefix.map' and pedfile to 'prefix.ped'
    kwargs: additional arguments passed to read_ped

    Returns: A PedigreeCollection object
    '''
    if prefix:
        pedfile = prefix + '.ped'
        mapfile = prefix + '.map'

    pop_handler = create_pop_handler_func(mapfile)
    return read_ped(pedfile, population_handler=pop_handler,
                    data_handler=plink_data_handler, connect_inds=False,
                    **kwargs)


def write_plink(pedigrees, filename_prefix, predicate=None, mapfile=False,
                compression=None):
    '''
    Write individual genotypes to a file in plink PED data format.
    Optionally outputs the genotype locations to the mapfile.

    Args:
    ------
    pedigrees: the pedigrees to be written
    filename_prefix: The output ped file ('.ped' will be appended)
    predicate: a callable that evaluates True when the individual 
      should be outputted to the file
    mapfile: True if the plink MAP file should be written
    compression: Compress the data? Options are bz2 and gzip2,
      otherwise data will be written uncompressed

    Returns: Nothing
    '''
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
                outline = [pedigree.label, ind.label,
                           ind.father.label if ind.father is not None else '0',
                           ind.mother.label if ind.mother is not None else '0',
                           1 if ind.sex == 0 else 2,
                           aff]
                # Get the genotypes in the format we need them
                g = []
                for chroma, chromb in ind.genotypes:
                    if ind.translated_alleles:
                        chroma = [reverse_translate_allele(a) for a in chroma]
                        chromb = [reverse_translate_allele(a) for a in chromb]
                    g.extend(chain.from_iterable(izip(chroma, chromb)))
                outline.extend(g)

                # Make strings
                outline = imap(str, outline)
                # Write it out
                f.write(delim.join(outline))
                f.write('\n')


def write_map(pedigrees, mapfile):
    '''
    Writes the genotype location data to a PLINK MAP file
    
    Arguments
    ------
    pedigrees: the population containing the data to be written
    mapfile: the name of the file to be output to

    Returns: Nothing
    '''
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
