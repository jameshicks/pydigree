from itertools import chain

from pydigree.common import grouper, cumsum
from pydigree.common import interleave
from pydigree.genotypes import ChromosomeTemplate, ChromosomeSet
from pydigree.io.base import read_ped, genotypes_from_sequential_alleles
from pydigree.exceptions import FileFormatError
from pydigree.io.smartopen import smartopen as open
import collections


def create_pop_handler_func(mapfile):
    chromosomes = read_map(mapfile)

    def pop_handler(pop):
        for chrom in chromosomes:
            pop.add_chromosome(chrom)

    return pop_handler


def plink_data_handler(ind, data):
    ''' A function to handle the data payload from a plink line '''
    ind.genotypes = genotypes_from_sequential_alleles(
        ind.chromosomes, data, missing_code='0')


def read_map(mapfile):
    """
    Reads a PLINK map file into a list of ChromosomeTemplate objects

    
    :param mapfile: Path of the file to be read
    :type mapfile: string

    :rtype: a list of ChromosomeTemplate objects
    """
    last_chr, last_pos = None, 0
    chroms = ChromosomeSet()
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
                    chromosome.finalize()
                    chroms.add_chromosome(chromosome)
                # Make the next chromosome
                chromosome = ChromosomeTemplate(label=chr)
            elif pos < last_pos:
                raise FileFormatError('Map file not sorted')
            chromosome.add_genotype(None, cm, label=label, bp=pos)
            last_chr, last_pos = chr, pos
    chromosome.finalize()
    chroms.add_chromosome(chromosome)
    return chroms


def read_plink(pedfile=None, mapfile=None, prefix=None, **kwargs):
    '''
    Read a plink file by specifying pedfile and mapfile directly,
    or by using a prefix. Pass additional arguments to 
    pydigree.io.base.read_ped with kwargs

    
    :param pedfile: a plink PED file to be read
    :param mapfile: a plink MAP file to be read
    :param prefix: sets mapfile to 'prefix.map' and pedfile to 'prefix.ped'
    :param kwargs: additional arguments passed to read_ped
    
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
                compression=None, output_chromosomes=None):
    '''
    Write individual genotypes to a file in plink PED data format.
    Optionally outputs the genotype locations to the mapfile.

    :param pedigrees: the pedigrees to be written
    :param filename_prefix: The output ped file ('.ped' will be appended)
    :param predicate: a callable that evaluates True when the individual 
      should be outputted to the file
    :param mapfile: True if the plink MAP file should be written
    :param compression: Compress the data? Options are bz2 and gzip,
      otherwise data will be written uncompressed

    :type pedigrees: IndividualContainer
    :type filename_prefix: string
    :type predicate: callable
    :type mapfile: string
    :type compression: 'gzip' or 'bz2'

    :rtype: void
    '''
    pedfile = filename_prefix + '.ped'
    if compression in {'gzip', 'gz'}:
        pedfile += '.gz'
    elif compression in {'bzip2', 'bz2'}:
        pedfile += '.bz2'
    write_ped(pedigrees, pedfile, predicate=predicate, 
                output_chromosomes=output_chromosomes)
    if mapfile:
        write_map(pedigrees, filename_prefix + '.map', output_chromosomes=output_chromosomes)


def write_ped(pedigrees, pedfile,  delim=' ', predicate=None,
              output_chromosomes=None):
    """
    write_ped writes data in a plink-format PED file, and optionally a
    plink-format map file.


    :param pedigrees: An object of class PedigreeCollection containing what you
        want to output
    :param pedfile: a string giving the name out the file to output to.
    :param mapfile: the name of a mapfile to output, if you want to output one.
        an object that evaluates as False or None will skip the mapfile
    :param genotypes: Should genotypes be output True/False
    :param delim: Field seperator
    :param predicate: Which inputs to include in the output file. If not 
        specified all are output. If the string is 'affected', only affected
        individuals are output. If the string is 'phenotyped', all individuals
        with phenotype information are output. Any other value of predicate
        must be a function to perform on the individual that evaluates to
        True/False for whether the individual should be output.

    Returns: Nothing
    """

    # Check if we're only supposed to be outputting certain chromosomes
    if output_chromosomes is not None:
        checkchroms = True
    else:
        checkchroms = False

    if not predicate:
        predicate = lambda x: True
    elif predicate == 'affected':
        predicate = lambda x: x.phenotypes['affected'] == 1
    elif predicate == 'phenotyped':
        predicate = lambda x: x.phenotypes['affected'] in set([0, 1])
    elif not isinstance(predicate, collections.Callable):
        raise ValueError('Not a valid predicate!')

    afflab = {1: '2', 0: '1', None: '-9'}

    with open(pedfile, 'w') as f:
        for pedigree in pedigrees.pedigrees:
            for ind in pedigree.individuals:
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
                # Make strings
                outline = list(map(str, outline))

                # Get the genotypes in the format we need them
                g = []
                for template, chromatids in zip(ind.chromosomes, ind.genotypes):
                    if checkchroms and template.outputlabel not in output_chromosomes:
                        continue
                    chroma, chromb = chromatids
                    ga = chroma.astype(str).tolist()
                    gb = chromb.astype(str).tolist()
                    gn = interleave(ga, gb)
                    g.extend(gn)
                outline.extend(g)

                # Write it out
                outline = delim.join(outline)
                f.write(outline)
                f.write('\n')


def write_map(pedigrees, mapfile, output_chromosomes=None):
    '''
    Writes the genotype location data to a PLINK MAP file

    :param pedigrees: the population containing the data to be written
    :param mapfile: the name of the file to be output to

    Returns: Nothing
    '''
    # Check if we're only supposed to be outputting certain chromosomes
    if output_chromosomes is not None:
        checkchroms = True
    else:
        checkchroms = False

    with open(mapfile, 'w') as f:
        for ci, chromosome in enumerate(pedigrees.chromosomes):
            if checkchroms and chromosome.outputlabel not in output_chromosomes:
                continue
            for mi, marker in enumerate(chromosome._iinfo()):
                label, cm, mb, frequency = marker
                if not mb:
                    mb = int(cm * 10e6)
                if not label:
                    label = 'SNP%s-%s' % (chromosome.outputlabel, mi)
                f.write('\t'.join(str(x) for x
                                  in [chromosome.outputlabel, label, cm, mb]) + '\n')
