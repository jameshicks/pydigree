import numpy as np

from pydigree.io.smartopen import smartopen
from pydigree.population import Population
from pydigree.individual import Individual
from pydigree.pedigree import Pedigree
from pydigree.pedigreecollection import PedigreeCollection
from pydigree.genotypes import Alleles, SparseAlleles

from collections import Callable
from functools import reduce

sex_codes = {'1': 0, '2': 1, 
             'M': 0, 'F': 1, 
             '0': None, '-9': None} 


class PEDRecord(object):
    def __init__(self, line, delimiter=' '):
        """
        Creates pedigree record from line 

        :param line: a line in the pedigree file
        :param delmiter: field separator
        :type line: string
        :type delimiter: string
        """

        split = line.strip().split(delimiter)
        self.fam, self.ind_id, self.fa, self.mo, self.sex = split[0:5]
        self.aff = split[5] if len(split) > 5 else None
        self.data = split[6:] if len(split) > 6 else None


    def create_individual(self, population=None):
        """
        Creates an Individual object from a Pedigree Record.

        The individual will have the id tuple of (fam_id, ind_id)
        
        :param population: Population for the individual to belong to
        :type population: Population

        :rtype: Individual
        """
        
        # Give a special ind_id for now to prevent overwriting duplicated
        # ind_ids between families
        temp_id = (self.fam, self.ind_id)    

        ind = Individual(population, temp_id, 
                         self.fa, self.mo, 
                         sex_codes[self.sex])    

        return ind


def connect_individuals(pop):
    """
    Makes the connections in the genealogy from parents to children
    and vice versa.

    :param population: Set of individuals to connect
    :type population: IndividualContainer

    :rtype void:
    """

    for ind in pop.individuals:   
        fam, _ = ind.label
        
        ind.father = pop[(fam, ind.father)] if ind.father != '0' else None
        ind.mother = pop[(fam, ind.mother)] if ind.mother != '0' else None

        ind.register_with_parents()


def sort_pedigrees(inds, population_handler):
    """
    Takes a set of individuals and sorts them into pedigrees.

    Individuals must have labels that are (pedid, indid) tuples

    :param inds: Individuals to be sorted
    :param population_handler: a function to set up the population
    :type inds: iterable
    :type population_handler: Callable

    :returns: Collection of pedigrees from the individuals
    :rtype: PedigreeCollections
    """
    if not isinstance(population_handler, Callable):
        population_handler = lambda *x: None

    pc = PedigreeCollection()

    inds = list(inds)
    # A dict that maps pedigree labels to the corresponding individuals
    pedigrees = {ind.label[0]: [] for ind in inds}

    for ind in inds:
        pedigrees[ind.label[0]].append(ind)


    for pedigree_label, ped_inds in pedigrees.items():
        ped = Pedigree(label=pedigree_label)

        population_handler(ped)
        
        # Fix the labels 
        for ind in ped_inds:
            ind.label = ind.label[1]
            ped[ind.label] = ind
            ind.population = ped
            ind.pedigree = ped
        
        pc[pedigree_label] = ped

    return pc

def read_ped(filename, population=None, delimiter=None, affected_labels=None,
             population_handler=None, data_handler=None, connect_inds=True,
             onlyinds=None):
    """
    Reads a plink format pedigree file, ie:
    
    ::    
        familyid indid father mother sex whatever whatever whatever
    
    into a pydigree pedigree object, with optional population to
    assign to pedigree members. If you don't provide a population
    you can't simulate genotypes!


    :param filename: The file to be read
    :param population: The population to assign individuals to
    :param delimiter: a string defining the field separator, 
        default: any whitespace
    :param affected_labels: The labels that determine affection status.
    :param population_handler: a function to set up the population 
    :param data_handler: a function to turn the 
        data into useful individual information
    :param connect_inds: build references between individuals. Requires all
        individuals be present in the file
    :param onlyinds: only include data for specified individuals 

    :type filename: string
    :type population: Population
    :type delimiter: string
    :type affected_labels: dict (str -> value)
    :type data_handler: callable
    :type connect_inds: bool
    :type onlyinds: iterable 


    :returns: individuals contained in the pedigree file 
    :rtype: PedigreeCollection
    """

    if not affected_labels:
        affected_labels = {'1': 0, '2': 1, 'A': 1, 'U': 0,
                           'X': None, '-9': None}

    if not isinstance(data_handler, Callable):
        data_handler = lambda *x: None
    
    if not isinstance(population_handler, Callable):
        population_handler = lambda *x: None

    population = Population() if population is None else population
    p = Pedigree()

    population_handler(p)

    # Step 1: Read the data and create the individuals
    with smartopen(filename) as f:
        # Parse the lines in the file
        for line in f:
            rec = PEDRecord(line, delimiter)
            
            if onlyinds and (rec.ind_id not in onlyinds):
                continue

            ind = rec.create_individual(population)
            ind.pedigree = p
            ind.phenotypes['affected'] = affected_labels.get(rec.aff, None)
            p[ind.label] = ind

            if rec.data:
                data_handler(p[ind.label], rec.data)

    # Step 2: Create the between-individual relationships

    # Fix the individual-level data: individuals currently only have parent-ids
    # in their parent fields and not references to actual individuals
    if connect_inds:
        connect_individuals(p)

    # Step 3: Separate the individuals into pedigrees
    pc = sort_pedigrees(p.individuals, population_handler)


    return pc


def read_phenotypes(pedigrees, csvfile, delimiter=',', missingcode='X'):
    """
    Reads a csv with header
    famid, ind, phen, phen, phen, phen etc etc

    Arguments
    :param pedigrees:   data to update
    :param csvfile:     the filename of the file containing phenotypes.
    :param delimiter:   the field delimiter for the file
    :param missingcode: the code for missing values
    :type pedigrees: PedigreeCollection
    :type csvfile: string
    :type missingcode: string

    :rtype: void
    """
    with smartopen(csvfile) as f:
        header = f.readline().strip().split(delimiter)
        for line in f:
            # Match columns to their column name
            d = dict(list(zip(header, line.strip().split(delimiter))))
            for k, v in list(d.items()):
                
                # Convert all phenotypes into floats
                try:
                    v = float(v)
                except ValueError:
                    if not v or v == missingcode:
                        v = None
                
                if k in set(['famid', 'id']):
                    continue
                
                fam, ind = d['famid'], d['id']
                pedigrees[fam][ind].phenotypes[k] = v


def write_pedigree(pedigrees, filename, delim=' '):
    ''' 
    Writes pedigree to a LINKAGE formatted pedigree file 

    :param pedigrees: Data to write
    :param filename: filename to write to
    :param delim: output field separator 

    :rtype: void
    '''
    sorting_key = lambda x: (x.population.label, x.depth, x.label)
    with smartopen(filename, 'w') as f:
        for ind in sorted(pedigrees.individuals, key=sorting_key):
            oline = [ind.population.label,
                     ind.label,
                     '0' if ind.is_founder() else ind.father.label,
                     '0' if ind.is_founder() else ind.mother.label,
                     '1' if ind.sex == 1 else '0',
                     '-9']
            
            oline = delim.join(oline)
            
            f.write(oline + '\n')


def write_phenotypes(pedigrees, filename, predicate=None,
                     missingcode='X', delim=','):
    """
    Writes phenotypes to a CSV (or other field delimited) file
    
    :param pedigrees: Data to write
    :param filename: filename to write to
    :param missingcode: code to use for missing values
    :param delim: output field separator 

    :type missingcode: string
    :type delim: string
    """
    inds = pedigrees.individuals

    if isinstance(predicate, Callable):
        inds = [x for x in inds if predicate(x)]

    available_phenotypes = reduce(set.union,
                                  [set(x.phenotypes.keys()) for x in inds])
    available_phenotypes = sorted(available_phenotypes)
    header = ['famid', 'id'] + available_phenotypes

    with smartopen(filename, 'w') as ofile:
        ofile.write(delim.join([str(x) for x in header]) + '\n')
        for ind in inds:
            row = [ind.population.label, ind.label]
            row += [ind.phenotypes.get(phenotype, missingcode)
                    for phenotype in available_phenotypes]
            row = delim.join([str(x) for x in row])
            ofile.write(row + '\n')


def genotypes_from_sequential_alleles(chromosomes, data, missing_code='0'):
    '''
    Takes a series of alleles and turns them into genotypes.

    For example: 
    The series '1 2 1 2 1 2' becomes 
    chrom1 = [1, 1, 1]
    chrom2 = [2, 2, 2]

    These are returned in the a list in the form:

    ::
        [(chroma, chromb), (chroma, chromb)...]


    :param chromosomes: genotype data
    :param data: The alleles to be turned into genotypes
    :param missing_code: value representing a missing allele

    :type chromosomes: list of ChromosomeTemplate
    :type missing_code: string
    :returns: A list of 2-tuples of Alleles objects
    '''

    genotypes = []

    data = np.array(data)

    if not np.issubdtype(type(missing_code), data.dtype):
        raise ValueError(
            'Invalid type for missing code: {}. Expected: {}'.format(
                type(missing_code), data.dtype))

    if np.issubdtype(data.dtype, str):
        data[data == missing_code] = ''
    else:
        data[data == missing_code] = 0

    strand_a = data[0::2]
    strand_b = data[1::2]

    start = 0
    for chrom in chromosomes:
        size = chrom.nmark()
        stop = start + size
        chroma = Alleles(strand_a[start:stop], template=chrom)
        chromb = Alleles(strand_b[start:stop], template=chrom)

        genotypes.append((chroma, chromb))
        start += size

    return genotypes
