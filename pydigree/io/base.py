from itertools import izip, chain, imap

import numpy as np

from pydigree.common import *
from pydigree.io.smartopen import smartopen as open
from pydigree.population import Population
from pydigree.individual import Individual
from pydigree.pedigree import Pedigree
from pydigree.pedigreecollection import PedigreeCollection
from pydigree.genotypes import Alleles, SparseAlleles


def read_ped(filename, population=None, delimiter=None, affected_labels=None,
             population_handler=None, data_handler=None, connect_inds=True,
             onlyinds=None):
    """
    Reads a plink format pedigree file, ie:
        familyid indid father mother sex whatever whatever whatever
    into a pydigree pedigree object, with optional population to
    assign to pedigree members. If you don't provide a population
    you can't simulate genotypes!

    Arguments
    -----
    filename: The file to be read
    population: The population to assign individuals to
    delimiter: a string defining the field separator, default: any whitespace
    affected_labels: The labels that determine affection status.
    population_handler: a function to set up the population 
    data_handler: a function to turn the data into useful individual information
    connect_inds: build references between individuals. Requires all
        individuals be present in the file
    onlyinds: a list of individuals to be processed, allows skipping parts
        of a file

    Returns: An object of class PedigreeCollection
    """
    sex_codes = {'1': 0, '2': 1, 'M': 0, 'F': 1, '0': None, '-9': None}
    if not affected_labels:
        affected_labels = {'1': 0, '2': 1,
                           'A': 1, 'U': 0,
                           'X': None,
                           '-9': None}

    # Tries to get a phenotype and returns unknown on failure
    def getph(ph):
        try:
            return affected_labels[ph]
        except KeyError:
            return None

    population = Population()

    p = Pedigree()
    if callable(population_handler):
        population_handler(p)

    pc = PedigreeCollection()

    with open(filename) as f:
        # Parse the lines in the file
        for line in f:
            split = line.strip().split(delimiter)
            if len(split) > 5:
                fam, id, fa, mo, sex, aff = split[0:6]
            elif len(split) == 5:
                fam, id, fa, mo, sex = split[0:5]
                aff = None
            # Give a special id for now, to prevent overwriting duplicated
            # ids between families
            id = (fam, id)

            if onlyinds and (id not in onlyinds):
                continue

            p[id] = Individual(population, id, fa, mo, sex)
            p[id].phenotypes['affected'] = getph(aff)
            p[id].pedigree = p
            p[id].sex = sex_codes[p[id].sex]

            if callable(data_handler) and len(split) > 6:
                data = split[6:]
                data_handler(p[id],  data)

    # Fix the individual-level data
    if connect_inds:
        for ind in p:
            fam, id = ind.label
            # Actually make the references instead of just pointing at strings
            ind.father = p[(fam, ind.father)] if ind.father != '0' else None
            ind.mother = p[(fam, ind.mother)] if ind.mother != '0' else None

            ind.register_with_parents()

    # Place individuals into pedigrees
    available = {x for x in p if type(x.label) is tuple}
    for pedigree_label in set(ind.label[0] for ind in p):
        ped = Pedigree(label=pedigree_label)
        if callable(population_handler):
            population_handler(ped)
        
        thisped = {x for x in available if x.label[0] == pedigree_label}
        available = available - thisped
        for ind in thisped:
            ind.label = ind.label[1]
            ped[ind.label] = ind
            ind.population = ped
            ind.pedigree = ped
        pc[pedigree_label] = ped

    return pc


def read_phenotypes(pedigrees, csvfile, delimiter=',', missingcode='X'):
    """
    Reads a csv with header
    famid,ind,phen,phen,phen,phen etc etc

    Arguments
    ------
    Pedigrees:   An object of class PedigreeCollection
    csvfile:     the filename of the file containing phenotypes.
    delimiter:   the field delimiter for the file
    missingcode: the code for missing values

    Returns: Nothing
    """
    with open(csvfile) as f:
        header = f.readline().strip().split(delimiter)
        for line in f:
            # Match columns to their column name
            d = dict(zip(header, line.strip().split(delimiter)))
            for k, v in d.items():
                # Convert all phenotypes into floats
                try:
                    v = float(v)
                except ValueError:
                    if not v:
                        v = None
                if k in set(['famid', 'id']):
                    continue
                fam, ind = d['famid'], d['id']
                pedigrees[fam][ind].phenotypes[k] = v


def write_pedigree(pedigrees, filename, missingcode='X', delim=' '):
    ''' Writes pedigree to a LINKAGE formatted pedigree file '''
    sorting_key = lambda x: (x.population.label, x.depth, x.label)
    with open(filename, 'w') as f:
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
    "Writes phenotypes to a CSV (or other delimited) file"
    inds = pedigrees.individuals

    if callable(predicate):
        inds = [x for x in inds if predicate(x)]

    available_phenotypes = reduce(set.union,
                                  [set(x.phenotypes.keys()) for x in inds])
    available_phenotypes = sorted(available_phenotypes)
    header = ['famid', 'id'] + available_phenotypes

    with open(filename, 'w') as ofile:
        ofile.write(delim.join([str(x) for x in header]) + '\n')
        for ind in inds:
            row = [ind.population.label, ind.label]
            row += [ind.phenotypes.get(phenotype, missingcode)
                    for phenotype in available_phenotypes]
            row = delim.join([str(x) for x in row])
            ofile.write(row + '\n')


def genotypes_from_sequential_alleles(chromosomes, data, missing_code='0',
                                      sparse=False):
    '''
    Takes a series of alleles and turns them into genotypes.

    For example: 
    The series '1 2 1 2 1 2' becomes 
    chrom1 = [1,1,1]
    chrom2 = [2,2,2]

    These are returned in the a list in the form [(chroma, chromb), (chroma, chromb)...]

    Arguments
    ------
    chromosomes: A list of ChromosomeTemplate objects corresponding to the 
    genotypes
    data: The alleles to be turned into genotypes
    sparse: Return SparseAlleless instead of non-sparse

    Returns: A list of 2-tuples of Alleles objects
    '''
    Chromobj = SparseAlleles if sparse else Alleles

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
    for i, chrom in enumerate(chromosomes):
        size = chrom.nmark()
        stop = start + size
        chroma = Chromobj(strand_a[start:stop], template=chrom)
        chromb = Chromobj(strand_b[start:stop], template=chrom)

        genotypes.append((chroma, chromb))
        start += size

    return genotypes
