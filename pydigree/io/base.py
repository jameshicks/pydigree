from itertools import chain

import numpy as np

from pydigree.common import *
from pydigree.io.smartopen import smartopen as open
from pydigree.population import Population
from pydigree.individual import Individual
from pydigree.pedigree import Pedigree
from pydigree.pedigreecollection import PedigreeCollection
from pydigree.genotypes import Alleles, SparseAlleles

from collections import Callable
from functools import reduce


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
    :param delimiter: a string defining the field separator, default: any whitespace
    :param affected_labels: The labels that determine affection status.
    :param population_handler: a function to set up the population 
    :param data_handler: a function to turn the data into useful individual information
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
    if isinstance(population_handler, Callable):
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

            if isinstance(data_handler, Callable) and len(split) > 6:
                data = split[6:]
                data_handler(p[id],  data)

    # Fix the individual-level data
    if connect_inds:
        for ind in p.individuals:
            fam, id = ind.label
            # Actually make the references instead of just pointing at strings
            ind.father = p[(fam, ind.father)] if ind.father != '0' else None
            ind.mother = p[(fam, ind.mother)] if ind.mother != '0' else None

            ind.register_with_parents()

    # Place individuals into pedigrees
    pedigrees = {}
    for ind in p.individuals:
        if ind.label[0] not in pedigrees:
            pedigrees[ind.label[0]] = []

        pedigrees[ind.label[0]].append(ind)

    for pedigree_label, ped_inds in list(pedigrees.items()):
        ped = Pedigree(label=pedigree_label)

        if isinstance(population_handler, Callable):
            population_handler(ped)
        
        for ind in ped_inds:
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
    :param pedigrees:   data to update
    :param csvfile:     the filename of the file containing phenotypes.
    :param delimiter:   the field delimiter for the file
    :param missingcode: the code for missing values
    :type pedigrees: PedigreeCollection
    :type csvfile: string
    :type missingcode: string

    :rtype: void
    """
    with open(csvfile) as f:
        header = f.readline().strip().split(delimiter)
        for line in f:
            # Match columns to their column name
            d = dict(list(zip(header, line.strip().split(delimiter))))
            for k, v in list(d.items()):
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
    ''' 
    Writes pedigree to a LINKAGE formatted pedigree file 

    :param pedigrees: Data to write
    :param filename: filename to write to
    :param missingcode: code to use for missing values
    :param delim: output field separator 

    :rtype: void
    '''
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

    with open(filename, 'w') as ofile:
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
    chrom1 = [1,1,1]
    chrom2 = [2,2,2]

    These are returned in the a list in the form [(chroma, chromb), (chroma, chromb)...]


    :param chromosomes: genotype data
    :param data: The alleles to be turned into genotypes
    :param missing_code: value representing a missing allele

    :type chromosomes: list of ChromosomeTemplate
    :type missing_code: string
    :returns: A list of 2-tuples of Alleles objects
    '''
    Chromobj = Alleles

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
