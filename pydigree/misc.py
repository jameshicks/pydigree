#!/usr/bin/env python

import random

from population import Population
from individual import Individual
from pedigree import Pedigree
from pedigreecollection import PedigreeCollection
from chromosome import Chromosome



def read_ped(filename, population=None, delimiter=None, affected_labels=None):
    """
    Reads a plink format pedigree file, ie:
        familyid indid father mother sex whatever whatever whatever
    into a pydigree pedigree object, with optional population to
    assign to pedigree members. If you don't provide a population
    you can't simulate genotypes!

    Arguements
    -----
    filename: The file to be read
    population: The population to assign individuals to
    delimiter: a string defining the field separator, default: any whitespace
    affected_labels: The labels that determine affection status.

    Returns: An object of class PedigreeCollection
    """
    p = Pedigree()
    pc = PedigreeCollection()
    sex_codes = {'1': 0, '2': 1, 'M': 0, 'F': 1}
    if not affected_labels:
        affected_labels = {'1': 0, '2': 1,
                           'A': 1, 'U': 0,
                           'X': None,
                           '-9': None}

    def getph(ph):
        if ph in affected_labels:
            return affected_labels[ph]
        else:
            return None
    with open(filename) as f:
        for line in f:
            fam, id, fa, mo, sex, aff = line.strip().split(delimiter, 5)
            # Give a special id for now, to prevent overwriting duplicated
            # ids between families
            id = (fam, id)
            p[id] = Individual(population, id, fa, mo, sex)
            p[id].phenotypes['aff'] = getph(aff)
            p[id].pedigree = p
        for ind in p:
            fam, id = ind.id
            # Actually make the references instead of just pointing at strings
            ind.father = p[(fam, ind.father)] if ind.father != '0' else None
            ind.mother = p[(fam, ind.mother)] if ind.mother != '0' else None
            ind.sex = sex_codes[ind.sex]
            if population:
                population.register_individual(p[id])
            ind.register_with_parents()
        for pedigree_label in set(ind.id[0] for ind in p):
            ped = Pedigree(pedigree_label)
            thisped = [x for x in p if x.id[0] == pedigree_label]
            for ind in thisped:
                ind.id = ind.id[1]
                ped[ind.id] = ind
            pc[pedigree_label] = ped
    return pc


def read_map(mapfile):
    last_chr = None
    chroms = []
    chromosome = Chromosome()
    with open(mapfile) as f:
        for line in f:
            line = line.strip().split()
            chr, label, cm, pos = line
            if chr != last_chr:
                chroms.append(chromosome)
                chromosome = Chromosome()
            chromosome.add_genotype(None, cm, label=label, bp=pos)
    return chroms


def read_phenotypes(pedigrees, csvfile, delimiter=',', missingcode='X'):
    """
    Reads a csv with header
    famid,ind,phen,phen,phen,phen etc etc

    Arguements
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
                if k in set(['famid', 'ind']):
                    continue
                fam, ind = d['famid'], d['ind']
                pedigrees[fam][ind].phenotypes[k] = v


### Extra genetics functions
def is_missing_genotype(g):
    return g == tuple([0, 0])


def ibs(g1, g2):
    """
    Returns the number of alleles identical by state between two genotypes
    Arguements: Two tuples
    Returns: an integer
    """
    g1, g2 = sorted(g1), sorted(g2)
    # Returns IBS state between two genotypes
    if is_missing_genotype(g1) or is_missing_genotype(g2):
        return None
    if g1 == g2:
        return 2
    g2s = set(g2)
    if g1[0] in g2s or g1[1] in g2s:
        return 1
    return 0
