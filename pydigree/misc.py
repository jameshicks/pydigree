#!/usr/bin/env python

import random

from population import Population,Pedigree
from individual import Individual

def read_ped(filename,population=None,affected_labels=None):
    """
    Reads a plink format pedigree file, ie:
        familyid indid father mother sex whatever whatever whatever
    into a pydigree pedigree object, with optional population to
    assign to pedigree members. If you don't provide a population
    you can't simulate genotypes!
    """
    p = Pedigree()
    if not affected_labels:
        affected_labels = {'1':0,'2':1,
                           'A':1,'U':0,
                           'X': None,
                           '-9': None}
    def getph(ph):
        if ph in affected_labels:
            return affected_labels[ph]
        else: return None
    with open(filename) as f:
        for line in f:
            fam,id,fa,mo,sex,aff = line.strip().split(None,5)
            p[id] = Individual(population,id,fa,mo,sex)
            p[id].phenotypes['aff'] = getph(aff)
            p[id].pedigree = p
        for ind in p:
            # Actually make the references instead of just pointing at strings
            ind.father = p[ind.father] if ind.father != '0' else None
            ind.mother = p[ind.mother] if ind.mother != '0' else None
            ind.sex = {'1':0,'2':1}[ind.sex]
            if population: population.register_individual(p[id])
            ind.register_with_parents()
    return p


### Extra genetics functions
def is_missing_genotype(g):
    return g == tuple([0,0])

def ibs(g1,g2):
    """ Returns the number of alleles identical by state between two genotypes """
    g1,g2 = sorted(g1),sorted(g2)
    # Returns IBS state between two genotypes
    if is_missing_genotype(g1) or is_missing_genotype(g2):
        return None
    if g1 == g2: return 2
    g2s = set(g2)
    if g1[0] in g2s or g1[1] in g2s: return 1
    return 0



