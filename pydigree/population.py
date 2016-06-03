#!/usr/bin/env python

# Packages we're going to use
import math
import numpy as np
import pandas as pd
from itertools import chain

# Abstract base class for population and pedigree
from collections import MutableMapping

# Other pydigree objects
from pydigree.individual import Individual
from pydigree.common import *
from pydigree.recombination import recombine
from pydigree.individualcontainer import IndividualContainer

missing_genotype = (0, 0)


def is_missing_genotype(g):
    return g == missing_genotype


# Population growth models.
# These are easy enough to supply your own if you wanted
def exponential_growth(p, r, t):
    """
    Models exponential growth over discrete generations.
    p: initial population
    r: growth rate
    t: number of generations
    """
    return p * math.exp(r * t)


def logistic_growth(p, r, k, t):
    """
    Models logistic growth over discrete generations.
    p: initial population
    r: growth rate
    k: final population
    t: number of generations
    """
    return (p * k) / (p + (k - p) * math.exp(-r * t))


# Classes
class Population(MutableMapping, IndividualContainer):
    # Methods for mapping types

    def __init__(self, intial_pop_size=0, name=None):
        self.chromosomes = []
        self.pool = None
        self.population = {}
        self.n0 = intial_pop_size
        self.name = name

    def __hash__(self):
        return id(self)

    # Things I have to implement for the ABC
    #
    #
    def __iter__(self):
        return iter(self.population.values())

    def __getitem__(self, key):
        return self.population[key]

    def __contains__(self, item):
        return item in self.population.values()

    def __len__(self):
        return len(self.population)

    def __setitem__(self, key, value):
        self.population[key] = value

    def __delitem__(self, key):
        del self.population[key]

    def size(self):
        """ Returns the number of individuals in the population. """
        return len(self.population)

    def remove_ancestry(self):
        """ Makes every individual in the population a founder """
        for x in self:
            x.remove_ancestry()

    # Adding and removing people
    #
    #
    def register_individual(self, ind):
        ''' Adds an individual to the population '''
        if ind.label in self.population:
            raise ValueError('ID %s already in population!' % ind.label)
        self.population[ind.label] = ind

    def remove_individual(self, ind):
        ''' Removes an individual from the population '''
        del self[ind.label]


    def update(self, other):
        '''
        Merges two datasets (i.e. performs Individual.update for each individual in the pedigree)

        Assumes unique individual IDs
        '''
        self.chromosomes = other.chromosomes
        self.clear_genotypes()

        selfids = {x.label for x in self}
        otherids = {x.label for x in other.individuals}
        overlap = set.intersection(selfids, otherids)

        if not overlap:
            return

        for x in overlap:
            self.population[x].update(other._getindividual(x))

    def _getindividual(self, label):
        return self[label]

    @property
    def individuals(self):
        ''' Returns a list of individuals in the population '''
        return [x for x in self]



    # Chromosome functions
    #
    #
    def add_chromosome(self, chrom):
        """ Adds a chromosome to the population """
        self.chromosomes.append(chrom)

    def chromosome_count(self):
        """ Returns the number of chromosomes in this population """
        return len(self.chromosomes)

    # Random mating
    #
    #
    def mate(self, ind1, ind2, id):
        """
        Creates an individual as the child of two specificied individual
        objects and randomly chooses a sex.
        """
        child = Individual(self, id, ind1, ind2, np.random.choice([0, 1]))
        return child

    def random_mating_generation(self, gensize):

        def rand_mate(pop, name=None):
            name = str(name) if not self.name else self.name + str(name)
            return self.mate(np.random.choice(self.males()),
                             np.random.choice(self.females()),
                             name)

        newpop = {}
        for i in xrange(gensize):
            newind = rand_mate(self, name=i)
            newpop[newind.label] = newind
        self.population = newpop

    def founder_individual(self, register=True, sex=None):
        ''' Creates a new founder individual and adds to the population '''
        if sex is not None:
            sex = sex.lower()
        sexd = {'m': 0, 'f': 1, None: np.random.choice([0, 1])}
        i = Individual(self, self.size(), None, None, sexd[sex])
        if register:
            self.register_individual(i)
        return i

    # Genotype Functions
    #
    #
    def get_founder_genotypes(self):
        ''' Gives genotypes to each founder in the population
        with chromosomes from the chromosome pool '''
        if not self.pool:
            raise ValueError('Nothing in chromosome pool')
        return self.pool.get_genotype_set()

    def get_genotypes(self):
        ''' 
        Causes each Individual object in the pedigree to request genotypes
        from its parents
        '''
        for x in self:
            x.get_genotypes()

    def get_linkage_equilibrium_genotypes(self):
        '''
        Returns a set of genotypes for an individual
        in linkage equilibrium
        '''
        return [[c.linkageequilibrium_chromosome(),
                 c.linkageequilibrium_chromosome()]
                for c in self.chromosomes]

    def ld(self, locusa, locusb, allelea=None, alleleb=None, method='r2'):
        """
        Returns a measure of linkage disquilibirum between alleles at two
        loci, for methods in {D, Dprime, r2}. If alleles aren't specified,
        the major alleles at each locus are chosen.

        LD metrics are functions of D, the measure of deviance from
        independence in sampling alleles. D = x11 - p1q1, where x11 is the
        frequency of allele 1 at locus p occuring with allele 1 at locus q.

        Lewontin's D' is a normalization of D, since D is dependent on
        frequency.

        D' = D / Dmax, where:
        Dmax = min(p1*q1, p2*q2) if D > 0
        Dmax = min(p2*q1,p1*q2) if D < 0

        r2 is a third metric, equal to
        D**2 / (p1 * p2 * q1 * q2)
        r2 ranges from 0 to 1

        Alleles are in linkage equilibrium if any of these metrics are 0.

        Arguments
        -----
        locusa:  The first locus
        locusb:  The second locus
        allelea: The allele at the first locus
                 (the major allele if not specified)
        alleleb: The allele at the second locus
                 (the major allele if not specified)
        method:  The metric of linkage disequilibrium
                 (r2 [default], Dprime, or D)

        Returns: a double
        """
        raise NotImplementedError('LD requires phase known data')
        method = method.lower()
        if method not in set(['r2', 'dprime', 'd']):
            raise ValueError("LD metric must be one of r2, Dprime or D")
        if not allelea:
            allelea = self.major_allele(locusa)
        if not alleleb:
            alleleb = self.major_allele(locusb)
        pop = [x for x in self if
               (not is_missing_genotype(x.get_genotype(locusa))) and
               (not is_missing_genotype(x.get_genotype(locusb)))]
        # Allele 1
        p1 = self.allele_frequency(locusa, allelea)
        p2 = 1 - p1
        # Allele 2
        q1 = self.allele_frequency(locusa, allelea)
        q2 = 1 - q1
        x = sum(1 for x in pop
                if x.has_allele(locusa, allelea) and
                x.has_allele(locusb, alleleb))
        D = x - p1 * q1
        # All of these metrics are normalizations of D by dividing it by
        # some expression. If D==0 all of them will equal 0.
        if method == 'd':
            return D
        elif method == 'dprime':
            if D > 0:
                Dmax = min(p1 * q1, p2 * q2)
            else:
                Dmax = min(p2 * q1, p1 * q2)
            return D / Dmax
        elif method == 'r2':
            return (D ** 2) / (p1 * p2 * q1 * q2)


