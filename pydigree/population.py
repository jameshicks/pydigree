#!/usr/bin/env python

# Packages we're going to use
import math
import numpy as np

# Other pydigree objects
from pydigree.individual import Individual

from pydigree.recombination import recombine
from pydigree.individualcontainer import IndividualContainer
from pydigree.genotypes import ChromosomeSet
from pydigree.simulation.mating import MatingStructure, RandomMating

missing_genotype = (0, 0)


def is_missing_genotype(g):
    return g == missing_genotype


# Population growth models.
# These are easy enough to supply your own if you wanted
def exponential_growth(p, r, t):
    """
    Models exponential growth over discrete generations.
    
    :param p: initial population
    :param r: growth rate
    :param t: number of generations
    :type p: numeric
    :type t: numeric
    :type r: numeric

    :returns: population size at time t
    :rtype: numeric
    """
    return p * math.exp(r * t)


def logistic_growth(p, r, k, t):
    """
    Models logistic growth over discrete generations.
    
    :param p: initial population
    :param r: growth rate
    :param k: final population
    :param t: number of generations

    :returns: population size at time t
    :rtype: numeric
    """
    return (p * k) / (p + (k - p) * math.exp(-r * t))


# Classes
class Population(IndividualContainer):
    # Methods for mapping types

    def __init__(self, intial_pop_size=0, name=None):
        self.chromosomes = ChromosomeSet()
        self.pool = None
        self.population = {}
        self.n0 = intial_pop_size
        self.name = name

    def __hash__(self):
        return id(self)

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
        for x in self.individuals:
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

    def add_founders(self, n):
        """
        Adds a number of founder individuals to the population

        :param n: number of individuals to add
        :type n: int

        :rtype: void
        """
        for _ in range(n):
            self.founder_individual(register=True)

    def update(self, other):
        '''
        Merges two datasets (i.e. performs Individual.update for each individual in the pedigree)

        Assumes unique individual IDs
        
        :param other: New data to merge in
        :type other: Population
        
        :return: void
        '''
        self.chromosomes = other.chromosomes
        self.clear_genotypes()

        selfids = {x.label for x in self.individuals}
        otherids = {x.label for x in other.individuals}
        overlap = set.intersection(selfids, otherids)

        if not overlap:
            return

        for x in overlap:
            self.population[x].update(other[x])

    def _getindividual(self, label):
        return self[label]

    @property
    def individuals(self):
        ''' Returns a list of individuals in the population '''
        return [x for x in self.population.values()]


    # Chromosome functions
    #
    #
    def add_chromosome(self, chrom):
        """ Adds a chromosome to the population """
        self.chromosomes.add_chromosome(chrom)

    def chromosome_count(self):
        """ Returns the number of chromosomes in this population """
        return len(self.chromosomes)

    # Random mating
    #
    #
    def mate(self, ind1, ind2, indlab, sex=None):
        """
        Creates an individual as the child of two specificied individual
        objects and randomly chooses a sex.

        :param ind1: The first parent
        :param ind2: The second parent
        :type ind1: Individual
        :type ind2: Individual
        :param indlab: ID label for the child
        :param sex: Sex of child, randomly chosen if not specified
        :type sex: {0,1}
        :return: An individual with ind1 and ind2 as parents
        :rtype: Individual
        """
        if sex is None:
            sex = np.random.choice([0, 1])
        child = Individual(self, indlab, ind1, ind2, sex)
        return child

    def advance_generation(self, gensize, mating=None):
        '''
        Simulates a generation of random mating.

        :param gensize: The size of the new generation
        :param mating: MatingScheme for the generation  
        :type gensize: numeric
        :type mating: MatingScheme
        '''
        if mating is None:
            mating = RandomMating()


        progeny = mating.next_generation(self, gensize)

        self.population = {x.label : x for x in progeny}

    def founder_individual(self, register=True, sex=None):
        "Creates a new founder individual and adds to the population"

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
        ''' 
        Gives genotypes to each founder in the population with chromosomes 
        from the chromosome pool. If there is no pool, genotypes are generated
        under linkage equilibrium
        '''
        for ind in self.individuals:
            if not self.pool:
                genotypes = self.get_linkage_equilibrium_genotypes()
            else: 
                genotypes = self.pool.get_genotype_set()

            ind.genotypes = genotypes

    def get_genotypes(self):
        ''' 
        Causes each Individual object in the pedigree to request genotypes
        from its parents
        '''
        for x in self.individuals:
            x.get_genotypes()

    def get_linkage_equilibrium_genotypes(self):
        '''
        Returns a set of genotypes for an individual in linkage equilibrium
        '''
        return [[c.linkageequilibrium_chromosome(),
                 c.linkageequilibrium_chromosome()]
                for c in self.chromosomes]


