"A pool of Chromosomes a simulation can draw from"

from itertools import chain

import numpy as np

from pydigree.recombination import recombine


def richards(A, C, M, B, T):
    '''
    Generates a population growth model function

    A: Lower Asymptope
    C: Upper Asymptope
    M: Maximum growth time
    B: Growth Rate
    T: Maximum growth position
    '''
    return lambda gen: A + (C / (1 + T * np.exp(-B * (gen - M)) ** (1 / T)))


class ChromosomePool(object):

    """
    A pool of Chromosomes a simulation can draw from 
    """

    def __init__(self, population=None, chromosomes=None, size=0):
        """
        Create the pool.

        :param population: 
        :param chromosomes: the set of chromosomes
        :param size: size of the pool
        :type size: int
        :type chromosomes: ChromosomeSet
        :type population: IndividualContainer
        """
        if population:
            self.chromosomes = population.chromosomes
        elif chromosomes:
            self.chromosomes = chromosomes
        self.pool = [[] * len(self.chromosomes)]
        self.n0 = size

    # Pool functions
    def size(self):
        """ Returns the size of the pool of available chromosomes """
        return len(self.pool[0])

    def initialize_pool(self, size=None):
        """ Initializes a pool of chromosomes for simulation """
        if self.n0 and not size:
            size = self.n0
        for i, q in enumerate(self.chromosomes):
            self.pool[i] = q.linkageequilibrium_chromosomes(2 * size)

    def fix(self, loc, value):
        ''' Sets all alleles at loc to value'''
        chromidx, posidx = loc
        p = self.pool[chromidx]
        for chrom in p:
            chrom[posidx] = value

    def iterate_pool(self, gensize):
        """
        Iterate pool simulates a generation of random mating
        between chromosomes instead of individuals. The pool of
        population chromosomes then contains chromosomes from the
        new generation.

        Arguements:
        gensize: The size of the next generation (rounded down to the integer)

        Returns: Nothing
        """
        # Generation sizes calculated from mathematical functions can have
        # non-integer values, which doesn't make much sense here.
        gensize = int(gensize)
        for i, c in enumerate(self.chromosomes):

            # Chromosomes have a 1/2 chance of being recombined with another
            def choose_chrom(pool, chrmap):
                """
                Get two random chromosomes, recombine them, return the result
                """
                # Since Alleles is a subclass of ndarray, numpy has been
                # treating pool as a multidimensional array. We'll generate
                # the indices ourself and get them that way. Eventually
                # I'll come back and fix the isinstancing of Alleles.
                qi, qw = np.random.randint(0, len(pool), 2)
                q, w = pool[qi], pool[qw]
                r = recombine(q, w, chrmap)
                return r

            newpool = [choose_chrom(self.pool[i], c.genetic_map)
                       for x in range(gensize)]
            self.pool[i] = newpool

    # Chromosome functions
    def chromosome(self, chromindex):
        """ 
        Get a random chromomsome from the pool
        :param chromindex: which chromosome are we looking for?
        :type chromindex: int

        :rtype: AlleleContainer
        """
        chidx = np.random.randint(0, len(self.pool[chromindex]))

        return self.pool[chromindex][chidx]

    def get_genotype_set(self):
        ''' Gives a full set of genotypes drawn from the chromosome pool '''
        return [[self.chromosome(i), self.chromosome(i)]
                for i, x in enumerate(self.chromosomes)]

    def evolve(self, growth_func, gens):
        ''' 
        Iterates the pool according to a popuation growth model.

        :param growth_func: A function that takes a generation number as an 
            argument and returns a generation size
        :param gens: number of generations to advance
        :type growth_func: Callable
        :type gens: int

        :rtype void: 
        '''
        for x in range(gens):
            self.iterate_pool(growth_func(x))

    @staticmethod
    def from_population(pop):
        """
        Creates a pool from an existing population.

        :param pop: Base population
        :type pop: Population

        :rtype: ChromosomePool
        """
        newpool = ChromosomePool(chromosomes=pop.chromosomes)
        newpool.n0 = len(pop.individuals) * 2
        for chridx, _ in enumerate(pop.chromosomes):
            poolchroms = chain.from_iterable(ind.genotypes[chridx]
                                             for ind in pop.individuals)
            thischrom = list(poolchroms)
            newpool.pool[chridx] = thischrom
        return newpool
