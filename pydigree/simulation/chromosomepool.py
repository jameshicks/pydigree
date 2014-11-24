import random

import numpy as np

from pydigree.recombination import recombine

class ChromosomePool(object):
    def __init__(self, population=None, chromosomes=None, size=0):
        if population:
            self.chromosomes = population.chromosomes
        elif chromosomes:
            self.chromosomes = chromosomes
        self.pool = [ [] * len(self.chromosomes) ]
        self.n0 = size
        self.generations = []

    # Pool functions
    def size(self):
        """ Returns the size of the pool of available chromosomes """
        return len(self.pool[0])

    def initialize_pool(self, size=None):
        """ Initializes a pool of chromosomes for simulation """
        if self.n0 and not size:
            size = self.n0
        for i, q in enumerate(self.chromosomes):
            self.pool[i] = q.linkageequilibrium_chromosome(2*size)
        self.generations.append(size)

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
                q, w = random.choice(pool), random.choice(pool)
                r = recombine(q,w, chrmap)
                return r

            newpool = [choose_chrom(self.pool[i], c.genetic_map)
                       for x in xrange(gensize)]
            self.pool[i] = newpool

        self.generations.append(gensize)
    
    # Chromosome functions
    def chromosome(self, chromindex):
        # Get a random chromomsome
        return random.choice(self.pool[chromindex])

    def genotypes(self):
        ''' Gives a full set of genotypes drawn from the chromosome pool '''
        return [ [self.chromosome(i), self.chromosome(i)] 
                 for i,x in enumerate(self.chromosomes)]
