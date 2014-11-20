import random

import numpy as np
from bitarray import bitarray

def __bitarray2ndarray(ba):
    return np.fromstring(ba.unpack(), dtype=np.int)

class ChromosomePool(object):
    def __init__(self, population=None, chromosomes=None, n0=0):
        if population:
            self.generations = []
            self.chromosome
        pass

    @property
    def chromosomes(self):
        if self.chromosomes:
            return self.chromosomes
        elif self.population:
            return self.population.chromosomes
        else:
            return []

    # Pool functions
    def size(self):
        """ Returns the size of the pool of available chromosomes """
        return len(self.pool[0])

    def initialize_pool(self, size=None):
        """ Initializes a pool of chromosomes for simulation """
        if self.n0 and not size:
            size = self.n0
        self.pool = [None] * len(self.chromosomes)
        for i, q in enumerate(self.chromosomes):
            self.pool[i] = [q.linkageequilibrium_chromosome(bitarray=True)
                            for x in xrange(2 * size)]
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
                return bitarray(recombine(q, w, chrmap))

            newpool = [choose_chrom(self.pool[i], c.genetic_map)
                  for x in xrange(gensize)]
            self.pool[i] = newpool

        self.generations.append(gensize)
    
    # Chromosome functions
    def chromosome(self, chromindex):
        # Get a random bitarray chromomsome
        c = random.choice(self.pool[chromindex])
        # Convert to a numpy array and return it
        return __bitarray2ndarray(c) + 1
        
    def genotypes(self):
        ''' Gives a full set of genotypes drawn from the chromosome pool '''
        return [ [self.chromosome(i), self.chromosome(i)] 
                 for i,x in enumerate(self.chromosomes)]
