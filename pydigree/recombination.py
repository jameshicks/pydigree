"""
Functions for recombining haploid chromosomes
"""

from bisect import bisect_left
from pydigree.genotypes import AlleleContainer

import numpy as np


def recombine(chr1, chr2, genetic_map):
    """
    Takes two chromatids and returns a simulated one by an exponential process
    
    :param chr1: first chrom
    :param chr2: second chrom
    :param genetic_map: map positions (in centiMorgans)
    :type chr1: AlleleContainer
    :type chr2: AlleleContainer
    :type genetic_map: sequence of floats

    :returns: Recombined chromosome
    :rtype: AlleleContainer
    """
    if not isinstance(chr1, AlleleContainer): 
        raise ValueError(
            'Invalid chromosome type for recombination: {}'.format(type(chr1))) 

    if type(chr1) is not type(chr2):
        raise ValueError("Can't mix chromosome types in recombination")

    if chr1.dtype != chr2.dtype:
        raise ValueError('Chromosomes have different data types')

    # An optimization for gene dropping procedures on IBD states.
    # If there is only one marker, choose one at random, and return that.
    # There's no need for searching through the map to find crossover points
    if len(genetic_map) == 1:
        return chr1 if np.random.randint(0, 2) else chr2

    newchrom = _recombine_haldane(chr1, chr2, genetic_map)
    return newchrom

def _recombine_haldane(chr1, chr2, genetic_map):    

    # The map is sorted list, and the last item will always be largest.
    maxmap = genetic_map[-1]
    nmark = len(genetic_map)
    
    # Return a new genotype container with the same specs as chr1
    newchrom = chr1.empty_like() 

    # Randomly pick a chromosome to start from
    # np.random.randint works on a half open interval, so the upper bound
    # specified is 2. We'll get zeros and ones out of it.
    flipped = np.random.randint(0, 2)

    last_crossover_index = 0
    crossover_position = 0
    while True:
        # Get from the next chromosome
        flipped = not flipped
        c = chr1 if flipped else chr2

        # Find the next crossover point
        # np.random.exponential is parameterized with the RECIPROCAL of the
        # rate parameter. With random.expovariate I would have used (0.01),
        # here I supply 100 as an argument.
        crossover_position += np.random.exponential(100)

        if crossover_position > maxmap:
            # We've reached the end of our journey here.
            newchrom.copy_span(c, last_crossover_index, None)
            break

        # Find the next crossover point in the chromosome by binary search
        nextidx = bisect_left(genetic_map, crossover_position, 
                              last_crossover_index, nmark)
        
        newchrom.copy_span(c, last_crossover_index, nextidx)

        # Get ready to do it all over again
        last_crossover_index = nextidx

    return newchrom
