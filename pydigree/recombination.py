#!/usr/bin/env python

import random

from array import array
from bisect import bisect_left


def recombine(chr1, chr2, map):
    newchrom = _recombine_haldane(chr1, chr2, map)
    if isinstance(chr1, array) and isinstance(chr2, array):
        if chr1.typecode != chr2.typecode:
            raise ValueError('Chromosomes have two different typecodes!')
        newchrom = array(chr1.typecode, newchrom)
    return newchrom


def _recombine_haldane(chr1, chr2, map):
    # The map is sorted list, and the last item will always be largest.
    maxmap = map[-1]

    newchrom = [None]*len(map)
    # Randomly pick a chromosome to start from
    flipped = random.choice([True, False])
    last_crossover_index = 0
    crossover_position = 0
    while True:
        # Get from the next chromosome
        flipped = not flipped
        c = chr1 if flipped else chr2

        # Find the next crossover point
        crossover_position += random.expovariate(0.01)

        if crossover_position > maxmap:
            # We've reached the end of our journey here.
            newchrom[last_crossover_index:] = c[last_crossover_index:]
            break 

        # Find the next crossover point in the chromosome by binary search
        nextidx = bisect_left(map, crossover_position)
        newchrom[last_crossover_index:nextidx] = c[last_crossover_index:nextidx]
        
        # Get ready to do it all over again
        last_crossover_index = nextidx

    return newchrom