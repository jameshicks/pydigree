#!/usr/bin/env python

from array import array
from bisect import bisect_left

import numpy as np

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
    nmark = len(map)

    newchrom = []
    # Randomly pick a chromosome to start from
    # np.random.randint works on a half open interval, so the upper bound 
    # specified is 2. We'll get zeros and ones out of it.
    flipped = np.random.randint(0,2)

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
            newchrom.extend(c[last_crossover_index:])
            break 

        # Find the next crossover point in the chromosome by binary search
        nextidx = bisect_left(map, crossover_position, last_crossover_index, nmark)
        newchrom.extend(c[last_crossover_index:nextidx])
        
        # Get ready to do it all over again
        last_crossover_index = nextidx

    return newchrom