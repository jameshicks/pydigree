#!/usr/bin/env python

import argparse
import random
import numpy as np
from itertools import izip

parser = argparse.ArgumentParser()
parser.add_argument('--size_mb', dest='size', help='Chromosome size (mb)', type=lambda x: int(float(x) * 10**6), default=150000000)
parser.add_argument('--nmark', help='Number of markers', type=int, default=int(1e5))
parser.add_argument('--min_freq', dest='minfreq', help='Minimum allele frequency', type=float, default=0.1)
parser.add_argument('--name', dest='name', help='Chromosome label', default='1')
args = parser.parse_args()

def random_float_range(low, high):
    ''' Returns a random float that meets the constraint low < x < high '''
    while True:
        r = random.random()
        if low < r < high:
            yield r

# Generate positions
positions_bp = sorted(np.random.randint(0, high=args.size+1, size=args.nmark))

# Generate allele freqs
freqs = random_float_range(args.minfreq, 1-args.minfreq)

for i, data in enumerate(izip(positions_bp, freqs)):
    bp, freq = data
    cm = bp / 1e6
    print args.name, 'marker{}'.format(i), cm, bp, freq, 1-freq
