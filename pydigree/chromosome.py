#!/usr/bin/env python

from array import array
from pydigree._pydigree import choice_with_probs


class Chromosome(object):
    """
    Chromsome is a class that keeps track of marker frequencies and distances.
    Not an actual chromosome with genotypes, which you would find under
    Individual.

    Markers are currently diallelic and frequencies are given for minor
    alleles. Marker frequencies must sum to 1. Major allele frequency
    is then f = 1 - f_minor.

    linkageequilibrium_chromosome generates chromsomes that are generated from
    simulating all markers with complete independence (linkage equilibrium).
    This is not typically what you want: you won't find any LD for association
    etc. linkageequilibrium_chromosome is used for 'seed' chromosomes when
    initializing a population pool or when simulating purely family-based
    studies for linkage analysis.
    """
    def __init__(self):
        self.genetic_map = []
        self.frequencies = []
        self.typecode = 'B'  # Unsigned char

    def __str__(self):
        return 'Chromosome object: %d markers, %f cM' % \
            (len(self.frequencies), sum(self.genetic_map) * 100)

    def size(self):
        return self.genetic_map[-1] - self.genetic_map[0]

    def add_genotype(self, frequency, map_position):
        self.genetic_map.append(map_position)
        self.frequencies.append(frequency)

    def linkageequilibrium_chromosome(self):
        # Returns a randomly generated chromosome
        return array(self.typecode,
                     [choice_with_probs([1, 2], [1-f, f])
                      for f in self.frequencies])
