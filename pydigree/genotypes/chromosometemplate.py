from bisect import bisect_right
from itertools import izip
import numpy as np

from pydigree.genotypes import Alleles, SparseAlleles
from pydigree.io.genomesimla import read_gs_chromosome_template

class ChromosomeTemplate(object):

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

    def __init__(self, label=None):
        # Chromosome name
        self.label = label
        # A list of floats that represent the position of the marker in cM
        self.genetic_map = []
        # A list of integers that doesnt do anything. Just for decoration
        self.physical_map = []
        # A list of floats representing minor allele frequencies
        self.frequencies = np.array([])
        # List of marker names
        self.labels = []
        # Reference Alleles
        self.reference = []
        # Alternates
        self.alternates = []

    def __str__(self):
        return 'ChromosomeTemplate object %s: %s markers, %s cM' % \
            (self.label if self.label is not None else 'object',
             len(self.frequencies),
             max(self.genetic_map) if self.genetic_map else 0)

    @property
    def outputlabel(self):
        ''' The label outputted when written to disk '''
        if self.label:
            return self.label
        else:
            return 0

    def __iter__(self):
        return izip(self.labels, self.genetic_map, self.physical_map)

    def _iinfo(self):
        return izip(self.labels, self.genetic_map, self.physical_map,
                    self.frequencies)

    @staticmethod
    def from_genomesimla(filename):
        return read_gs_chromosome_template(filename)

    def nmark(self):
        ''' Returns the number of markers on the chromosome '''
        return len(self.genetic_map)

    def size(self):
        ''' Returns the size of the chromosome in centimorgans '''
        return self.genetic_map[-1] - self.genetic_map[0]

    def add_genotype(self, frequency=None, map_position=None, label=None,
                     bp=None, reference=None, alternates=None):
        try:
            frequency = float(frequency) if frequency is not None else -1
        except TypeError:
            raise ValueError('Invalid value for frequency %s' % frequency)
        self.genetic_map.append(map_position)
        self.frequencies = np.append(self.frequencies, frequency)
        self.physical_map.append(bp)
        self.labels.append(label)
        self.reference.append(reference)
        self.alternates.append(alternates)

    def set_frequency(self, position, frequency):
        """ Manually change an allele frequency """
        self.frequencies[position] = frequency

    def empty_chromosome(self, dtype=np.uint8):
        return Alleles(np.zeros(self.nmark(), dtype=dtype))

    def closest_marker(self, position, map_type='physical'):
        "Returns the index of the closest marker to a position"
        if map_type == 'physical':
            mp = self.physical_map
        elif map_type == 'genetic':
            mp = self.genetic_map
        else:
            raise ValueError("Map type must be 'physical' or 'genetic'")

        # Find the index in mp with value lte to position
        left_idx = bisect_right(mp, position) - 1

        if left_idx == self.nmark() - 1:
            # If we're already at the last marker, we know to just return
            # left_idx
            return left_idx

        right_idx = left_idx + 1

        right_pos = mp[right_idx]
        left_pos = mp[left_idx]

        if abs(right_pos - position) < abs(left_pos - position):
            return right_idx
        else:
            return left_idx

    def linkageequilibrium_chromosome(self, sparse=False):
        """ Returns a randomly generated chromosome """
        if (self.frequencies < 0).any():
            raise ValueError('Not all frequencies are specified')
        r = np.random.random(self.nmark())
        r = np.array(r < self.frequencies, dtype=np.int8) + 1

        if sparse:
            return SparseAlleles(r, refcode=1, template=self)
        else:
            return Alleles(r, template=self)

    def linkageequilibrium_chromosomes(self, nchrom):
        """ Returns a numpy array of many randomly generated chromosomes """
        chroms = np.random.random((nchrom, self.nmark()))
        chroms = np.int8((chroms < self.frequencies) + 1)
        return [Alleles(r) for r in chroms]
