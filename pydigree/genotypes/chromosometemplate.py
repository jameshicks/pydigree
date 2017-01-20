from bisect import bisect_right

import numpy as np

from pydigree.genotypes import Alleles, SparseAlleles
from pydigree.io.genomesimla import read_gs_chromosome_template

class ChromosomeSet(object):
    """
    An object representing the full complement of variants in a population
    """

    def __init__(self):
        self.chroms = []

    def __iter__(self):
        for c in self.chroms:
            yield c

    def __getitem__(self, idx):
        return self.chroms[idx]

    def __len__(self):
        return len(self.chroms)

    def add_chromosome(self, template):
        """
        Add a chromosome to the set.

        :param template: the chromosome to add
        :type template: ChromosomeTemplate
        """
        self.chroms.append(template)

    def finalize(self):
        "Finalize each template in the set"
        for c in self.chroms:
            c.finalize()

    def frequency(self, chrom, variant):
        """
        Get the frequency of a variant

        :param chrom: index of chromosome
        :param variant: index of marker on the chromosome
        :type chrom: int
        :type variant: int
        """
        return self.chroms[chrom].frequencies[variant]

    def physical_map(self, chrom, variant):
        """
        Get the physical position of a variant

        :param chrom: index of chromosome
        :param variant: index of marker on the chromosome
        :type chrom: int
        :type variant: int
        """
        return self.chroms[chrom].physical_map[variant]

    def marker_label(self, chrom, variant):
        """
        Get the label of a variant

        :param chrom: index of chromosome
        :param variant: index of marker on the chromosome
        :type chrom: int
        :type variant: int
        """
        return self.chroms[chrom].labels[variant]

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
        self.final = False
        # Chromosome name
        self.label = label
        # A list of floats that represent the position of the marker in cM
        self.genetic_map = []
        # A list of integers that doesnt do anything. Just for decoration
        self.physical_map = []
        # A list of floats representing minor allele frequencies
        self.frequencies = []
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
        return zip(self.labels, self.genetic_map, self.physical_map)

    def _iinfo(self):
        return zip(self.labels, self.genetic_map, self.physical_map,
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
        if self.final:
            raise ValueError('Chromosome already finalized')
        try:
            frequency = float(frequency) if frequency is not None else -1
        except TypeError:
            raise ValueError('Invalid value for frequency %s' % frequency)
        self.genetic_map.append(map_position if map_position else 0)
        self.frequencies.append(frequency if frequency else 0)
        self.physical_map.append(bp if bp else 0)
        self.labels.append(label)
        self.reference.append(reference)
        self.alternates.append(alternates)

    def set_frequency(self, position, frequency):
        """ Manually change an allele frequency """
        self.frequencies[position] = frequency

    def empty_chromosome(self, dtype=np.uint8, sparse=False, refcode=None):
        if sparse:
            
            return SparseAlleles(size=self.nmark(), template=self, refcode=0)
        else:
            return Alleles(np.zeros(self.nmark(), dtype=dtype), template=self)

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

    def finalize(self):
        self.final = True
        self.frequencies = np.array(self.frequencies)
        self.physical_map = np.array(self.physical_map, dtype=np.int)
        self.genetic_map = np.array(self.genetic_map)

    def linkageequilibrium_chromosome(self, sparse=False):
        """ Returns a randomly generated chromosome """
        if (self.frequencies < 0).any():
            raise ValueError('Not all frequencies are specified')
        r = np.random.random(self.nmark())
        r = np.array(r < self.frequencies, dtype=np.int8) + 1

        if sparse:
            return SparseAlleles(r-1, refcode=0, template=self)
        else:
            return Alleles(r, template=self)

    def linkageequilibrium_chromosomes(self, nchrom):
        """ Returns a numpy array of many randomly generated chromosomes """
        chroms = np.random.random((nchrom, self.nmark()))
        chroms = np.int8((chroms < self.frequencies) + 1)
        return [Alleles(r) for r in chroms]
