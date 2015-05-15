from itertools import izip

import numpy as np

from pydigree.exceptions import NotMeaningfulError
from pydigree.cyfuncs import fastfirstitem


class GenotypedChromosome(np.ndarray):

    ''' A class for holding genotypes '''
    def __new__(cls, data):
        obj = np.asarray(data).view(cls)
        return obj

    def __lt__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __gt__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __lte__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __gte__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    @property
    def missingcode(self):
        return 0 if np.issubdtype(self.dtype, np.integer) else ''

    @property
    def missing(self):
        ''' Returns a numpy array indicating which markers have missing data '''
        return self == self.missingcode

    def nmark(self):
        '''
        Return the number of markers represented by the
        GenotypedChromosome object
        '''
        return self.shape[0]


class SparseGenotypedChromosome(object):

    '''
    An object representing a set of haploid genotypes efficiently by 
    storing allele differences from a reference. Useful for manipulating
    genotypes from sequence data (e.g. VCF files)
    '''

    def __init__(self, data):
        data = np.array(data)
        self.dtype = data.dtype
        self.size = data.shape[0]
        self.non_refalleles = self.__array2nonref(data)
        self.missingindices = self.__array2missing(data)

    def __array2nonref(self, data):
        refcode = 0 if np.issubdtype(self.dtype, np.integer) else '0'
        return {i: x for i, x in enumerate(data)
                if x != refcode and x != ''}

    def __array2missing(self, data):
        return [i for i, x in enumerate(data) if x == self.missingcode]

    @property
    def missingcode(self):
        return 0 if np.issubdtype(self.dtype, np.integer) else ''

    @property
    def missing(self):
        ''' Returns a numpy array indicating which markers have missing data '''
        base = np.zeros(self.size, dtype=np.bool_)
        base[self.missingindices] = 1
        return base

    def __eq__(self, other):
        if self.size != other.size:
            raise ValueError('Trying to compare different-sized chromosomes')

        # SparseGenotypedChromosome saves differences from a reference,
        # so all reference sites are equal, and we mark everything True
        # to start, and go through and set any differences to False
        base = np.ones(self.size, dtype=np.bool_)

        nonref_a = self.non_refalleles.viewitems()
        nonref_b = other.non_refalleles.viewitems()

        # Get the alleles that are in nonref_a or nonref_b but not both
        neq_alleles = (nonref_a ^ nonref_b)
        neq_sites = fastfirstitem(neq_alleles)

        base[neq_sites] = 0

        return base

    def __neq__(self, other):
        return np.logical_not(self == other)

    def nmark(self):
        '''
        Return the number of markers (both reference and non-reference)
        represented by the SparseGenotypedChromosome object
        '''
        return self.size

    def todense(self):
        '''
        Returns a non-sparse GenotypeChromosome equivalent
        to a SparseGenotypedChromosome object.
        '''

        arr = np.zeros(self.size, dtype=np.uint8).astype(self.dtype)
        for loc, allele in self.non_refalleles.iteritems():
            arr[loc] = allele

        arr[self.missing] = self.missingcode

        return arr


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

    def __str__(self):
        return 'Chromosome %s: %s markers, %s cM' % \
            (self.label if self.label is not None else 'object',
             len(self.frequencies), max(self.genetic_map) if self.genetic_map else 0)

    def __iter__(self):
        return izip(self.labels, self.genetic_map, self.physical_map)

    def _iinfo(self):
        return izip(self.labels, self.genetic_map, self.physical_map,
                    self.frequencies)

    def nmark(self):
        ''' Returns the number of markers on the chromosome '''
        return len(self.genetic_map)

    def size(self):
        ''' Returns the size of the chromosome in centimorgans '''
        return self.genetic_map[-1] - self.genetic_map[0]

    def add_genotype(self, frequency=None, map_position=None, label=None, bp=None):
        try:
            frequency = float(frequency) if frequency is not None else -1
        except TypeError:
            raise ValueError('Invalid value for frequency %s' % frequency)
        self.genetic_map.append(map_position)
        self.frequencies = np.append(self.frequencies, frequency)
        self.physical_map.append(bp)
        self.labels.append(label)

    def set_frequency(self, position, frequency):
        """ Manually change an allele frequency """
        self.frequencies[position] = frequency

    def empty_chromosome(self, dtype=np.uint8):
        return GenotypedChromosome(np.zeros(self.nmark(), dtype=dtype))

    def linkageequilibrium_chromosome(self):
        """ Returns a randomly generated chromosome """
        if (self.frequencies < 0).any():
            raise ValueError('Not all frequencies are specified')
        r = np.random.random(self.nmark())
        r = np.array(r > self.frequencies, dtype=np.uint8) + 1
        return GenotypedChromosome(r)

    def linkageequilibrium_chromosomes(self, nchrom):
        """ Returns a numpy array of many randomly generated chromosomes """
        chroms = np.random.random((nchrom, self.nmark()))
        chroms = np.uint8((chroms > self.frequencies) + 1)
        return [GenotypedChromosome(r) for r in chroms]
