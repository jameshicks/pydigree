from itertools import izip

import numpy as np

from pydigree.exceptions import NotMeaningfulError
from pydigree.cyfuncs import fastfirstitem
from pydigree.io.genomesimla import read_gs_chromosome_template
from pydigree.common import spans, all_same_type


class AlleleContainer(object):

    " A base class for the interface *Alleles object must implement"

    def empty_like(self):
        raise NotImplementedError

    def copy_span(self, template, start, stop):
        raise NotImplementedError

    def dtype(self):
        raise NotImplementedError

    def __eq__(self, other):
        raise NotImplementedError


class LabelledAlleles(AlleleContainer):

    def __init__(self, spans=None, chromobj=None, nmark=None):
        if not (chromobj or nmark):
            raise ValueError('One of chromobj or nmark must be specified')
        self.spans = spans if spans is not None else []
        self.chromobj = chromobj
        self.nmark = nmark if self.chromobj is None else self.chromobj.nmark()

    def __eq__(self, other):
        if not isinstance(other, LabelledAlleles):
            return False
        return all(x == y for x, y in izip(self.spans, other.spans))

    def empty_like(self):
        return LabelledAlleles([], chromobj=self.chromobj, nmark=self.nmark)

    @property
    def dtype(self):
        return type(self)

    @staticmethod
    def founder_chromosome(ind, chromidx, hap, chromobj=None, nmark=None):
        n = nmark if not chromobj else chromobj.nmark()
        spans = [InheritanceSpan(ind, chromidx, hap, 0, n)]
        return LabelledAlleles(spans=spans, chromobj=chromobj, nmark=nmark)

    def add_span(self, new_span):
        if any(new_span.stop < x.stop for x in self.spans):
            raise ValueError('Overwriting not supported for LabelledAlleles')
        self.spans.append(new_span)

    def copy_span(self, template, copy_start, copy_stop):
        if not isinstance(template, LabelledAlleles):
            raise ValueError(
                'LabelledAlleles can only copy from other LabelledAlleles')

        if copy_stop is None:
            copy_stop = self.nmark

        for span in template.spans:
            template_start, template_stop = span.interval
            # There are three possible things that can happen here:
            # 1) The template span is before the copy region, so we ignore it
            # 2) The template span is after the copy region, so we're done
            # 3) The template span overlaps the copy region.
            #
            # We're going to bail out early if the end of the copy region is in
            # the template span, so #2 should never be encountered. I've left
            # The test in, to cover the bases. #1 is trivial to check for.
            # Scenario #3 is where actual copying occurs.
            if template_stop < copy_start:
                # Happens before the copy region we're looking for
                continue

            if template_start > copy_stop:
                # Happens after the end of the copy region so we're done
                return

            if (template_start <= copy_start <= template_stop
                    or template_start <= copy_stop <= template_stop):
                # This this span of template overlaps the copy region, so
                # we have something to do

                # We take the maxes of the starts and the mins of the stops
                this_copystart = max(copy_start, template_start)
                this_copystop = min(copy_stop, template_stop)

                copy_span = InheritanceSpan(span.ancestor,
                                            span.chromosomeidx,
                                            span.haplotype,
                                            this_copystart,
                                            this_copystop)

                self.add_span(copy_span)

                # The early bailout if we're done copying
                if copy_stop <= template_stop:
                    return

    def delabel(self):
        # Check to make sure all the founders are delabeled
        if not all_same_type(self.spans, InheritanceSpan):
            for span in self.spans:
                if isinstance(span.ancestral_chromosome, LabelledAlleles):
                    raise ValueError('Ancestral chromosome {} {} {}'
                                     'has not been delabeled'.format(
                                         self.individual,
                                         self.chromosomeidx,
                                         self.haplotype))

        nc = self.spans[0].ancestral_chromosome.empty_like()
        for span in self.spans:
            nc.copy_span(span.ancestral_chromosome, span.start, span.stop)
        return nc


class InheritanceSpan(object):
    __slots__ = ['ancestor', 'chromosomeidx', 'haplotype', 'start', 'stop']

    def __init__(self, ancestor, chromosomeidx, haplotype, start, stop):
        self.ancestor = ancestor
        self.chromosomeidx = chromosomeidx
        self.haplotype = haplotype
        self.start = start
        self.stop = stop

    def __repr__(self):
        return 'InheritanceSpan{}'.format(self.to_tuple())

    def __eq__(self, other):
        return (self.ancestor == other.ancestor and
                self.chromosomeidx == other.chromosomeidx and
                self.haplotype == other.haplotype and
                self.start == other.start and
                self.stop == other.stop)

    @property
    def interval(self):
        return self.start, self.stop

    def to_tuple(self):
        return (self.ancestor, self.chromosomeidx, self.haplotype,
                self.start, self.stop)

    @property
    def ancestral_chromosome(self):
        return self.ancestor.genotypes[self.chromosomeidx][self.haplotype]


class AncestralAllele(object):
    __slots__ = ['ancestor', 'haplotype']

    def __init__(self, anc, hap):
        self.ancestor = anc
        self.haplotype = hap

    def __repr__(self):
        return 'AncestralAllele: {}: {}'.format(self.ancestor, self.haplotype)

    def __eq__(self, other):
        return (self.ancestor is other.ancestor and
                self.haplotype is other.haplotype)


class Alleles(np.ndarray, AlleleContainer):

    ''' A class for holding genotypes '''
    def __new__(cls, data, template=None, **kwargs):
        obj = np.asarray(data, **kwargs).view(cls)
        obj.template = template
        return obj

    def __array__finalize__(self, obj):
        if obj is None:
            return
        self.template = getattr(obj, 'template', None)

    def __lt__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __gt__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __le__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __ge__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    @property
    def missingcode(self):
        return 0 if np.issubdtype(self.dtype, np.integer) else ''

    @property
    def missing(self):
        " Returns a numpy array indicating which markers have missing data "
        return self == self.missingcode

    def nmark(self):
        '''
        Return the number of markers represented by the
        Alleles object
        '''
        return self.shape[0]

    def copy_span(self, template, copy_start, copy_stop):
        self[copy_start:copy_stop] = template[copy_start:copy_stop]

    def empty_like(self, blank=True):
        ''' Returns an empty Alleles object like this one '''
        z = np.zeros(self.nmark(), dtype=self.dtype)
        
        return Alleles(z, template=self.template)


class SparseAlleles(AlleleContainer):

    '''
    An object representing a set of haploid genotypes efficiently by 
    storing allele differences from a reference. Useful for manipulating
    genotypes from sequence data (e.g. VCF files)
    '''

    def __init__(self, data, template=None):
        self.template = template

        data = np.array(data)
        self.dtype = data.dtype
        self.size = data.shape[0]

        refcode = 0 if np.issubdtype(self.dtype, np.integer) else '0'
        self.non_refalleles = self._array2nonref(data, refcode)
        self.missingindices = self._array2missing(data, self.missingcode)

    def __lt__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __gt__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __le__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __ge__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    @staticmethod
    def _array2nonref(data, refcode):
        '''
        Returns a dict of the form index: value where the data is different than a reference
        '''
        return {i: x for i, x in enumerate(data)
                if x != refcode and x != ''}

    @staticmethod
    def _array2missing(data, missingcode):
        ''' Returns a list of indices where there are missingvalues '''
        return [i for i, x in enumerate(data) if x == missingcode]

    @property
    def missingcode(self):
        return 0 if np.issubdtype(self.dtype, np.integer) else ''

    @property
    def missing(self):
        " Returns a numpy array indicating which markers have missing data "
        base = np.zeros(self.size, dtype=np.bool_)
        base[self.missingindices] = 1
        return base

    def __eq__(self, other):
        if isinstance(other, SparseAlleles):
            return self.__speq__(other)
        elif isinstance(other, Alleles):
            return (self.todense() == other)
        elif np.issubdtype(type(other), self.dtype):
            if self.template is None:
                raise ValueError(
                    'Trying to compare values to sparse without reference')

            eq = np.array(self.template.reference, dtype=self.dtype) == other
            neq_altsites = [k for k, v in self.non_refalleles if k != other]
            eq_altsites = [k for k, v in self.non_refalleles if k == other]
            eq[neq_altsites] = False
            eq[eq_altsites] = True
            return eq
        else:
            raise ValueError(
                'Uncomparable types: {} and {}'.format(self.dtype,
                                                       type(other)))

    def __speq__(self, other):
        if self.size != other.size:
            raise ValueError('Trying to compare different-sized chromosomes')

        # SparseAlleles saves differences from a reference,
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

    def __ne__(self, other):
        return np.logical_not(self == other)

    def nmark(self):
        '''
        Return the number of markers (both reference and non-reference)
        represented by the SparseAlleles object
        '''
        return self.size

    def todense(self):
        '''
        Returns a non-sparse GenotypeChromosome equivalent
        to a SparseAlleles object.
        '''

        arr = np.zeros(self.size, dtype=np.uint8).astype(self.dtype)
        for loc, allele in self.non_refalleles.iteritems():
            arr[loc] = allele

        arr[self.missing] = self.missingcode

        return Alleles(arr, template=self.template)


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
        return 'Chromosome %s: %s markers, %s cM' % \
            (self.label if self.label is not None else 'object',
             len(self.frequencies),
             max(self.genetic_map) if self.genetic_map else 0)

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

    def linkageequilibrium_chromosome(self):
        """ Returns a randomly generated chromosome """
        if (self.frequencies < 0).any():
            raise ValueError('Not all frequencies are specified')
        r = np.random.random(self.nmark())
        r = np.array(r > self.frequencies, dtype=np.uint8) + 1
        return Alleles(r)

    def linkageequilibrium_chromosomes(self, nchrom):
        """ Returns a numpy array of many randomly generated chromosomes """
        chroms = np.random.random((nchrom, self.nmark()))
        chroms = np.uint8((chroms > self.frequencies) + 1)
        return [Alleles(r) for r in chroms]
