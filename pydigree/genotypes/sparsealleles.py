from pydigree.genotypes import AlleleContainer, Alleles
import numpy as np
from itertools import izip
from pydigree.exceptions import NotMeaningfulError
from pydigree.cyfuncs import fastfirstitem

class SparseAlleles(AlleleContainer):

    '''
    An object representing a set of haploid genotypes efficiently by 
    storing allele differences from a reference. Useful for manipulating
    genotypes from sequence data (e.g. VCF files)
    '''

    def __init__(self, data, refcode=None, template=None):
        self.template = template

        data = np.array(data)
        self.dtype = data.dtype
        self.size = data.shape[0]
        if refcode is not None:
            self.refcode = refcode
        else:
            self.refcode = 0 if np.issubdtype(self.dtype, np.integer) else '0'
        self.non_refalleles = self._array2nonref(data,
                                                 self.refcode,
                                                 self.missingcode)
        self.missingindices = self._array2missing(data,
                                                  self.missingcode)

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
    def _array2nonref(data, refcode, missingcode):
        '''
        Returns a dict of the form index: value where the data is 
        different than a reference value
        '''
        idxes = np.where(np.logical_and(data != refcode,
                                        data != missingcode))[0]
        nonref_values = data[idxes]
        return dict(izip(idxes, nonref_values))

    @staticmethod
    def _array2missing(data, missingcode):
        ''' Returns a list of indices where there are missingvalues '''
        return list(np.where(data == missingcode)[0])

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

    def empty_like(self):
        if not np.issubdtype(self.dtype, np.int):
            raise ValueError 
        raw = np.zeros(self.nmark(), dtype=self.dtype) + self.refcode
        return SparseAlleles(raw, refcode=self.refcode, template=self.template)
