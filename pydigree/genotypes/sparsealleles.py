from collections import Sequence

import numpy as np

from pydigree.cydigree.datastructures import SparseArray
from pydigree.exceptions import NotMeaningfulError
from pydigree.genotypes import AlleleContainer, Alleles
from pydigree.exceptions import NotMeaningfulError
from pydigree.common import mode

class SparseAlleles(AlleleContainer):

    '''
    An object representing a set of haploid genotypes efficiently by 
    storing allele differences from a reference. Useful for manipulating
    genotypes from sequence data (e.g. VCF files)
    '''

    def __init__(self, data=None, refcode=None, missingcode='.', template=None):
        self.template = template

        if refcode is None:
            if data is None:
                raise IndexError
            else:
                refcode = mode(data)

        self.refcode = refcode

        if isinstance(refcode, str):
            self.dtype = np.dtype("S")
        elif isinstance(refcode, np.int):
            self.dtype = np.int
        else:
            raise IndexError

        if data is None:
            if template is None:
                raise ValueError('No template')
            self.container = SparseArray(len(data), refcode) 
            return 

        if type(data) is SparseArray:
            raise NotImplementedError
        
        else:    
            if not isinstance(data, np.ndarray):
                data = np.array(data)
            try:
                missingidx = np.where(data == missingcode)[0]
            except FutureWarning:
                assert 0
            # assert 0
            data[missingidx] = refcode
            self.container = SparseArray.from_dense(data, refcode)
            self.missingindices = set(missingidx)

        self.size = len(self.container)

    def __getitem__(self, key):
        return self.container[key]

    def __setitem__(self, key, value):
        self.container[key] = value

    @property
    def missingcode(self):
        return 0 if np.issubdtype(self.dtype, np.integer) else ''

    @property
    def missing(self):
        " Returns a numpy array indicating which markers have missing data "
        base = np.zeros(self.size, dtype=np.bool_)
        base[list(self.missingindices)] = 1
        return base

    def __eq__(self, other):
        if type(other) is SparseAlleles:
            return self.container == other.container
        else:
            return self.container == other

    def __ne__(self, other):
        if type(other) is SparseAlleles:
            return self.container != other.container
        else:
            return self.container != other

    def nmark(self):
        '''
        Return the number of markers (both reference and non-reference)
        represented by the SparseAlleles object
        '''
        return self.container.size

    def todense(self):
        raise NotImplementedError

    def empty_like(self):
        output = SparseAlleles(template=self.template,
                               missingcode=self.missingcode,
                               refcode=self.refcode)
        return output

    def copy_span(self, template, copy_start, copy_stop):
        raise NotImplementedError
