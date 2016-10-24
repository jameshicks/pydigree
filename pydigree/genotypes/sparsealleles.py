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

    In the interest of conserving memory for sequencing data, all alleles must
    be represented by a signed 8-bit integer (i.e. between -128 and 127). 
    Negative values are interpreted as missing.
    '''

    def __init__(self, data=None, refcode=0, missingcode=-1, size=None, template=None):
        self.template = template

        self.refcode = refcode

        if data is None:
            if template is None and size is None:
                raise ValueError('No template or size')
            elif template is not None and size is None:
                size = self.template.nmark()
            self.container = SparseArray(size, refcode) 
            return 

        if type(data) is SparseArray:
            raise NotImplementedError
        
        else:    
            if not isinstance(data, np.ndarray):
                data = np.array(data)
            
            self.container = SparseArray.from_dense(data, refcode)


        self.size = len(self.container)

    def __getitem__(self, key):
        return self.container[key]

    def __setitem__(self, key, value):
        self.container[key] = value

    def keys(self):
        return self.container.keys()
    
    def values(self):
        return self.container.values()

    @property
    def missingcode(self):
        return -1

    @property
    def missing(self):
        " Returns a numpy array indicating which markers have missing data "
        missingindices = [i for i,v in self.container.items() if v == self.missingcode]
        base = np.zeros(self.size, dtype=np.bool_)
        base[missingindices] = 1
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
        dense = Alleles(self.container.tolist(), template=self.template)
        return dense

    def empty_like(self):
        output = SparseAlleles(template=self.template,
                               missingcode=self.missingcode,
                               refcode=self.refcode, size=self.nmark())
        return output

    def copy_span(self, template, copy_start, copy_stop):
        if isinstance(template, SparseAlleles):
            self.container[copy_start:copy_stop] = template.container[copy_start:copy_stop]
        else:
            self.container = template[copy_start:copy_stop]

    def copy(self):
        outp = self.empty_like()
        outp.container = self.container.copy()
        return outp

    @staticmethod
    def empty(reference=None, template=None, missingcode=''):
        out = SparseAlleles(size, template=template, missingcode=missingcode)

        return out 