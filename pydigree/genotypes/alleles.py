import numpy as np

from pydigree.exceptions import NotMeaningfulError
from .genoabc import AlleleContainer

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
