import numpy as np

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
        """
        Returns a numpy array indicating which markers have missing data
        
        :returns: missingness array
        :rtype: np.array
        """
        return np.array(self == self.missingcode)

    def nmark(self):
        '''
        Return the number of markers represented by the Alleles object

        :returns: number of markers
        :rtype: int
        '''
        return self.shape[0]

    def copy_span(self, template, copy_start, copy_stop):
        """
        Copies a span of another AlleleContainer to this one

        :param template: Container to copy from
        :type template: AlleleContainer
        :param copy_start: start point for copy (inclusive)
        :type copy_start: int
        :param copy_stop: end_point for copy (exclusive)
        :type copy_stop: int

        :rtype: void
        """
        self[copy_start:copy_stop] = template[copy_start:copy_stop]

    def empty_like(self):
        ''' Returns an empty Alleles object like this one '''
        z = np.zeros(self.nmark(), dtype=self.dtype)

        return Alleles(z, template=self.template)
