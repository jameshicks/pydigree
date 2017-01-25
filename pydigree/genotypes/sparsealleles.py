import numpy as np

from pydigree.cydigree.sparsearray import SparseArray
from pydigree.genotypes import AlleleContainer, Alleles
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

    def __init__(self, data=None, refcode=0, size=None, template=None):
        self.template = template

        if refcode is None:
            refcode = 0

        if data is None:
            if template is None and size is None:
                raise ValueError('No template or size')
            elif template is not None and size is None:
                size = self.template.nmark()
            self.container = SparseArray(size, refcode) 
            return 

        elif type(data) is SparseArray:
            self.container = data.copy()
        
        else:    
            if not isinstance(data, np.ndarray):
                data = np.array(data)
            ref = refcode if refcode is not None else mode(data)
            self.container = SparseArray.from_dense(data, ref)


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
    def refcode(self):
        """ 
        Returns the sparse value in the container 
        
        :rtype: int8_t
        """
        return self.container.ref

    @property
    def missingcode(self):
        "Returns the code used for missing values"
        return -1

    @property
    def dtype(self):
        return int

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

        :returns: markercount
        :rtype: int
        '''
        return self.container.size

    def todense(self):
        """
        Converts to a dense representation of the same genotypes (Alleles).

        :returns: dense version
        :rtype: Alleles
        """
        dense = Alleles(self.container.tolist(), template=self.template)
        return dense

    def empty_like(self):
        """
        Creates a blank SparseAlleles with same parameters

        :returns: empty SparseAlleles
        """
        output = SparseAlleles(template=self.template,
                               refcode=self.refcode, size=self.nmark())
        return output

    def copy_span(self, template, copy_start, copy_stop):
        """
        Copies one segment of a chromosome over to the other

        :param template: the data to be copied from
        :param copy_start: where to start copying (inclusive)
        :param copy_stop: where to stop copying (exclusive)
        :type template: AlleleContainer
        :type copy_start: int
        :type copy_stop: int
        :rtype void:
        """
        if isinstance(template, SparseAlleles):
            self.container[copy_start:copy_stop] = template.container[copy_start:copy_stop]
        else:
            self.container[copy_start:copy_stop] = template[copy_start:copy_stop]

    def copy(self):
        """
        Creates a copy of the current data

        :returns: cloned allele set
        :rtype: SparseAlleles
        """
        outp = self.empty_like()
        outp.container = self.container.copy()
        return outp

    @staticmethod
    def empty(template):
        """
        Creates an empty SparseAlleles (everybody is wild-type)

        :param template: The chromosome info associated with this set of alleles
        :type template: ChromosomeTemplate

        :returns: Empty container
        :rtype: SparseAlleles
        """
        out = SparseAlleles(template.nmark(), template=template)

        return out 