import numpy as np

class NotMeaningfulError(Exception):
    pass

class GenotypedChromosome(np.ndarray):
    ''' A class for holding genotypes '''
    def __new__(cls, data):
        obj = np.asarray(data).view(cls)
        return obj

    def __lt__(self, other):
        raise NotMeaningfulError('Value comparisions not meaningful for genotypes')
        
    def __gt__(self, other):
        raise NotMeaningfulError('Value comparisions not meaningful for genotypes')

    def __lte__(self, other):
        raise NotMeaningfulError('Value comparisions not meaningful for genotypes')
    
    def __gte__(self, other):
        raise NotMeaningfulError('Value comparisions not meaningful for genotypes')

    @property
    def missing(self):
        if np.issubdtype(self.dtype, np.integer):
            return self == 0 
        else:
            return self == '0'
