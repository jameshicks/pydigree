''' A mixin class for building mixed models '''

from itertools import izip

import numpy as np

from scipy.linalg import inv
from scipy.sparse import issparse
from scipy.sparse.linalg import inv as spinv

class MixedModelMixin(Object):
    def __getmatrix(self, mtype):
        if mtype == 'additive': return self.additive_relationship_matrix()
        elif mtype == 'dominance': return self.dominance_relationship_matrix()
        elif mtype == 'mitochrondrial': return self.mitochondrial_relationship_matrix()
        else:
            raise ValueError('Invalid pedigree covariance matrix type: {}'.format(mtype))

    def __sort_inds_in_ped(self, indlist):
        ''' Takes a list of individuals, filters out the ones in the current pedigree and sorts those '''
        return sorted(x for x in indlist if x.pedigree.label == self.label,
                      key=lambda x: (x.pedigree.label, x.id))

    def vcov_matrix(self, incidence_matrices, covariance_matrices, variance_components):
        """ Returns the variance-covariance matrix (V) for individuals in a mixed model """
        return sum(sigma * Z * __getmatrix(G) * Z.T
                   for Z, G, sigma in izip(incidence_matrices, covariance_matrices, variance_components))

    def inv_vcov_matrix(self, incidence_matrices, covariance_matrices, variance_components):
        V = vcov_matrix(incidence_matrices, covariance_matrices, variance_components)
        return spinv(V) if issparse(V) else inv(V)
