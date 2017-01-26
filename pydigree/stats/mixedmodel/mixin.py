''' A mixin class for building mixed models '''

import numpy as np

from scipy.linalg import inv
from scipy.sparse import issparse
from scipy.sparse.linalg import inv as spinv


class MixedModelMixin(object):
    "Mixin class for populations to implement mixed model functions"

    def __getmatrix(self, mtype):
        "Get a genetic relationship matrix"
        if mtype == 'additive':
            return self.additive_relationship_matrix()
        elif mtype == 'dominance':
            return self.dominance_relationship_matrix()
        elif mtype == 'mitochrondrial':
            return self.mitochondrial_relationship_matrix()
        else:
            raise ValueError(
                'Invalid pedigree covariance matrix type: {}'.format(mtype))

    def __sort_inds_in_ped(self, indlist):
        """ 
        Takes a list of individuals, filters out the ones in the current 
        pedigree and sorts those.

        :param indlist: individuals to sort 
        """
        return sorted((x for x in indlist if x.pedigree.label == self.label),
                      key=lambda x: (x.pedigree.label, x.label))

    def incidence_matrix(self, variable=None, inds=None, onlylevels=None):
        """
        Generates an incidence matrix for random effect in a mixed model based 
        on variable. If no variable is given, the individual is used (i.e. an 
        identity matrix).
        
        :param variable: phenotype to form the matrix for
        :param inds: only use these individuals
        :param onlylevels: only use these levels of the random effect
        """
        if variable is None:
            getvar = lambda ind: ind.label
        else:
            getvar = lambda ind: ind.phenotypes[variable]

        levels = sorted({getvar(ind) for ind in self.individuals})

        if onlylevels is not None:
            onlylevels = set(onlylevels)
        else:
            # If we're not reducing the number of levels, we'll set onlylevels 
            # to levels so the intersection of levels and onlylevels reduces 
            # to just levels
            onlylevels = levels

        levels = sorted(levels & onlylevels)

        # Test for cases incompatible for the mixed model
        if not levels:
            raise ValueError('No valid levels for variable!')
        elif len(levels) == 1:
            raise ValueError('Variable only has one level!')

        if not inds:
            inds = self.individuals

        Z = []
        for ind in inds:
            Z.append([getvar(ind) == level for level in levels])

        # Make sure everyone has a 1 in their row for the incidence matrix
        if not all(sum(row) for row in Z):
            raise ValueError(
                'Individuals are missing values in random effects!')

        return np.matrix(Z, dtype=np.int8)
