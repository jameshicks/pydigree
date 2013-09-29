#!/usr/bin/env python

import numpy as np

class MixedModel(object):
    def __init__(self):
        pass
    def _makeX(self):
        """
        Builds the design matrix for the fixed effects in the mixed model.
        Includes a column of ones for the intercept.
        
        Arguements: None
        Returns: A numpy matrix
        """
        if len(set(len(x) for x in self.fixed_effects)) > 1:
            raise ValueError('Not all fixed effects have same number of observations')
        nobs = len(self.fixed_effects[0])
        mat = [[1] * nobs] + self.fixed_effects
        return np.matrix(mat).transpose()
    def maximize(self):
        pass
    def likelihood(self,metric='REML'):
        pass
    def breeding_values(self):
        pass
    
