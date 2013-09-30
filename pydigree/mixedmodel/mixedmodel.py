#!/usr/bin/env python

import numpy as np
from blup import blup
from likelihood import makeV, restricted_loglikelihood

class MixedModel(object):
    """
    Fits linear models in the form of y = X * b + sum(Z_i * u_i) + e, where:
      y is the vector of outcomes
      X is a design matrix of fixed effects
      b is the vector of coefficients correstponding to those fixed effects
      Z_i is an incidence matrix corresponding to random effect i
      u_i is a vector of values corresponding to random effect i
      e is a vector of errors
    """
    def __init__(self,pedigrees=None,outcome=None,fixed_effects=None,random_effects=None):
        self.maximized = False
        self.obs = None
    def observations(self):
        def has_all_fixefs(ind,effects):
            return all(ind[effect] is not None for effect in effects)
        return [x for x in self.pedigrees.individual_chain() if has_all_fixefs(x,self.fixed_effects)]
    def nobs(self): return len(self.observations())
    def _makeX(self):
        """
        Builds the design matrix for the fixed effects in the mixed model.
        Includes a column of ones for the intercept.
        
        Arguements: None
        Returns: A numpy matrix
        """
        obs = self.observations()
        xmat = [[1] * len(obs)]
        for ob in obs: xmat.append([x.phenotypes[phen] for phen in self.fixed_effects])
        return np.matrix(zip(*xmat))
    def _makeZs(self):
        pass 
    def add_fixed_effects(self):
        pass
    def maximize(self,method='anneal'):
        if self.maximized == method: return 
        pass
    def likelihood(self,metric='REML'):
        pass
    def blup(self):
        if not self.maximized:
            raise ValueError('Model not maximized!')
        fe,re = blup(self.y,self.X,self.Zlist,self.covmats,self.variance_components)
        self.fixef_blues = fe
        self.ranef_blups = re 
        return fe,re 
