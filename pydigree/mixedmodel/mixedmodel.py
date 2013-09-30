#!/usr/bin/env python

import numpy as np
from blup import blup
from likelihood import makeV, restricted_loglikelihood

def is genetic_effect(effect):
    return effect in set(['additive','dominance','mitochondrial'])

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
    def __init__(self, pedigrees=None, outcome=None, fixed_effects=None,
                 random_effects=None, covariance_matrices):
        self.maximized = False
        self.variance_components = [None] * len(random_effects)
        if not len(random_effects) == len(covariance_matrices):
            raise ValueError('Each random effect needs a covariance matrix, specify None for genetic effects')
        self.fit_model()
    def fit_model(self):
        self.makeX()
        self._makeZlist()
    def clear_model(self):
        """ Clears all parameters from the model """
        pass
    def observations(self):
        def has_all_fixefs(ind,effects):
            return all(ind[effect] is not None for effect in effects)
        obs = [x for x in self.pedigrees.individuals() if has_all_fixefs(x,self.fixed_effects)]
        if not self.obs: self.obs = obs
        return obs
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
        X = np.matrix(zip(*xmat))
        self.X = X
        return X
    def _makeZs(self):
        """
        Makes the incidence matrix for random effects

        Arguements: None
        Returns: A list of numpy matrixes
        """
        Zlist = []
        for effect_name in random_effects:
            allinds = self.pedigrees.individuals()
            obs = frozenset(self.observations())
            obsidx = [i for i,x in allinds if x in obs]
            if is_genetic_effect(effect_name):
                incidence_matrix = np.matrix(np.eye(len(allinds)))[obsidx,:]
                Zlist.append(incidence_matrix)
            else: raise NotImplementedError('Arbitrary random effects not yet implemented')
        self.Zlist = Zlist
    def add_fixed_effects(self,effect):
        pass
    def add_random_effect(self,effect,covariance_matrix):
        pass
    def set_variance_components(self,variance_components):
        if not all(x is not None for x in variance_components):
            raise ValueError('Not all variance components are specified')
        self.variance_components = variance_components
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
