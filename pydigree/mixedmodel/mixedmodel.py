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
    def __init__(self,pedigrees=None,outcome=None,fixed_effects=None,random_effects=None):
        self.maximized = False
        self.obs = None
    def clear_model(self):
        """ Clears all parameters from the model """
        self.maximized = False
        self.obs = []
        self.X = None
        self.Zlist = []
    def observations(self):
        def has_all_fixefs(ind,effects):
            return all(ind[effect] is not None for effect in effects)
        obs = [x for x in self.pedigrees.individual_chain() if has_all_fixefs(x,self.fixed_effects)]
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
            allinds = self.pedigrees.individual_chain()
            obs = frozenset(self.observations())
            obsidx = [i for i,x in allinds if x in obs]
            if is_genetic_effect(effect_name):
                incidence_matrix = np.matrix(np.eye(len(allinds)))[obsidx,:]
                Zlist.append(incidence_matrix)
            else: raise NotImplementedError('Arbitrary random effects not yet implemented')
        self.Zlist = Zlist
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
