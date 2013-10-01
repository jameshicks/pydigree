#!/usr/bin/env python

import numpy as np
from blup import blup
from likelihood import makeV, restricted_loglikelihood

def is_genetic_effect(effect):
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
    def __init__(self, pedigrees, outcome=None, fixed_effects=None,
                 random_effects=None, covariance_matrices=None):
        self.maximized = False
        self.pedigrees = pedigrees
        self.outcome = outcome
        if not random_effects:
            self.random_effects = []
            self.covariance_matrices = []
        else: self.random_effects = random_effects
        self.fixed_effects = fixed_effects if fixed_effects else []
        self.variance_components = [None] * len(self.random_effects)
        self.obs = []
    def fit_model(self):
        self._makey()
        self._makeX()
        self._makeZs()
    def clear_model(self):
        """ Clears all parameters from the model """
        pass
    def observations(self):
        def has_all_fixefs(ind,effects):
            if not effects: return True
            return all(ind.phenotypes[effect] is not None for effect in effects)
        def has_outcome(ind):
            return ind.phenotypes[self.outcome] is not None
        obs = [x for x in self.pedigrees.individuals() \
               if has_all_fixefs(x,self.fixed_effects) and has_outcome(x)]
        if not self.obs: self.obs = obs
        return obs
    def nobs(self): return len(self.observations())
    def _makey(self):
        obs = self.observations()
        self.y = np.matrix([x.phenotypes[self.outcome] for x in obs]).transpose()
    def _makeX(self):
        """
        Builds the design matrix for the fixed effects in the mixed model.
        Includes a column of ones for the intercept.
        
        Arguements: None
        Returns: A numpy matrix
        """
        obs = self.observations()
        xmat = [[1] * len(obs)]
        for phen in self.fixed_effects: xmat.append([ob.phenotypes[phen] for ob in obs])
        import pdb; pdb.set_trace()
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
        for effect_name in self.random_effects:
            allinds = self.pedigrees.individuals()
            obs = frozenset(self.observations())
            obsidx = [i for i,x in enumerate(allinds) if x in obs]
            if is_genetic_effect(effect_name):
                incidence_matrix = np.matrix(np.eye(len(allinds)))[obsidx,:]
                Zlist.append(incidence_matrix)
            else: raise NotImplementedError('Arbitrary random effects not yet implemented')
        self.Zlist = Zlist
    def set_outcome(self,outcome):
        self.outcome = outcome
        self._makey()
    def add_fixed_effects(self,effect):
        self.fixed_effects.append(effect)
        self._makeX()
    def add_random_effect(self,effect,covariance_matrix):
        self.random_effects.append(effect)
        self.covariance_matrices.append(covariance_matrix)
        self._makeZs()
    def add_genetic_effect(self,type='additive'):
        self.add_random_effect(type,self.pedigrees.additive_relationship_matrix())
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
        if not all(self.variance_components):
            raise ValueError('Variance components not specified! Maximize the model or set them yourself.')
        fe,re = blup(self.y,self.X,self.Zlist,self.covmats,self.variance_components)
        self.fixef_blues = fe
        self.ranef_blups = re 
        return fe,re 
