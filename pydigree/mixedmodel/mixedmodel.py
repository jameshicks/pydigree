#!/usr/bin/env python

from itertools import izip

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import eye as sparseeye
from scipy.optimize import fmin_l_bfgs_b

from blup import blup
from likelihood import restricted_loglikelihood

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
        self.V = None
    def fit_model(self):
        self._makeR()
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
    def residual_variance(self):
        return np.var(self.y) - sum(self.variance_components) 
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
        self.Zlist = [csc_matrix(Z) for Z in Zlist]
    def _makeR(self):
        self.R = sparseeye(self.nobs())
    def _makeV(self,vcs=None):
        if (not vcs) and (not self.variance_components):
            raise ValueError('Variance components not set')
        if not vcs: variance_components = self.variance_components
        else: variance_components = vcs
        V = sum(sigma * Z * A * Z.T for sigma,Z,A in \
                izip(variance_components, self.Zlist, self.covariance_matrices))
        V = V + (np.var(self.y) - sum(variance_components)) * self.R
        if vcs is not None: return V
        else: self.V = V
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
    def maximize(self,method='L-BFGS-B'):
        if self.maximized == method: return
        starts = self.__starting_variance_components()
        b = [(0,np.var(self.y))] * len(self.random_effects)
        def cb(x): print x
        r = fmin_l_bfgs_b(self.__reml_optimization_target,starts,bounds=b,approx_grad=1,callback=cb)
        import pdb; pdb.set_trace()
        self.variance_components = r[0].tolist()
    def likelihood(self,vmat=None):
        if vmat is None: V = self.V
        else: V = vmat
        return restricted_loglikelihood(self.y, V,self.X)
    def blup(self):
        if not all(self.variance_components):
            raise ValueError('Variance components not specified! Maximize the model or set them yourself.')
        fe,re = blup(self.y,self.X,self.Zlist,self.covariance_matrices,self.variance_components)
        self.fixef_blues = fe
        self.ranef_blups = re 
        return fe,re 
    def __reml_optimization_target(self,vcs):
        Q = self._makeV(vcs=vcs.tolist())
        return self.likelihood(vmat=Q)
    def __starting_variance_components(self):
        v = np.var(self.y)
        n = float(len(self.random_effects))
        return [v/(n+1) for r in self.random_effects]
