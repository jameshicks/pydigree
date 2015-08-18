#!/usr/bin/env python

from itertools import izip
import copy

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import eye as sparseeye
from scipy.linalg import inv, pinv
from scipy.optimize import fmin_l_bfgs_b

from pydigree.mixedmodel.blup import blup
from pydigree.mixedmodel.likelihood import restricted_loglikelihood
from pydigree.mixedmodel.likelihood import full_loglikelihood


def is_genetic_effect(effect):
    return effect in set(['additive', 'dominance', 'mitochondrial'])


def _make_incidence_matrix(individuals, effect_name):
    if is_genetic_effect(effect_name):
        incidence_matrix = np.matrix(np.eye(len(individuals)))
    else:
        raise NotImplementedError('Arbitrary random effects not implemented')
    return incidence_matrix


class RandomEffect(object):
    __slots__ = ['variance_component', 'incidence_matrix', 'covariance_matrix']

    def __init__(self, individuals, label, variance,
                 incidence, covariance_matrix):
        self.label = label
        self.variance = variance_component
        if not incidence:
            m = _make_incidence_matrix(individuals,
                                       self.label)
            self.incidence_matrix = m
        else:
            self.incidence_matrix = incidence
        if not covariance_matrix:
            nobs = len(individuals)
            self.covariance_matrix = sparseeye(nobs, nobs)
        else:
            self.covariance_matrix = covariance_matrix

    # Convenience properties for linear algebra
    @property
    def sigma(self):
        "Convenience property for returning the variance of the component"
        return self.variance_component

    @property
    def Z(self):
        "Convenience property for returning the incidence matrix"
        return self.incidence

    @property
    def G(self):
        "Convenience property for returning the covariance_matrix"
        return self.covariance_matrix


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
        else:
            if not all(isinstance(RandomEffect, x) for x in random_effects):
                raise ValueError(
                    'Random effects must be of class RandomEffect')
            self.random_effects = random_effects
        
        # Extra variance component for residual variance. Works like
        # any other random effect.
        residual = RandomEffect(self.observations(), 'Residual')
        self.random_effects.append(residual)

        self.fixed_effects = fixed_effects if fixed_effects else []
        self.obs = []
        self.V = None

    def copy(self):
        "Return a copy of the model"
        # We want to avoid copying pedigree and individual data, so
        # we'll set the pedigrees attribute to None for a sec, and then
        # change it back
        peds = self.pedigrees
        self.pedigrees = None

        newmm = copy.deepcopy(self)

        newmm.pedigrees = peds
        self.pedigrees = peds

        return newmm

    def fit_model(self):
        """
        Builds X, Z, Y, and R for the model

        Arguements: None
        Returns: Nothing
        """
        self.y = self._makey()
        self.X = self._makeX()
        self.Zlist = self._makeZs()

    def _fit_results(self):
        self.V = self._makeV()
        self.beta = self._makebeta()

    def clear_model(self):
        """ Clears all parameters from the model """
        self.random_effects = []
        self.fixed_effects = []
        self.Zlist = []
        self.X = None
        self.y = None
        self.beta = None
        self.V = None
        self.obs = None

    def observations(self):
        """
        Returns a list of the fully observed individuals in the model
        Fully observed individuals have observations for each fixed effect and
        and observation for the outcome variable.

        Arguements: None
        Returns: The list of fully observed individuals
        """

        def has_all_fixefs(ind, effects):
            if not effects:
                return True
            return all(ind.phenotypes[effect] is not None
                       for effect in effects)

        def has_outcome(ind):
            return ind.phenotypes[self.outcome] is not None

        obs = [x for x in self.pedigrees.individuals
               if has_all_fixefs(x, self.fixed_effects) and has_outcome(x)]
        if not self.obs:
            self.obs = obs
        return obs

    def nobs(self):
        """
        Returns the number of fully observed individuals in the pedigree
        Arguements: None
        Returns: An integer
        """
        return len(self.observations())

    @property
    def variance_components(self):
        return [x.sigma for x in self.random_effects]

    @property
    def covariance_matrices(self):
        return [x.covariance_matrix for x in self.random_effects]

    @property
    def R(self):
        "Covariance matrix of residual_variance"
        return self.random_effects[-1].covariance_matrix

    def residual_variance(self):
        """ Returns the variance in y not accounted for by random effects """
        return self.random_effects[-1].sigma

    def _makey(self):
        """ Prepares the vector of outcome variables for model estimation """
        obs = self.observations()
        return np.matrix([x.phenotypes[self.outcome] for x in obs]).transpose()

    def _makeX(self):
        """
        Builds the design matrix for the fixed effects in the mixed model.
        Includes a column of ones for the intercept.

        Arguements: None
        Returns: A numpy matrix
        """
        obs = self.observations()
        xmat = [[1] * len(obs)]
        for phen in self.fixed_effects:
            xmat.append([ob.phenotypes[phen] for ob in obs])
        X = np.matrix(zip(*xmat))
        return X

    def _makeZs(self):
        """
        Makes the incidence matrix for random effects

        Arguements: None
        Returns: A list of numpy matrixes
        """
        Zlist = [ranef.Z for ranef in self.random_effects]

        return [csc_matrix(Z) for Z in Zlist]



    def _makeV(self, vcs=None):
        if (not vcs) and (not self.variance_components):
            raise ValueError('Variance components not set')
        if not vcs:
            variance_components = self.variance_components
        else:
            variance_components = vcs
        
        V = sum(sigma * Z * A * Z.T for sigma, Z, A in
                izip(variance_components,
                     self.Zlist,
                     self.covariance_matrices))
        
        return V

    def _makebeta(self):
        """
        Calculates BLUEs for the fixed effects portion of the model

        Reference:
        McCulloch & Seale. Generalized, Linear, and Mixed Models. (2001)
        Equation 6.24
        """
        vinv = inv(self.V.todense())
        return pinv(self.X.T * vinv * self.X) * self.X.T * vinv * self.y

    def set_outcome(self, outcome):
        """ Sets the outcome for the mixed model """
        self.outcome = outcome
        self.y = self._makey()

    def add_fixed_effects(self, effect):
        """ Adds a fixed effect to the model """
        self.fixed_effects.append(effect)
        self.X = self._makeX()

    def add_random_effect(self, effect):
        """ Adds a random effect to the model """
        if not isinstance(RandomEffect, effect):
            raise ValueError('Random effect must be of type RandomEffect')
        self.random_effects.insert(-1, effect)
        self.Zlist = self._makeZs()

    def add_genetic_effect(self, type='additive'):
        if type.lower() != 'additive':
            raise NotImplementedError(
                'Nonadditive genetic effects not implemented')
        inds = [x.label for x in self.observations()]
        peds = self.pedigrees
        covmat = self.add_random_effect(type,
                               peds.additive_relationship_matrix(inds))
        effect = RandomEffect(self.observations, type, covariance_matrix=covmat)
        self.add_random_effect(effect)

    def set_variance_components(self, variance_components):
        """
        Manually set variance components for each random effect in the model.
        Useful if you know a priori, say a heritability, and just want to
        predict breeding values for the trait.
        """
        if not all(x is not None for x in variance_components):
            raise ValueError('Not all variance components are specified')
        for sigma, ranef in izip(variance_components, self.random_effects):
            ranef.variance_component = sigma

    def maximize(self, method='L-BFGS-B', verbose=False):
        """
        Finds the optimal values for variance components of the model by
        restricted maximum likelihood estimation.
        """
        if self.maximized == method:
            return
        starts = self.__starting_variance_components()

        def cb(x):
            if verbose:
                print 'Iteration VC estimates: %s' % \
                    ', '.join(str(y) for y in x.tolist())

        b = [(0, np.var(self.y))] * len(self.random_effects)
        r = fmin_l_bfgs_b(self.__reml_optimization_target,
                          starts, bounds=b, approx_grad=1, callback=cb)
        if verbose:
            print r
        self.set_variance_components(r[0].tolist())

    def likelihood(self, estimator='full', vmat=None):
        """
        Returns the likelihood of the model with the current model parameters
        """
        estimator = estimator.lower()
        if vmat is None:
            V = self.V
        else:
            V = vmat
        if estimator == 'full':
            return full_loglikelihood(self.y, V, self.X)
        elif estimator == 'restricted':
            return restricted_loglikelihood(self.y, V, self.X)

    def blup(self):
        if not all(self.variance_components):
            raise ValueError('Varcomps not specified! Maximize or set them')
        b = blup(self.y, self.X, self.Zlist,
                 self.covariance_matrices,
                 self.variance_components).transpose().tolist()[0]
        return b

    def summary(self):
        """ Prints a summary of the current model """
        self._fit_results()
        print 'Fixed effects:'
        fixefnames = ['Intercept'] + self.fixed_effects
        betas = self.beta.T.tolist()[0]
        print '\t'.join(['Name', 'Estimate'])
        for name, beta in zip(fixefnames, betas):
            print '\t'.join(str(q) for q in [name, beta])
        print
        print 'Variance components:'
        print '\t'.join(['Component', 'Variance', '% Variance'])
        for name, vc in zip(self.random_effects, self.variance_components):
            print '\t'.join(str(val) for val in [name,
                                                 vc,
                                                 100 * vc / np.var(self.y)])
        print
        print 'Loglikelihood: %s' % self.likelihood()

    def __reml_optimization_target(self, vcs):
        """ Optimization target for maximization. """
        Q = self._makeV(vcs=vcs.tolist())
        return -1.0 * self.likelihood(estimator='restricted', vmat=Q)

    def __starting_variance_components(self):
        """
        Starting variance components in optimization.
        Chooses all variance components (including residual) to be equal.
        """
        v = np.var(self.y)
        n = float(len(self.random_effects))
        return [v / (n + 1) for r in self.random_effects]
