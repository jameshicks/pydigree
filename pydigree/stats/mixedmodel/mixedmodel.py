#!/usr/bin/env python
from __future__ import division

from itertools import izip
import copy

import numpy as np
from scipy.sparse import csc_matrix, issparse
from scipy.sparse import eye as sparseeye
from scipy.linalg import inv as scipy_inv
from scipy.linalg import pinv
from scipy.optimize import minimize

from pydigree.stats.mixedmodel.blup import blup, makeLHS
from pydigree.stats.mixedmodel.likelihood import makeP
from pydigree.stats.mixedmodel.likelihood import restricted_loglikelihood
from pydigree.stats.mixedmodel.likelihood import reml_gradient, reml_hessian
from pydigree.stats.mixedmodel.likelihood import full_loglikelihood

from pydigree.stats.mixedmodel.maximization import newtonlike_maximization
from pydigree.stats.mixedmodel.maximization import expectation_maximization_reml
from pydigree.stats.mixedmodel.maximization import MLEResult


def is_genetic_effect(effect):
    return effect in set(['additive', 'dominance', 'mitochondrial'])


def inv(M):
    if issparse(M):
        M = M.todense()
    return scipy_inv(M)


def _make_incidence_matrix(individuals, effect_name):
    if effect_name.lower() == 'residual':
        incidence_matrix = sparseeye(len(individuals))
    elif is_genetic_effect(effect_name):
        incidence_matrix = sparseeye(len(individuals))
    else:
        raise NotImplementedError('Arbitrary random effects not implemented')
    return incidence_matrix


class RandomEffect(object):
    __slots__ = ['label', 'variance_component',
                 'incidence_matrix', 'covariance_matrix']

    def __init__(self, individuals, label, variance=None,
                 incidence_matrix=None, covariance_matrix=None):
        nobs = len(individuals)

        self.label = label
        self.variance_component = variance
        if incidence_matrix is None:
            m = _make_incidence_matrix(individuals,
                                       self.label)
            self.incidence_matrix = m
        elif incidence_matrix == 'eye':
            self.incidence_matrix = sparseeye(nobs, nobs)
        else:
            self.incidence_matrix = incidence_matrix
        if covariance_matrix is None:

            self.covariance_matrix = sparseeye(nobs, nobs)
        else:
            self.covariance_matrix = covariance_matrix

    def __repr__(self):
        return 'Random Effect: {}'.format(self.label)

    # Convenience properties for linear algebra
    @property
    def sigma(self):
        "Convenience property for returning the variance of the component"
        return self.variance_component

    @property
    def Z(self):
        "Convenience property for returning the incidence matrix"
        return self.incidence_matrix

    @property
    def G(self):
        "Convenience property for returning the covariance_matrix"
        return self.covariance_matrix

    @property
    def V_i(self):
        return self.Z * self.G * self.Z.T


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
        self.mle = None
        self.pedigrees = pedigrees

        self.outcome = outcome
        self.fixed_effects = fixed_effects if fixed_effects else []
        self.obs = []

        if not random_effects:
            self.random_effects = []
        else:
            if not all(isinstance(x, RandomEffect) for x in random_effects):
                raise ValueError(
                    'Random effects must be of class RandomEffect')
            self.random_effects = random_effects

        # Extra variance component for residual variance. Works like
        # any other random effect.
        residual = RandomEffect(self.observations(), 'Residual')
        self.random_effects.append(residual)

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

    def _centery(self):
        self.y = self.y - self.y.mean()

    def fit_model(self):
        """
        Builds X, Z, Y, and R for the model

        Arguements: None
        Returns: Nothing
        """
        self.y = self._makey()
        self.X = self._makeX()
        self.Zlist = self._makeZs()

        if not all(x is not None for x in self.variance_components):
            self.set_variance_components(self.__starting_variance_components())

        self.V = self._makeV()
        self.beta = self._makebeta()

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
            try:
                return ind.phenotypes[self.outcome] is not None
            except KeyError:
                return False

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

    @property
    def P(self):
        "Projection matrix"
        return makeP(self.X, inv(self.V))

    def residual_variance(self):
        """ Returns the variance in y not accounted for by random effects """
        return self.random_effects[-1].sigma

    def coefficient_matrix(self):
        ginvlist = [(1 / rf.sigma) * inv(rf.G) for rf in self.random_effects]
        return makeLHS(self.X, self.Zlist, ginvlist, inv(self.R))

    def _fixef_vcv(self):
        "Variance-covariance matrix of fixed effect estimates"
        # Add 1 to include intercept
        nfixef = len(self.fixed_effects) + 1
        C = self.coefficient_matrix()
        vcv_beta = inv(C)[0:nfixef, 0:nfixef]
        return vcv_beta

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
        if (not vcs) and (not all(x is not None for x in self.variance_components)):
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
        if not isinstance(effect, RandomEffect):
            raise ValueError('Random effect must be of type RandomEffect')
        self.random_effects.insert(-1, effect)
        self.Zlist = self._makeZs()

    def add_genetic_effect(self, kind='additive'):
        inds = [x.full_label for x in self.observations()]
        peds = self.pedigrees

        if kind.lower() == 'additive':
            covmat = peds.additive_relationship_matrix(inds)
        elif kind.lower() == 'dominance':
            covmat = peds.dominance_relationship_matrix(inds)
        else:
            raise NotImplementedError(
                'Nonadditive/dominance genetic effects not implemented')
        effect = RandomEffect(
            self.observations(), kind, covariance_matrix=covmat)
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

    def maximize(self, method="Average Information",
                 starts=None, verbose=False):
        """
        Finds the optimal values for variance components of the model by
        restricted maximum likelihood estimation.
        """

        if (isinstance(self.mle, MLEResult) and
                self.maximized.method == method):
            return

        if starts is None:
            starts = self.__starting_variance_components()

        if method == 'scipy':
            mle = self._maximize_scipy(verbose=verbose)
        elif method.lower() in {'em', 'emreml', 'expectation-maximization'}:
            mle = expectation_maximization_reml(self, starts, verbose=verbose)
        else:
            mle = newtonlike_maximization(self, starts, method,
                                          verbose=verbose)

        self.mle = mle
        self.set_variance_components(mle.parameters)
        self.fit_model()

        # Get the full loglikelihood at the REML maximimum so we
        # can use it later
        self.mle.full_loglikelihood = full_loglikelihood(self.y, self.V,
                                                         self.X, self.beta)

    @property
    def maximized(self):
        return isinstance(self.mle, MLEResult)

    def _maximize_scipy(self, method='L-BFGS-B', verbose=True):
        """
        Finds the optimal values for variance components of the model by
        restricted maximum likelihood estimation using scipy minimization.

        Only for debugging use at the moment.
        """

        starts = self.__starting_variance_components()

        def cb(x):
            if verbose:
                print 'Iteration VC estimates: %s' % \
                    ', '.join(str(y) for y in x.tolist())

        r = minimize(self._reml_optimization_target, x0=starts,
                     method='CG',
                     jac=self._reml_gradient,
                     callback=cb)

        if not r.success:
            raise np.linalg.LinAlgError('Model fit did not converge!')

        if verbose:
            print r

        mle = MLEResult(r.x, r.fun, 'scipy', r.jac, r.hess)
        return mle

    def loglikelihood(self, restricted=False, vmat=None):
        """
        Returns the loglikelihood of the model with the current model parameters
        """
        if self.mle is not None:
            if restricted:
                return self.mle.restricted_loglikelihood
            else:
                return self.mle.full_loglikelihood

        if self.V is None:
            self.V = self._makeV()
        if vmat is None:
            V = self.V
        else:
            V = vmat
        if not restricted:
            return full_loglikelihood(self.y, V, self.X, self.beta)
        else:
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
        if not self.maximized:
            raise ValueError('Model not maximized!')

        self._fit_results()

        print
        print 'Linear mixed model fit by {}'.format(self.maximized)
        print
        print 'Fixed effects:'
        fixefnames = ['(Intercept)'] + self.fixed_effects
        betas = self.beta.T.tolist()[0]
        print '\t'.join(['Name', 'Estimate'])
        for name, beta in zip(fixefnames, betas):
            print '\t'.join(q for q in [name, '{:5.3f}'.format(beta)])
        print
        print 'Variance components:'
        print '\t'.join(['Component', 'Variance', '% Variance'])
        totalvar = sum(self.variance_components)
        for effect, vc in zip(self.random_effects, self.variance_components):
            print '\t'.join(v for v in [effect.label,
                                        '{:5.3f}'.format(vc),
                                        '{:5.3f}'.format(100 * vc / totalvar)])
        print
        print 'Loglikelihood: {:10.2f}'.format(self.loglikelihood())
        print 'Observations: {}'.format(self.nobs())
        print

    def _reml_optimization_target(self, vcs):
        """ Optimization target for maximization. """
        Q = self._makeV(vcs=vcs.tolist())
        return -1.0 * self.loglikelihood(restricted=True, vmat=Q)

    def _reml_gradient(self, vcs):
        Q = self._makeV(vcs=vcs.tolist())
        return -1.0 * reml_gradient(self.y, self.X, Q, self.random_effects)

    def _reml_hessian(self, vcs):
        Q = self._makeV(vcs.tolist())
        return -1.0 * reml_hessian(self.y, self.X, Q, self.random_effects)

    def __starting_variance_components(self):
        """
        Starting variance components in optimization.
        Chooses all variance components (including residual) to be equal.
        """
        v = np.var(self.y)
        n = float(len(self.random_effects))
        return [v / (n + 1) for r in self.random_effects]
