"""
Functions for computing likelihoods of linear mixed models
"""
from __future__ import division

from math import log, pi

import numpy as np
from scipy.sparse import bsr_matrix, issparse
from scipy.linalg import pinv
from scipy.linalg import inv as scipy_inv
from scipy import matrix
# np.seterr(invalid='ignore')


l2pi = log(2 * pi)

def inv(M):
    if issparse(M):
        M = M.todense()
    return scipy_inv(M)

def logdet(M):
    """ Returns the (positive) log determinant of a matrix. """
    sign, logdet = np.linalg.slogdet(M.todense() if issparse(M) else M)
    return logdet


def makeP(X, Vinv):
    """ Makes the P matrix commonly found in mixed model estimation """
    return Vinv - Vinv * X * pinv(X.T * Vinv * X) * X.T * Vinv


def makeVinv(V):
    if issparse(V):
        V = V.todense()
    return bsr_matrix(inv(V))


def full_loglikelihood(y, V, X, beta, Vinv=None):
    """
    Returns the full loglikelihood of a mixed model

    Ref: SAS documentation for PROC MIXED
    """
    if Vinv is None:
        Vinv = makeVinv(V)
    n = X.shape[0]
    fixefresids = y - X * beta
    llik = -0.5 * (n * l2pi + logdet(V) + fixefresids.T * Vinv * fixefresids)
    return matrix.item(llik)


class MixedModelLikelihood(object):
    """
    A class describing the state of a mixed model likelihood function being 
    maximized
    """
    def __init__(self, mm, starts=None, info='fisher'):
        self.mm = mm
        
        if starts is not None:
            self.set_parameters(starts)
        if starts is None and all(vc is None for vc in mm.variance_components):
            raise ValueError('No variance components!')
        else:
            self.set_parameters(mm.variance_components)
        

        if info.lower() in {'fs', 'fisher', 'fisher scoring'}:
            self.set_info('fs')
        elif info.lower() in {'newton-raphson', 'newton', 'nr'}:
            self.set_info('nr')
        elif info.lower() in {'average information', 'aireml', 'ai'}:
            self.set_info('ai')
        elif info.lower() in {'em', 'emreml', 'expectation-maximization', 'grid'}:
            pass
        else:
            raise ValueError('Unknown maximization method')

        self.method = info
    
    def set_info(self, info):
        self.method = info

    def set_parameters(self, params):
        self.parameters = np.array(params)

        self.V = self.mm._makeV(params)
        self.Vinv = inv(self.V)
        self.P = makeP(self.mm.X, self.Vinv)

        # We need beta for non-reml computations
        X = self.mm.X
        self.beta = pinv(X.T * self.Vinv * X) * X.T * self.Vinv * self.mm.y
        self.resid = self.mm.y - X * self.beta

    def gradient(self):
        "The gradient of the likelihood function w/r/t each variance component"
        ranefs = self.mm.random_effects
        nabla = [self.gradient_element(rf.V_i) for rf in ranefs]
        return np.array(nabla)

    def info_matrix(self, kind=None):

        if not kind:
            kind = self.method
        
        if kind.lower() == 'fs':
            information_element = self.fisher_element
        elif kind.lower() == 'nr':
            information_element = self.observed_element
        elif kind.lower() == 'ai':
            information_element = self.ai_element
        elif kind.lower() == 'hessian':
            information_element = self.hessian_element
        else:
            raise ValueError('Unknown information matrix: {}'.format(kind))

        varmats = [x.V_i for x in self.mm.random_effects]
        nrf = len(varmats)

        mat = np.zeros((nrf, nrf))

        for i, V_i in enumerate(varmats):
            for j, V_j in enumerate(varmats):
                if j < i: 
                    continue

                element = information_element(V_i, V_j)

                mat[i, j] = element
                mat[j, i] = element

        return np.matrix(mat)

class ML(MixedModelLikelihood):

    def loglikelihood(self):
        n = self.mm.nobs()

        y, X, beta = self.mm.y, self.mm.X, self.beta
        resid = y - X * beta
        llik = -0.5 * (n*l2pi + logdet(self.V) + resid.T * self.Vinv * resid)
        return matrix.item(llik)

    def gradient_element(self, V_i):
        resid, Vinv = self.resid, self.Vinv

        return -0.5 * np.trace(Vinv * V_i) + 0.5 * resid.T * Vinv * V_i * Vinv * resid

    def hessian_element(self, V_i, V_j):
        resid, Vinv = self.resid, self.Vinv
        term1 = 0.5 * np.trace(Vinv * V_i * Vinv * V_j)
        term2 = resid.T * Vinv * V_i * Vinv * V_j * Vinv * resid

        return matrix.item(term1 - term2)

    def observed_element(self, V_i, V_j):
        return -1 * self.hessian_element(V_i, V_j)

    def fisher_element(self, V_i, V_j):
        Vinv = self.Vinv
        return 0.5 * np.trace(Vinv * V_i * Vinv * V_j)
    
    def ai_element(self, V_i, V_j):
        resid, Vinv = self.resid, self.Vinv
        return 0.5 * resid.T * Vinv * V_i * Vinv * V_j * Vinv * resid

    def expectation_maximization(self):
        "Performs a round of Expectation-Maximization ML"
        resid, Vinv = self.resid, self.Vinv

        n = self.mm.nobs()
        coefficients = np.array([
            matrix.item(resid.T * Vinv * rf.V_i * Vinv * resid - np.trace(Vinv * rf.V_i))
                        for rf in self.mm.random_effects])

        levelsizes = np.array([x.nlevels for x in self.mm.random_effects])

        delta = (self.parameters ** 2 / levelsizes) * coefficients
        return self.parameters + delta

class REML(MixedModelLikelihood):

    def loglikelihood(self):
        """
        Returns the restricted loglikelihood for mixed model variance component
        estimation.

        References:

        Harville. 'Maximum Likelihood Approaches to Variance Component
        Estimation and to Related Problems' Journal of the American Statistical
        Association. (1977) (72):258
        """
        y, V, X, P, Vinv = self.mm.y, self.V, self.mm.X, self.P, self.Vinv
        n = X.shape[0]

        rank = np.linalg.matrix_rank(X)
        llik_restricted = -0.5 * (logdet(V.todense())
                                  + logdet(X.transpose() * Vinv * X)
                                  + y.T * P * y
                                  + (n - rank) * l2pi)
        return matrix.item(llik_restricted)

    def gradient_element(self, V_i):
        y, P = self.mm.y, self.P
        term1 = -0.5 * np.trace(P * V_i)
        term2 = y.T * P * V_i * P * y
        return matrix.item(term1+term2)

    def hessian_element(self, V_i, V_j):
        y, P = self.mm.y, self.P
        common_term = P * V_i * P * V_j
        a = 0.5 * np.trace(common_term)
        b = y.T * common_term * P * y
        return matrix.item(a - b)

    def observed_element(self, V_i, V_j):        
        return -1.0 * self.hessian_element(V_i, V_j)
  
    def fisher_element(self, V_i, V_j):
        return .5 * np.trace(self.P * V_i * self.P * V_j)
    
    def ai_element(self, V_i, V_j):
        y, P = self.mm.y, self.P
        return .5 * y.T * P * V_i * P * V_j * P * y

    def expectation_maximization(self):
        "Performs a round of Expectation-Maximization REML"
        y, P = self.mm.y, self.P

        n = self.mm.nobs()
        coefficients = np.array([
            matrix.item(y.T * P * rf.V_i * P * y - np.trace(P * rf.V_i))
                        for rf in self.mm.random_effects])

        levelsizes = np.array([x.nlevels for x in self.mm.random_effects])

        delta = (self.parameters ** 2 / levelsizes) * coefficients
        return self.parameters + delta
