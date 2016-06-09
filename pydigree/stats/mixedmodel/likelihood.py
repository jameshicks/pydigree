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
        elif info.lower() in {'em', 'emreml', 'expectation-maximization'}:
            pass
        else:
            raise ValueError('Unknown maximization method')

        self.method = info
        
    def set_parameters(self, params):
        self.parameters = np.array(params)

        self.V = self.mm._makeV(params)
        self.Vinv = inv(self.V)
        self.P = makeP(self.mm.X, self.Vinv)

        # We need beta for non-reml computations
        X = self.mm.X
        self.beta = pinv(X.T * self.Vinv * X) * X.T * self.Vinv * self.mm.y

class ML(MixedModelLikelihood):
    def set_info(self, info):
        d = {'fs': self.ml_fisher_information_matrix,
             'nr': self.ml_observed_information_matrix,
             'ai': self.ml_average_information_matrix}
    
        self.hessian = d[info]

    def gradient(self):
        "The gradient of the ML function w/r/t each variance component"

        def dML_dsigma(ranef):
            "The ML derivative with regard to a variance component"

            Z, G = ranef.Z, ranef.G
            y, X, beta, Vinv = self.mm.y, self.mm.X, self.beta, self.Vinv
            V_i = Z * G *Z.T

            resid = y - X * beta
            term1 = -0.5 * np.trace(Vinv * V_i)
            term2 = 0.5 * resid.T * Vinv * V_i * Vinv * resid
            return term1 + term2 

        ranefs = self.mm.random_effects
        nabla = [dML_dsigma(rf) for rf in ranefs]
        return np.array(nabla)

    def loglikelihood(self):
        n = self.mm.nobs()

        y, X, beta = self.mm.y, self.mm.X, self.beta
        resid = y - X * beta
        return -0.5 * (n*l2pi + logdet(self.V) + resid.T * self.Vinv * resid)

    def ml_hessian(self):
        ranefs = self.mm.random_effects
        n_ranefs = len(ranefs)
        mat = np.zeros((n_ranefs, n_ranefs))

        resid = self.mm.y - self.mm.X * self.beta
        
        def ml_hessian_element(resid, Vinv, V_i, V_j):
            common_term = Vinv * V_i * Vinv * V_j
            a = 0.5 * np.trace(common_term)
            b = resid.T * common_term * resid
            return matrix.item(a - b)

        for i, ranef_a in enumerate(ranefs):
            dV_dsigma_a=ranef_a.V_i
            for j, ranef_b in enumerate(ranefs):

                if j < i:
                    # Already set when we did the other side of the matrix
                    continue

                dV_dsigma_b=ranef_b.V_i
                element=ml_hessian_element(resid, self.Vinv, dV_dsigma_a, dV_dsigma_b)

                mat[i, j]=element
                mat[j, i]=element

        return np.array(mat)

    def ml_observed_information_matrix(self):
        return -self.ml_hessian()

    def ml_fisher_information_matrix(self):
        ranefs = self.mm.random_effects
        n_ranefs = len(ranefs)
        mat = np.zeros((n_ranefs, n_ranefs))

        resid = self.mm.y - self.mm.X * self.beta
        
        def ml_fisher_element(Vinv, V_i, V_j):
            return (0.5 * np.trace(Vinv * V_i * Vinv * V_j))

        for i, ranef_a in enumerate(ranefs):
            dV_dsigma_a=ranef_a.V_i
            for j, ranef_b in enumerate(ranefs):

                if j < i:
                    # Already set when we did the other side of the matrix
                    continue

                dV_dsigma_b=ranef_b.V_i
                element=ml_fisher_element(self.Vinv, dV_dsigma_a, dV_dsigma_b)

                mat[i, j]=element
                mat[j, i]=element

        return np.array(mat)

    def ml_average_information_matrix(self):
        ranefs = self.mm.random_effects
        n_ranefs = len(ranefs)
        mat = np.zeros((n_ranefs, n_ranefs))

        resid = self.mm.y - self.mm.X * self.beta
        
        def ml_fisher_element(Vinv, V_i, V_j):
            return (0.5 * resid.T * (Vinv * V_i * Vinv * V_j) * resid)

        for i, ranef_a in enumerate(ranefs):
            dV_dsigma_a=ranef_a.V_i
            for j, ranef_b in enumerate(ranefs):

                if j < i:
                    # Already set when we did the other side of the matrix
                    continue

                dV_dsigma_b=ranef_b.V_i
                element=ml_fisher_element(self.Vinv, dV_dsigma_a, dV_dsigma_b)

                mat[i, j]=element
                mat[j, i]=element

        return np.array(mat)

class REML(MixedModelLikelihood):

    def set_info(self, info):
        d = {'fs': self.reml_fisher_information_matrix,
             'nr': self.reml_observed_information_matrix,
             'ai': self.reml_average_information_matrix}
    
        self.hessian = d[info]

    def gradient(self):
        "The gradient of the REML function w/r/t each variance component"

        def dREML_dsigma(y, Z, G, P):
            "The REML derivative with regard to a variance component"
            PZGZt = P * Z * G * Z.T
            dl_dsig = -.5 * np.trace(PZGZt) + .5 * (y.T * PZGZt * P * y)
            return matrix.item(dl_dsig)

        ranefs = self.mm.random_effects
        nabla = [dREML_dsigma(self.mm.y, rf.Z, rf.G, self.P) for rf in ranefs]
        return np.array(nabla)

    def loglikelihood(self):
        """
        Returns the restricted loglikelihood for mixed model variance component
        estimation.

        References:

        Harville. 'Maximum Likelihood Approaches to Variance Component
        Estimation and to Related Problems' Journal of the American Statistical
        Association. (1977) (72):258

        Lange, Westlake, & Spence. 'Extensions to pedigree analysis III.
        Variance components by the scoring method.' Ann Hum Genet. (1976).
        39:4,485-491
        DOI: 10.1111/j.1469-1809.1976.tb00156.x

        SAS documentation for PROC MIXED
        """
        y, V, X, P, Vinv = self.mm.y, self.V, self.mm.X, self.P, self.Vinv
        n = X.shape[0]

        rank = np.linalg.matrix_rank(X)
        llik_restricted = -0.5 * (logdet(V.todense())
                                  + logdet(X.transpose() * Vinv * X)
                                  + y.T * P * y
                                  + (n - rank) * l2pi)
        return matrix.item(llik_restricted)


    def reml_hessian(self):
        y, V, X, P, Vinv = self.mm.y, self.V, self.mm.X, self.P, self.Vinv
        ranefs = self.mm.random_effects

        n_ranefs = len(ranefs)
        mat = np.zeros((n_ranefs, n_ranefs))

        def reml_hessian_element(y, P, dV_dsigma_a, dV_dsigma_b):
            common_term=P * dV_dsigma_a * P * dV_dsigma_b
            a=.5 * np.trace(common_term)
            b=y.T * common_term * P * y
            return matrix.item(a - b)

        for i, ranef_a in enumerate(ranefs):
            dV_dsigma_a=ranef_a.V_i
            for j, ranef_b in enumerate(ranefs):

                if j < i:
                    # Already set when we did the other side of the matrix
                    continue

                dV_dsigma_b=ranef_b.V_i
                element=reml_hessian_element(y, P, dV_dsigma_a, dV_dsigma_b)

                mat[i, j]=element
                mat[j, i]=element

        return np.array(mat)


    def reml_observed_information_matrix(self):        
        return -self.reml_hessian()

    def reml_fisher_information_matrix(self):
        y, V, X, P, Vinv = self.mm.y, self.V, self.mm.X, self.P, self.Vinv
        ranefs=self.mm.random_effects

        mat=np.zeros((len(ranefs), len(ranefs)))

        def reml_fisher_element(P, dV_dsigma_a, dV_dsigma_b):
            return .5 * np.trace(P * dV_dsigma_a * P * dV_dsigma_b)

        for i, ranef_a in enumerate(ranefs):
            dV_dsigma_a=ranef_a.V_i

            for j, ranef_b in enumerate(ranefs):
                if j < i:
                    # Already set when we did the other side of the matrix
                    continue

                element=reml_fisher_element(P, dV_dsigma_a, ranef_b.V_i)

                mat[i, j]=element
                mat[j, i]=element

        return mat


    def reml_average_information_matrix(self):
        y, V, X, P, Vinv = self.mm.y, self.V, self.mm.X, self.P, self.Vinv
        ranefs=self.mm.random_effects
        
        mat=np.zeros((len(ranefs), len(ranefs)))

        def reml_average_information_element(y, P, dV_dsigma_a, dV_dsigma_b):
            return .5 * y.T * P * dV_dsigma_a * P * dV_dsigma_b * P * y

        for i, ranef_a in enumerate(ranefs):
            dV_dsigma_a=ranef_a.V_i

            for j, ranef_b in enumerate(ranefs):
                if j < i:
                    # Already set when we did the other side of the matrix
                    continue
                dV_dsigma_b=ranef_b.V_i
                element=reml_average_information_element(y, P,
                                                           dV_dsigma_a,
                                                           dV_dsigma_b)
                mat[i, j]=element
                mat[j, i]=element

        return mat

    def expectation_maximization(self):
        "Performs a round of Expectation-Maximization REML"
        y, P = self.mm.y, self.P

        n = self.mm.nobs()
        coefficients = np.array([
            matrix.item(y.T * P * cov * P * y - np.trace(P * cov))
            for cov in self.mm.covariance_matrices])

        delta = (self.parameters ** 2 / n) * coefficients
        return self.parameters + delta
