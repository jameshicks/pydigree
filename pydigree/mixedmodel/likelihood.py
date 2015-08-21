"""
Functions for computing likelihoods of linear mixed models
"""

from math import log, pi

import numpy as np
from scipy.sparse import bsr_matrix, issparse
from scipy.linalg import pinv, inv
from scipy import matrix
np.seterr(invalid='ignore')


l2pi = log(2 * pi)


def logdet(M):
    """ Returns the (positive) log determinant of a matrix. """
    sign, logdet = np.linalg.slogdet(M.todense() if issparse(M) else M)
    return logdet


def makeR(y, X, V=None, Vinv=None):
    """ Makes the R matrix commonly found in REML estimation """
    if V is None and Vinv is None:
        raise ValueError('Variance matrix not specified')
    elif Vinv is None and V is not None:
        Vinv = makeVinv(V)
    return y - X * pinv(X.transpose() * Vinv * X) * X.transpose() * Vinv * y


def makeP(X, Vinv):
    """ Makes the P matrix commonly found in mixed model estimation """
    return Vinv - Vinv * X * inv(X.T * Vinv * X) * X.T * Vinv


def makeVinv(V):
    if issparse(V):
        V = V.todense()
    return bsr_matrix(inv(V))


def full_loglikelihood(y, V, X, P=None, Vinv=None):
    """
    Returns the full loglikelihood of a mixed model

    Ref: SAS documentation for PROC MIXED
    """
    if Vinv is None:
        Vinv = makeVinv(V)
    R = makeR(y, X, Vinv=Vinv)
    n = X.shape[0]
    llik = -0.5 * (logdet(V.todense()) + R.transpose() * Vinv * R + n * l2pi)
    return matrix.item(llik)


def restricted_loglikelihood(y, V, X, P=None, Vinv=None):
    """
    Returns the restricted loglikelihood for mixed model variance component
    estimation.

    References:

    Harville. 'Maximum Likelihood Approaches to Variance Component Estimation
    and to Related Problems' Journal of the American Statistical Association.
    (1977) (72):258

    Lange, Westlake, & Spence. 'Extensions to pedigree analysis III. Variance
    components by the scoring method.' Ann Hum Genet. (1976). 39:4,485-491
    DOI: 10.1111/j.1469-1809.1976.tb00156.x

    SAS documentation for PROC MIXED
    """
    if Vinv is None:
        Vinv = makeVinv(V)
    if P is None:
        P = makeP(X, Vinv=Vinv)
    n = X.shape[0]
    rank = np.linalg.matrix_rank(X)
    llik_restricted = -0.5 * (logdet(V.todense())
                              + logdet(X.transpose() * Vinv * X)
                              + y.T * P * y
                              + (n - rank) * l2pi)
    return matrix.item(llik_restricted)


def reml_gradient(y, X, V, ranefs, P=None, Vinv=None):
    if Vinv is None:
        Vinv = makeVinv(V)
    if P is None:
        P = makeP(X, Vinv)
    nabla = [dREML_dsigma(y, rf.Z, rf.G, P) for rf in ranefs]
    return np.array(nabla)


def dREML_dsigma(y, Z, G, P):
    "The REML derivative of V with regard to sigma"
    PZGZt = P * Z * G * Z.T
    dl_dsig = -.5 * np.trace(PZGZt) + .5 * (y.T * PZGZt * P * y)
    return matrix.item(dl_dsig)


def reml_hessian_element(y, P, dV_dsigma_a, dV_dsigma_b):
    common_term = P * dV_dsigma_a * P * dV_dsigma_b
    a = .5 * np.trace(common_term)
    b = y.T * common_term * P * y
    return matrix.item(a - b)


def reml_hessian(y, X, V, ranefs, P=None, Vinv=None):
    if Vinv is None:
        Vinv = makeVinv(V)
    if P is None:
        P = makeP(X, Vinv)

    n_ranefs = len(ranefs)
    mat = np.zeros((n_ranefs, n_ranefs))

    for i, ranef_a in enumerate(ranefs):
        dV_dsigma_a = ranef_a.V_i
        for j, ranef_b in enumerate(ranefs):

            if j < i:
                # Already set when we did the other side of the matrix
                continue

            dV_dsigma_b = ranef_b.V_i
            element = reml_hessian_element(y, P, dV_dsigma_a, dV_dsigma_b)

            mat[i, j] = element
            mat[j, i] = element

    return np.array(mat)


def reml_observed_information_matrix(y, X, V, ranefs, P=None, Vinv=None):
    return -reml_hessian(y, X, V, ranefs, P, Vinv)


def reml_fisher_element(P, dV_dsigma_a, dV_dsigma_b):
    return .5 * np.trace(P * dV_dsigma_a * P * dV_dsigma_b)


def reml_fisher_information_matrix(y, X, V, ranefs, P=None, Vinv=None):
    if Vinv is None:
        Vinv = makeVinv(V)
    if P is None:
        P = makeP(X, Vinv)

    mat = np.zeros((len(ranefs), len(ranefs)))

    for i, ranef_a in enumerate(ranefs):
        dV_dsigma_a = ranef_a.V_i

        for j, ranef_b in enumerate(ranefs):
            if j < i:
                # Already set when we did the other side of the matrix
                continue

            element = reml_fisher_element(P, dV_dsigma_a, ranef_b.V_i)

            mat[i, j] = element
            mat[j, i] = element

    return mat


def reml_average_information_element(y, P, dV_dsigma_a, dV_dsigma_b):
    return y.T * P * dV_dsigma_a * P * dV_dsigma_b * P * y


def reml_average_information_matrix(y, X, V, ranefs, P=None, Vinv=None):
    if Vinv is None:
        Vinv = makeVinv(V)
    if P is None:
        P = makeP(X, Vinv)
    mat = np.zeros((len(ranefs), len(ranefs)))
    for i, ranef_a in enumerate(ranefs):
        dV_dsigma_a = ranef_a.V_i

        for j, ranef_b in enumerate(ranefs):
            if j < i:
                # Already set when we did the other side of the matrix
                continue
            dV_dsigma_b = ranef_b.V_i
            element = reml_average_information_element(y, P,
                                                       dV_dsigma_a,
                                                       dV_dsigma_b)
            mat[i, j] = element
            mat[j, i] = element

    return mat
