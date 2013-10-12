"""
Functions for computing likelihoods of linear mixed models
"""

from math import log, pi

import numpy as np
from scipy.sparse import csc_matrix
from scipy.linalg import pinv, inv
from scipy import matrix
np.seterr(invalid='ignore')


l2pi = log(2 * pi)


def logdet(M):
    """ Returns the (positive) log determinant of a matrix. """
    sign, logdet = np.linalg.slogdet(M)
    return logdet


def makeP(y, X, V=None, Vinv=None):
    """ Makes the P matrix commonly found in REML estimation """
    if V is None and Vinv is None:
        raise ValueError('Variance matrix not specified')
    elif Vinv is None and V is not None:
        Vinv = csc_matrix(inv(V.todense()))
    return y - X * pinv(X.transpose() * Vinv * X) * X.transpose() * Vinv * y


def full_loglikelihood(y, V, X):
    """
    Returns the full loglikelihood of a mixed model

    Ref: SAS documentation for PROC MIXED
    """
    Vinv = csc_matrix(inv(V.todense()))
    P = makeP(y, X, Vinv=Vinv)
    n = X.shape[0]
    llik = -0.5 * (logdet(V.todense()) + P.transpose() * Vinv * P + n * l2pi)
    return matrix.item(llik)


def restricted_loglikelihood(y, V, X):
    """
    Returns a value proportional to the restricted loglikelihood for mixed
    model estimation.

    References:

    Harville. 'Maximum Likelihood Approaches to Variance Component Estimation
    and to Related Problems' Journal of the American Statistical Association.
    (1977) (72):258

    Lange, Westlake, & Spence. 'Extensions to pedigree analysis III. Variance
    components by the scoring method.' Ann Hum Genet. (1976). 39:4,485-491
    DOI: 10.1111/j.1469-1809.1976.tb00156.x

    SAS documentation for PROC MIXED
    """

    Vinv = csc_matrix(inv(V.todense()))
    P = makeP(y, X, Vinv=Vinv)
    # This value is only proportional to the restricted likelihood.
    # I'm saving some needless expense by skipping some terms that are constant
    # across all formulations of the variance components.
    # To get the actual estimate of the likelihood, calculate:
    #   -0.5(w + (n-p)*ln(2pi)), where:
    #   n is the number of rows of X
    #   p is the rank of X.
    #
    # I've left these terms out for two main reasons:
    # 1) They don't change across different values of V. For optimization
    #    purposes, they're irrelevant.
    # 2) The term (n-p)ln(2pi) requires numpy.linalg.matrix_rank, which is
    #    not found in older versions of numpy. At some point I might require
    #    newer versions of numpy but, because of (1) I don't think I'll bother.
    llik_restricted = logdet(V.todense()) + logdet(X.transpose() * Vinv * X) \
        + P.transpose() * Vinv * P
    return matrix.item(llik_restricted)
