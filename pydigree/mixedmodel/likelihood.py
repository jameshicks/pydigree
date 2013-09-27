#!/usr/bin/env python

from itertools import izip
from math import log as ln

import numpy as np
from numpy.linalg import inv,pinv,det

def makeV(covariance_mats, incidence_mats, variance_components, residual_variance,R):
    """
    Returns the V matrix for the restricted loglikelihood function. For a
    decent definition, see

    Reference:
    Lange, Westlake, & Spence. Extensions to pedigree analysis III. Variance components by the scoring method.
    Ann Hum Genet. (1976). 39:4,485-491 DOI: 10.1111/j.1469-1809.1976.tb00156.x
    """
    V = sum(sigma * Z * A * Z.transpose()  \
            for sigma,A,Z in \
            izip(variance_components, covariance_mats,incidence_mats))
    V = V + residual_variance * R 
    return V
def restricted_loglikelihood(y,V,X):
    """
    Returns a value proportional to the restricted loglikelihood for mixed model estimation.

    References:
    
    Harville. 'Maximum Likelihood Approaches to Variance Component Estimation and to Related Problems'
    Journal of the American Statistical Association. (1977) (72):258

    Lange, Westlake, & Spence. Extensions to pedigree analysis III. Variance components by the scoring method.
    Ann Hum Genet. (1976). 39:4,485-491 DOI: 10.1111/j.1469-1809.1976.tb00156.x

    SAS documentation for PROC MIXED
    """
    
    Vinverse = inv(V)
    r = y - X * pinv(X.transpose() * Vinverse * X) * X.transpose() * Vinverse * y
    # This value is only proportional to the restricted likelihood. I'm saving some needless expense by
    # skipping some terms that are constant across all formulations of the variance components.
    # To get the actual estimate of the likelihood, calculate -0.5(w + (n-p)*ln(2pi)), where:
    #   n is the number of rows of X
    #   p is the rank of X.
    #
    # I've left these terms out for two main reasons:
    # 1) They don't change across different values of V. For optimization purposes, they're irrelevant.
    # 2) The term (n-p)ln(2pi) requires numpy.linalg.matrix_rank, which is not found in older
    #    versions of numpy. At some point I might require newer versions of numpy but, because of (1)
    #    I don't think i'll bother.
    llik_restricted = ln(det(V)) + ln(det(X.transpose() * Vinverse * X)) + r.transpose() * Vinverse * r
    return llike_restricted
