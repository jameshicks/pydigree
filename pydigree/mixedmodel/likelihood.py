#!/usr/bin/env python

from itertools import izip
from math import log as ln

import numpy as np
from numpy.linalg import inv,pinv,det

def makeV(covariance_mats, variance_components, residual_variance):
    """
    Returns the V matrix for the restricted loglikelihood function. For a
    decent definition, see

    Amos. 'Robust Variance-Components Approach for Assessing Genetic Linkage in Pedigrees'
    Am J Hum Genet (1994) 54:535-543
    """
    V = sum(partial_variance * covariance_matrix \
            for partial_variance, covariance_matrix in \
            izip(variance_components, covariance_mats))
    variance = sum(variance_components) + residual_variance
    np.fill_diagonal(V,variance)
    return V
def restricted_loglikelihood(y,v,x):
    """
    Returns the restricted loglikelihood for mixed model estimation.

    Reference:
    Harville. 'Maximum Likelihood Approaches to Variance Component Estimation and to Related Problems'
    Journal of the American Statistical Association. (1977) (72):258
    "
    vinv = inv(v)
    p = vinv - vinv * x * pinv(x.T * vinv * x) * x.T * vinv
    return - 0.5 * ( ln(det(v) + ln( x.T * vinv * x) + y.T * p * y )
