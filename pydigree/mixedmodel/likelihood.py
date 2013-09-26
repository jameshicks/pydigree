#!/usr/bin/env python

from itertools import izip
from math import log as ln

import numpy as np
from numpy.linalg import inv,pinv,det

def makeV(covariance_mats, variance_components, residual_variance):
    """
    Returns the V matrix for the restricted loglikelihood function. For a
    decent definition, see

    Reference:
    Lange, Westlake, & Spence. Extensions to pedigree analysis III. Variance components by the scoring method.
    Ann Hum Genet. (1976). 39:4,485-491 DOI: 10.1111/j.1469-1809.1976.tb00156.x
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

    References:
    Harville. 'Maximum Likelihood Approaches to Variance Component Estimation and to Related Problems'
    Journal of the American Statistical Association. (1977) (72):258

    Lange, Westlake, & Spence. Extensions to pedigree analysis III. Variance components by the scoring method.
    Ann Hum Genet. (1976). 39:4,485-491 DOI: 10.1111/j.1469-1809.1976.tb00156.x
    """
    
    vinv = inv(v)
    p = vinv - vinv * x * pinv(x.T * vinv * x) * x.T * vinv
    return -0.5 * ( ln(det(v) + ln( x.T * vinv * x) + y.T * p * y ))
                    
