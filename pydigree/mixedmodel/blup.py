#!/usr/bin/env python

import numpy as np
from itertools import izip

def makeC(X,Zlist,Ginvlist,Rinv):
    predictors = [X] + Zlist
    Ginvlist = [0]+Ginvlist
    cmat = []
    for i,a in enumerate(predictors):
        row = []
        for j,b in enumerate(predictors):
            element = a.transpose() * Rinv * b
            if i == j:
                # This is the diagonal
                element = element + Ginvlist[i]
            row.append(element)
        cmat.append(row)
    return np.bmat(cmat)
def makeE(y,X,Zlist,Rinv):
    m = []
    effects = [X] + Zlist
    for effect in effects:
        m.append( [ effect.transpose() * Rinv * y ] )
    return np.bmat(m)

def blup(y,X,Zlist,covariance_matrices,variance_components,R=None):
    """
    Solves mixed model equations for one or more uncorrelated random effects.

    Where: n is the number of individuals to be estimated (including those not
    directly observed for breeding value estimation), m is the number of
    observations, and p is the number of fixed effects, the arguments to
    this function are:
    
    y: a p x 1 matrix of responses  
    X: an m x p design matrix of fixed effects
    Z: a list of n x m design matrices for random effects
    covariance_matrices: a list of n x n covariance matrices for the random effects
    variance_components: a list of scalar ratios in the form sigma_e / sigma_a 
    R: Matrix size n x n of residual errors for GLS approximation (not implemented!)
    """
    if R is None:
        R = np.matrix(np.eye(X.shape[0]))
        Rinv = R
    else: 
        return NotImplementedError('Structured R matrices are not supported at this time')
    residual_variance = 1 - sum(variance_components)
    # Make the covariance matrices into G matrices
    Ginvmats = [ (residual_variance/lambd) * np.linalg.inv(covmat) for lambd,covmat in izip(variance_components,covariance_matrices)]

    # Hendersons Mixed Model equation looks like
    # Cp = E
    # Where p is a column vector of BLUEs for fixed effect and BLUPs for random effects and
    #     /                                                            \
    #     |  (X' Rinv X)       (X' Rinv Z1)          (X' Rinv Z2)      |  
    # C = | (Z1' Rinv X)   (Z1' Rinv Z1 + G1inv)     (Z1' Rinv Z2 )    |  
    #     | (Z2' Rinv X)       (Z2' Rinv Z1)     (Z2' Rinv Z1 + G2inv) |     
    #     \                                                            /
    #
    #     /            \
    #     | X'  Rinv y |                               
    # E = | Z1' Rinv y |
    #     | Z2' Rinv y |
    #     \            /
    C = makeC(X,Zlist,Ginvmats,Rinv)
    E = makeE(y,X,Zlist,Rinv)
    predictions = np.linalg.solve(C,E)
    return predictions
