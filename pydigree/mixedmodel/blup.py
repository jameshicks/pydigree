"""
Functions for computing Henderson's mixed model equationss
"""

import numpy as np
from itertools import izip
from scipy.linalg import inv, solve


def makeLHS(X, Zlist, Ginvlist, Rinv):
    """
    Makes the left-hand side of the mixed model equations. See
    notes on pydigree.mixedmodel.blup
    """
    predictors = [X] + Zlist
    Ginvlist = [0] + Ginvlist
    cmat = []
    for i, a in enumerate(predictors):
        row = []
        for j, b in enumerate(predictors):
            element = a.transpose() * Rinv * b
            if i == j:
                # This is the diagonal
                element = element + Ginvlist[i]
            row.append(element)
        cmat.append(row)
    return np.bmat(cmat)


def makeRHS(y, X, Zlist, Rinv):
    """
    Makes the right-hand side (RHS) of the mixed model equations. See
    notes on pydigree.mixedmodel.blup
    """
    effects = [X] + Zlist
    m = [[effect.transpose() * Rinv * y] for effect in effects]
    return np.bmat(m)


def blup(y, X, Zlist, covariance_matrices, variance_components, R=None):
    """
    Solves Hendersons mixed model equations for one or more
    uncorrelated random effects.

    Value
    ------
    Returns the column vector corresponding to the best linear unbiased
    estimator for fixed effects and best linear unbiased predictors for
    the random effects.

    Arguements
    ------
    Where: n is the number of individuals to be estimated (including those
    not directly observed for breeding value estimation), m is the number
    of observations, and p is the number of fixed effects, the arguments to
    this function are:

    y:                   a p x 1 matrix of responses
    X:                   an m x p design matrix of fixed effects
    Z:                   a list of n x m design matrices for
                         random effects
    covariance_matrices: a list of n x n covariance matrices
                         for the random effects
    variance_components: a list of variance components associated with
                         each random effect
    R:                   Matrix size n x n of residual errors for GLS
                         approximation (not implemented!)

    Additional information
    -----
    Hendersons Mixed Model equation looks like
    Cp = RHS
        /                                                            \
        |  (X' Rinv X)       (X'  Rinv Z1)         (X'  Rinv Z2)     |
    C = | (Z1' Rinv X)   (Z1' Rinv Z1 + G1inv)     (Z1' Rinv Z2)     |
        | (Z2' Rinv X)       (Z2' Rinv Z1)     (Z2' Rinv Z1 + G2inv) |
        \                                                            /

          /            \      /     \   X,R,y: described in docstring
          | X'  Rinv y |      |  b  |   Zx: incidence matrix of random eff. x
    RHS = | Z1' Rinv y |  p = |  a  |   Gx: covariance matrix for random eff. x
          | Z2' Rinv y |      |  d  |   b: vector of fixed effect estimates
          \            /      \     /   a: vector of BLUPs for random effect 1
                                        d: vector of BLUPs for random effect 2

    """
    if not Zlist:
        raise ValueError('No random effects. Do you want an OLS/GLS solver?')
    m, nfixef = X.shape
    n = Zlist[0].shape[0]
    if R is None:
        Rinv = np.matrix(np.eye(n))
    else:
        Rinv = np.linalg.inv(R)
    residual_variance = np.var(y) - sum(variance_components)
    # Make the covariance matrices into G matrices
    Ginvmats = [1 / v * inv(A.todense())
                for v, A in izip(variance_components,  covariance_matrices)]
    LHS = makeLHS(X, Zlist, Ginvmats, Rinv)
    RHS = makeRHS(y, X, Zlist, Rinv)
    predictions = solve(LHS, RHS)
    return predictions
