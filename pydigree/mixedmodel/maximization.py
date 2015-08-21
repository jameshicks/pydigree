import numpy as np
from numpy.linalg import inv
from scipy.sparse import csc_matrix

from likelihood import reml_gradient
from likelihood import reml_observed_information_matrix
from likelihood import reml_fisher_information_matrix
from likelihood import reml_average_information_matrix
from likelihood import restricted_loglikelihood
from likelihood import makeP, makeVinv


def iterative_scoring_method(mm, starts, method='Fisher',
                             tol=1e-4, verbose=False):
    """
    Updates variance components for a mixed model in an iterative scheme.
    """
    if method.lower() in {'newton-raphson', 'newton', 'nr'}:
        information_mat = reml_observed_information_matrix
    elif method.lower() in {'fisher scoring', 'fisher', 'fs'}:
        information_mat = reml_fisher_information_matrix
    elif method.lower() in {'average information', 'aireml', 'ai'}:
        information_mat = reml_average_information_matrix
    else:
        raise ValueError('Unknown maximization method')

    if verbose:
        print 'Maximizing model with {} method'.format(method)

    vcs = np.array(starts)

    i = 0
    while True:
        V = mm._makeV(vcs.tolist())

        # Complicated things we only want to calculate once
        Vinv = makeVinv(V)
        P = makeP(mm.X, Vinv)

        # Makes the information matrix and gradient, then performs an iteration
        grad = reml_gradient(mm.y, mm.X, V, mm.random_effects, P=P, Vinv=Vinv)
        mat = information_mat(mm.y, mm.X, V, mm.random_effects, P=P, Vinv=Vinv)
        delta = scoring_iteration(mat, grad)

        
        vcs = vcs - delta
        llik = restricted_loglikelihood(mm.y, V, mm.X, P, Vinv)
        
        if verbose:
            print i, llik, vcs
        
        if (np.abs(delta) < tol).all():
            break

        i += 1

    mm.set_variance_components(vcs.tolist())


def scoring_iteration(info_mat, gradient):
    info_mat = np.matrix(info_mat)
    gradient = np.matrix(gradient)
    return -1.0 * np.array(info_mat.I * gradient.T).T[0]
