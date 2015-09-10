import numpy as np
from numpy.linalg import inv, LinAlgError
from scipy.sparse import csc_matrix
from scipy import matrix

from likelihood import reml_gradient
from likelihood import reml_observed_information_matrix
from likelihood import reml_fisher_information_matrix
from likelihood import reml_average_information_matrix
from likelihood import restricted_loglikelihood
from likelihood import makeP, makeVinv


def is_invertible(m):
    return np.isfinite(np.linalg.cond(m.todense()))


def iterative_scoring_method(mm, starts, method='Fisher', maxiter=250,
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

        if not np.isfinite(delta).all():
            raise LinAlgError('NaNs in scoring update')

        vcs = vcs - delta
        llik = restricted_loglikelihood(mm.y, V, mm.X, P, Vinv)

        if verbose:
            print i, llik, vcs

        if (np.abs(delta) < tol).all():
            break

        i += 1
        if i > maxiter:
            raise LinAlgError('Ran out of scoring iterations')

    mm.set_variance_components(vcs.tolist())


def scoring_iteration(info_mat, gradient):
    try:
        info_mat = np.matrix(info_mat)
        gradient = np.matrix(gradient)
        return -1.0 * np.array(info_mat.I * gradient.T).T[0]
    except LinAlgError:
        raise LinAlgError('Information matrix not invertible!')


def expectation_maximization_reml(mm, starts, maxiter=100, tol=1e-4,
                                  verbose=False):
    i = 0

    if verbose:
        print 'Maximizing model by expectation-maximization'

    n = mm.nobs()
    y = mm.y
    vcs = np.array(starts)
    while True:
        V = mm._makeV(vcs.tolist())

        # Complicated things we only want to calculate once
        Vinv = makeVinv(V)
        P = makeP(mm.X, Vinv)

        coefficients = np.array([
            matrix.item(y.T * P * cov * P * y - np.trace(P * cov))
            for cov in mm.covariance_matrices])

        delta = (vcs ** 2 / n) * coefficients
        vcs += delta

        llik = restricted_loglikelihood(mm.y, V, mm.X, P, Vinv)

        if (np.abs(delta) < tol).all():
            break

        if verbose:
            print i, llik, vcs

        if i > maxiter:
            raise LinAlgError('Ran out of scoring iterations')

    mm.set_variance_components(vcs.tolist())
