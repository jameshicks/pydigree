import numpy as np
from numpy.linalg import LinAlgError

from pydigree.stats.mathfuncs import is_positive_definite, grid

# from .likelihood import makeP


class MLEResult(object):
    " An object representing the result of maximization of a mixed model "

    def __init__(self, parameters, restricted_loglikelihood, method,
                 jacobian=None, hessian=None, full_loglikelihood=None):
        self.restricted = True
        self.restricted_loglikelihood = restricted_loglikelihood
        self.full_loglikelihood = full_loglikelihood
        self.method = method
        self.parameters = parameters
        self.jacobian = jacobian
        self.hessian = hessian


def newtonlike_maximization(mm, likelihood, maxiter=250, tol=1e-4, 
                            scoring=10, verbose=False):
    """
    Updates variance components for a linear mixed model by an
    iterative scheme to find the restricted maximum likelihood estimates
    of the variance components in the model.

    Models are fit by iterating the following equation:
        theta_(i+1) = theta_i - inverse(J(theta_i)) * S(theta_i)
    Where:
        theta_i = A vector of the estimated variance components at
            iteration i
        S(theta): The score function (gradient) of the REML loglikelihood
        J(theta): The information matrix (the matrix of second derivatives
            of the REML loglikelihood with regard to each of the variance
            compents in theta)

    In all optimization schemes in this function S(theta) is the gradient
    of the REML loglikelihood function evaluated at the current estimates
    of theta.

    The information matrix J(theta) is a q x q square matrix for the q
    predictors in the model. There are a few information matrices available
    for use, specified by the argument `method`. 'Newton-Raphson' uses the
    observed information matrix (the negative Hessian). This is the most
    complicated to calculate, is very sensitive to starting values, and can
    occasionally have numerical issues. The method value 'Fisher scoring' uses
    the Fisher Information Matrix (expected value of negative hessian),
    which is simpler to calculate. This is a common way of fitting mixed model
    and is the default method for this optimizer. The last ('Average
    Information') uses the average of the Fisher Information matrix and
    Observed Information matrix. This is a common approach in the animal
    breeding literature. Averaging the two results in the elimination
    of a time consuming trace term, making this the fastest method in terms
    of time per iteration, though it may require a few more iterations than
    Fisher scoring or Newton-Raphson.

    When the change in the proportion of variance explained by each variance
    component after iteration falls below `tol` for every variance component,
    iteration stops and the estimated variance components are returned. Setting
    the tolerance based on the proportion has the effect of standardizing
    tolerances over any amount of variance.

    Occasionally, Newton-type algorithms will push the estimated values of
    the variance components to invalid levels (i.e. below zero, above the
    total variance of the outcome variable). Outside the valid range, the
    loglikelihood surface becomes ill-conditioned and the optimizer may not
    return back to valid parameter estimates. This is especially true for
    Fisher scoring and AIREML when the true value of a variance component is
    close to the border of valid estimates. The information matrices used by
    Fisher Scoring and AIREML, being approximations to the Hessian, can put the
    iteration estimates of the variance components outside the valid space. The
    parameter `constrained` enforces validity of variance component estimates
    in two ways: likelihoods must be monotonically increasing, and variance
    component estimates must be in the valid range. If the change in estimates
    violates either of these, a line search between that change and changes in
    the same direction but exponentially decreasing in magnitude is performed
    until a valid set of estimates is met.

    The scoring argument allows you to run a few iterations of Fisher scoring
    or AI-REML to get close to the maximum, then switch over to Newton-Raphson
    to end quicker.

    :param mm: model to be maximized
    :param likelihood: a Likelihood object
    :param maxiter: The maximum number of iterations of scoring before raising
        an error
    :param tol: The minimum amount of change in the proportion of variance by 
        any of the variance components to continue iterating.
    :param scoring: Number of iterations of Fisher Scoring or AI-REML before
        switching to Newton-Raphson. If already using Newton-Raphson,
        this argument has no effect.
    :param verbose: Print likelihood, variance component, and relative variance
        estimates at each iteration. Useful for debugging or watching the
        direction the optimizer is taking.

    :type mm: MixedModel
    :type tol: float
    :type likelihood: Likelihood
    :type maxiter: integer
    :type scoring: integer
    :type verbose: bool

    Returns: An MLE object of the variance components at the MLE
    """
    if verbose:
        print(('Maximizing model by {}'.format(likelihood.method)))

    vcs = np.array(likelihood.parameters)

    llik = likelihood.loglikelihood()

    if verbose:
        print(('{} {} {} {}'.format(0, llik, vcs, vcs / vcs.sum())))

    for i in range(maxiter):
        if (i - 1) == scoring:
            likelihood.set_info('nr')

        # Make the information matrix and gradient
        grad = likelihood.gradient()
        mat = likelihood.info_matrix()
        delta = scoring_iteration(mat, grad)

        if not is_positive_definite(mat):
            raise LinAlgError('Information matrix not positive definite')

        if np.linalg.cond(mat) > 1e4:
            raise LinAlgError(
                'Condition number of information matrix too high')
        
        if not np.isfinite(delta).all():
            raise LinAlgError('NaNs in scoring update')


        new_vcs = vcs - delta

        likelihood.set_parameters(new_vcs)
        new_llik = likelihood.loglikelihood()

        if new_vcs.sum() / mm._variance_after_fixefs() > 10:
            raise LinAlgError('Optimizer left parameter space')
        
        relative_changes = (new_vcs / new_vcs.sum()) - (vcs / vcs.sum())

        if verbose:
            print((i+1, new_llik, new_vcs, \
                new_vcs / new_vcs.sum(), relative_changes))

        if (abs(relative_changes) < tol).all():
            mle = MLEResult(new_vcs.tolist(), new_llik, likelihood.method,
                            jacobian=grad, hessian=mat)
            return mle
        vcs = new_vcs
        llik = new_llik

    raise LinAlgError('Ran out of scoring iterations')


def scoring_iteration(info_mat, gradient):
    """
    Performs an iteration for a Newton-type maximization algorithm

    :param info_mat: A matrix of second derivatives (or approximations of it) 
        at the current parameter estimates
    :param gradient: A vector containing the gradient at the current parameter
        estimates

    :type info_mat: 2d numpy array
    :type gradient: 1d numpy array

    :returns: change in parameters for the current iteration.
    :rtype: numpy array

    Raises: LinAlgError if the information matrix is singular

    """
    try:
        info_mat = np.matrix(info_mat)
        gradient = np.matrix(gradient)
        return -1.0 * np.array(info_mat.I * gradient.T).T[0]
    except LinAlgError:
        raise LinAlgError('Information matrix not invertible!')


def expectation_maximization(mm, likelihood, starts=None, maxiter=10000, 
                             tol=1e-4, return_after=1e30, verbose=False):
    '''
    Maximizes a linear mixed model by Expectation-Maximization

    Formulas for EM-REML are given in Lynch & Walsh, Ch 27, Example 5 (pg. 799)

    Unlike the Newton-type algorithms, EM only makes use of the first 
    derivative of the loglikelihood function. The presence of second 
    derivatives means that Newton-type maximization will converge very quickly,
    since it works on a better approximation of the likelihood surface. 

    EM tends to converge VERY slowly, because the changes at every step
    are so small. For example, a model that took 3 iterations/0m3.927s to 
    converge with AI-REML took 52 iterations/0m32.803s with EM-REML. 
    Individual EM iterations are relatively fast because you don't have to
    compute the Hessian (or an approximation of it). But since you have to 
    invert the variance-covariance matrix of observations each iteration
    regardless, time adds up quickly.

    However, EM has the nice property that it monotonically converges to
    a maximum, and avoids the parameter esimtate out-of-range problems 
    occasionally found with Newton-type methods and have to be remedied as 
    a special case. 

    :param mm: a model to be maximized
    :param likelihood: A likelihood object
    :param starts: starting values
    :param maxiter: The maximum number of iterations of scoring before raising 
        an error
    :param tol: The minimum amount of change in the proportion of variance by 
        any of the variance components to continue iterating. 
    :param verbose: Print likelihood, variance component, and relative variance 
        estimates at each iteration. Useful for debugging or watching the 
        direction the optimizer is taking.
    :param return_after: return estimates after n iterations, regardless of 
        whether or not the model converged.
    :type mm: MixedModel
    :type maxiter: integer
    :type tol: float
    :type verbose: bool
    :type starts: numpy array

    :returns:  variance components at the MLE
    :rtype: MLEResult
    '''
    i = 0

    if verbose:
        print('Maximizing model by Expectation-Maximization')

    vcs = likelihood.parameters if starts is None else starts 
    llik = likelihood.loglikelihood()
    while i < return_after:
        if verbose:
            print((i, llik, vcs))

        new_vcs = likelihood.expectation_maximization()
        likelihood.set_parameters(new_vcs)

        llik = likelihood.loglikelihood()

        if (np.abs(new_vcs - vcs) < tol).all():
            break

        vcs = new_vcs

        if i > maxiter:
            raise LinAlgError('Ran out of scoring iterations')

        i += 1
    mle = MLEResult(vcs, llik, 'Expectation-Maximization')
    return mle


# def minque(mm, starts=None, value=0, maxiter=200, tol=1e-4,
#            verbose=False, return_after=1e300, return_vcs=False):
#     """ 
#     MINQUE (MInimum Norm Quadratic Unbiansed Estimation). Only used for 
#     historical purposes or getting starting variance components for another
#     maximization scheme.

#     MINQUE gets variance component estimates by solving the equation Cz=t

#     For d random effects 
#     z is a vector of variance compnents
#     C is a dxd matrix with element  C_ij trace(P * V_i * P * V_j)
#     t is a column vector with row element i = y' * P * V_i * P * y

#     M = I_n - (1/n) * ONES_n * ONES_n'
#     (Ones_n is a row vector of all ones)


#     Useful reference: 
#     J.W. Keele & W.R. Harvey (1988) "Estimation of components of variance and
#     covariance by symmetric difference squaredand minimum norm quadratic 
#     unbiased estimation: a comparison" Journal of Animal Science
#     Vol 67. No.2 p348-356
#     doi:10.2134/jas1989.672348x
#     """
#     d = len(mm.random_effects)  # the number of random effects
#     if verbose:
#         print('Maximizing model by MINQUE')

#     if starts is not None:
#         weights = np.array(starts)

#     elif value == 0:
#         # MINQUE(0)
#         weights = np.zeros(d)
#         weights[-1] = 1

#     elif value == 1:
#         # MINQUE(1)
#         weights = np.ones(d)

#     vcs = np.var(mm.y - mm.X * mm.beta) * weights
#     n = mm.nobs()

#     y = mm.y

#     if verbose:
#         print(vcs)
#     for i in range(maxiter):

#         if i + 1 > return_after:
#             return vcs

#         V = mm._makeV(vcs=vcs.tolist())
#         Vinv = makeVinv(V)
#         P = makeP(mm.X, Vinv)

#         t = [matrix.item(y.T * P * ranef.V_i * P * y)
#              for ranef in mm.random_effects]
#         t = np.matrix(t).T

#         # Make C
#         C = []
#         for ranef_i in mm.random_effects:
#             row = [np.trace(P * ranef_i.V_i * P * ranef_j.V_i)
#                    for ranef_j in mm.random_effects]
#             C.append(row)
#         C = np.matrix(C)
#         new_vcs = scipy.linalg.solve(C, t).T[0]

#         delta = (new_vcs / new_vcs.sum()) - (vcs / vcs.sum())
#         llik = restricted_loglikelihood(mm.y, V, mm.X, P, Vinv)

#         if all(delta < tol):
#             if return_vcs:
#                 return new_vcs
#             mle = MLEResult(new_vcs, llik, 'MINQUE')
#             return mle

#         if verbose:
#             print((i, llik, vcs))
#         vcs = new_vcs
#         weights = vcs

def grid_search(mm, likelihood, nevals=10, oob=False):
    ''' 
    Grid searches a likelihood function over a range of variance 
    component values

    :param mm: A model to be maximized
    :param likelihood: a MixedModelLikelihood object (i.e. ML, or REML)
    :param nevals: number of evaluations over each component
    :param oob: Evaluate out of bounds arguments (e.g. sum(vcs) > var(y))

    :type mm: MixedModel
    :type nevals: integer
    :type oob: bool
    
    :returns: variance components at the best node in the grid
    :rtype: MLEResult
    '''
    totvar =  mm.y.var()
    
    if not oob:
        pred = lambda *x: sum(x) <= totvar
    else:
        pred = None


    low, high = 0, totvar

    best = None, -np.inf

    def likefunc(*args):
        likelihood.set_parameters(args)
        return likelihood.loglikelihood()

    for test in grid(likefunc, len(mm.random_effects), low, high, nevals, predicate=pred):

        if test[1] > best[1]:
            best = test

    mle = MLEResult(best[0], best[1], 'GRID')
    return mle