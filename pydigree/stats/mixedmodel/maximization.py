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


def newtonlike_maximization(mm, starts, method='Fisher', maxiter=250,
                            tol=1e-4, constrained=True, verbose=False):
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
    close to the border of valid estimates. The parameter `constrained`
    enforces validity of variance component estimates in two ways: likelihoods
    must be monotonically increasing, and variance component estimates must be
    in the valid range. If the change in estimates violates either of these, a
    line search between that change and changes in the same direction but
    exponentially decreasing in magnitude is performed until a valid set of 
    estimates is met. 

    Arguments:
    mm: a MixedModel object to be maximized
    starts: Starting values for the variance components
    method: The method to use to fit the model. 
        * Options: 'Newton-Raphson', 'Fisher Scoring', 'Average Information'
    maxiter: The maximum number of iterations of scoring before raising 
        an error
    tol: The minimum amount of change in the proportion of variance by any of 
        the variance components to continue iterating. 
    constrained: Force optimizer to keep variance component estimates
        in the range 0 <= vc <= variance(y), by performing a line search
    verbose: Print likelihood, variance component, and relative variance 
        estimates at each iteration. Useful for debugging or watching the 
        direction the optimizer is taking.

    Returns: A numpy array of the variance components at the MLE
    """
    if method.lower() in {'newton-raphson', 'newton', 'nr'}:
        information_mat = reml_observed_information_matrix
        method = 'Newton-Raphson'
    elif method.lower() in {'fisher scoring', 'fisher', 'fs'}:
        information_mat = reml_fisher_information_matrix
        method = 'Fisher Scoring'
    elif method.lower() in {'average information', 'aireml', 'ai'}:
        information_mat = reml_average_information_matrix
        method = 'Average Information REML'
    else:
        raise ValueError('Unknown maximization method')

    if verbose:
        print 'Maximizing model by {}'.format(method)

    vcs = np.array(starts)

    # Get the loglikelihood at the start
    V = mm._makeV(vcs.tolist())
    Vinv = makeVinv(V)
    P = makeP(mm.X, Vinv)
    llik = restricted_loglikelihood(mm.y, V, mm.X, P, Vinv)

    if verbose:
        print '{} {} {} {}'.format(0, llik, vcs, vcs / vcs.sum())

    for i in xrange(maxiter):
        # Make the information matrix and gradient
        grad = reml_gradient(mm.y, mm.X, V, mm.random_effects, P=P, Vinv=Vinv)
        mat = information_mat(mm.y, mm.X, V, mm.random_effects, P=P, Vinv=Vinv)
        delta = scoring_iteration(mat, grad)

        if not np.isfinite(delta).all():
            raise LinAlgError('NaNs in scoring update')

        # Newton-Raphson type optimization methods like these sometimes
        # rattle around a bit, and overshoot the boundary the valid parameter
        # space (i.e. all variance components meet 0 <= sigma_i <= sigma_total)
        # and then run out of control in the ill-conditioned space outside the
        # valid parameters (e.g Hessians with condition numbers in the 1e200
        # range).
        #
        # This is especially true for AI-REML in the case that the true value
        # of a variance components is close to the boundary of valid values.
        # This is covered in detail in:
        #
        # Mishchenko, Holmgren, & Ronnegard (2007).
        # "Newton-Type Methods for REML Estimation in Genetic Analysis
        #   of Quantitative Traits"  arXiv:0711.2619 [q-bio.OT]
        #
        # They suggest a line search to shrink the step size (alpha)
        # when the loglikelihood starts decreasing. The step size shrinks
        # until it finds a step that improves the loglikelihood. This
        # is what we're doing here. Additionally we go into the line-search
        # portion
        #
        # Most of the time we'll only run one iteration of this loop, and break
        # out of it, but when you need it, you need it.
        new_vcs = np.zeros(vcs.shape[0])
        for n in xrange(25):
            alpha = 2 ** (-n)
            new_vcs = vcs - alpha * delta

            if constrained and (any(new_vcs < 0) or
                                any(new_vcs > mm.y.var())):
                # Don't bother evaluating variance components
                # if theyre not valid
                continue
            # If we're not changing the the parameters in any meaningful
            # way, we can leave too, beacuse we've probably found a maximum
            relative_changes = abs(delta) / vcs.sum()
            if (abs(relative_changes) < tol).all():
                break

            V = mm._makeV(new_vcs.tolist())
            Vinv = makeVinv(V)
            P = makeP(mm.X, Vinv)

            new_llik = restricted_loglikelihood(mm.y, V, mm.X, P, Vinv)

            # If we have an improvement in loglikelihood, we can leave the
            # search
            improvement = new_llik - llik
            if constrained and improvement > 0:
                break

            
        else:
            # If we've shrunk the change in variance components down by
            # factor of 2**(-25) = 2.98-08 and we still haven't gotten a
            # an improvement, we're probably in a maximum, so we'll step back
            # to our last variance component values
            new_vcs = vcs

            V = mm._makeV(new_vcs.tolist())
            Vinv = makeVinv(V)
            P = makeP(mm.X, Vinv)

            new_llik = restricted_loglikelihood(mm.y, V, mm.X, P, Vinv)
            grad = reml_gradient(
                mm.y, mm.X, V, mm.random_effects, P=P, Vinv=Vinv)
            mat = information_mat(
                mm.y, mm.X, V, mm.random_effects, P=P, Vinv=Vinv)
            delta = scoring_iteration(mat, grad)

        if new_vcs.sum() / np.var(mm.y) > 10:
            raise LinAlgError('Optimizer left parameter space')

        if verbose:
            print i+1, new_llik, new_vcs, new_vcs / new_vcs.sum()

        relative_changes = abs(delta) / vcs.sum()
        if (abs(relative_changes) < tol).all():
            return new_vcs.tolist()

        vcs = new_vcs
        llik = new_llik

    raise LinAlgError('Ran out of scoring iterations')


def scoring_iteration(info_mat, gradient):
    """
    Performs an iteration for a Newton-type maximization algorithm

    Arguments:
    info_mat: A matrix of second derivatives (or approximations of it) at 
        the current parameter estimates
    gradient: A vector containing the gradient at the current parameter
        estimates

    Returns: A numpy array of the change in parameters for the current 
        iteration.

    Raises: LinAlgError if the information matrix is singular

    """
    try:
        info_mat = np.matrix(info_mat)
        gradient = np.matrix(gradient)
        return -1.0 * np.array(info_mat.I * gradient.T).T[0]
    except LinAlgError:
        raise LinAlgError('Information matrix not invertible!')


def expectation_maximization_reml(mm, starts, maxiter=10000, tol=1e-4,
                                  verbose=False, vcs_after_maxiter=False):
    '''
    Maximizes a linear mixed model by Expectation-Maximization

    Formulas for EM-REML are given in Lynch & Walsh, Ch 27, Example 5 (pg. 799)

    EM-REML tends to converge VERY slowly, because the changes at every step
    are so small. For example, a model that took 3 iterations/0m3.927s to 
    converge with AI-REML took 52 iterations/0m32.803s with EM-REML. 
    Individual EM iterations are relatively fast but since you have to 
    invert the variance-covariance matrix of observations each iteration, it
    adds up quickly.

    However, EM-REML has the nice property that it monotonically converges to
    a maximum, and avoids the parameter esimtate out-of-range problems 
    occasionally found with Newton-type methods and have to be remedied as 
    a special case. 

    Since it does take so much real time to solve the REML equations with 
    EM-REML, it's better used 

    Arguments:
    mm: a MixedModel object to be maximized
    starts: Starting values for the variance components
    maxiter: The maximum number of iterations of scoring before raising 
        an error
    tol: The minimum amount of change in the proportion of variance by any of 
        the variance components to continue iterating. 
    verbose: Print likelihood, variance component, and relative variance 
        estimates at each iteration. Useful for debugging or watching the 
        direction the optimizer is taking.
    vcs_after_maxiter: Returns estimates of variance components after reaching
        the maximum number of iterations. Useful for getting starting values 
        for variance components.

    Returns: A numpy array of the variance components at the MLE
    '''
    i = 0

    if verbose:
        print 'Maximizing model by Expectation-Maximization'

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
        new_vcs = vcs + delta

        llik = restricted_loglikelihood(mm.y, V, mm.X, P, Vinv)

        if (np.abs(delta / vcs.sum()) < tol).all():
            break

        vcs = new_vcs
        if verbose:
            print i, llik, vcs

        if i > maxiter and vcs_after_maxiter:
            break

        if i > maxiter:
            raise LinAlgError('Ran out of scoring iterations')

        i += 1
    return vcs.tolist()
