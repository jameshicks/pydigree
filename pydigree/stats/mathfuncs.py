import itertools

import numpy as np
from numpy.linalg import LinAlgError
from scipy.optimize import approx_fprime

def is_positive_definite(mat):
    "Returns true if a matrix is positive definite (all eigvals > 0)"
    return all(np.linalg.eigvals(mat) > 0)

def approx_hessian(x0, func, epsilon=1.e-5, linear_approx=False, *args):
    """
    A numerical approximation to the Hessian matrix of cost function at
    location x0 (hopefully, the minimum)

    Adapted from: https://gist.github.com/jgomezdans/3144636
    """
    # ``calculate_cost_function`` is the cost function implementation
    # The next line calculates an approximation to the first
    # derivative
    f1 = approx_fprime(x0, func, epsilon, *args)

    # This is a linear approximation. Obviously much more efficient
    # if cost function is linear
    if linear_approx:
        f1 = np.matrix(f1)
        return f1.transpose() * f1
    # Allocate space for the hessian
    n = x0.shape[0]
    hessian = np.zeros((n, n))
    # The next loop fill in the matrix
    xx = x0
    for j in xrange(n):
        xx0 = xx[j]  # Store old value
        xx[j] = xx0 + epsilon  # Perturb with finite difference
        # Recalculate the partial derivatives for this new point
        f2 = approx_fprime(x0, func, epsilon, *args)
        hessian[:, j] = (f2 - f1)/epsilon  # scale...
        xx[j] = xx0  # Restore initial value of x0
    return hessian


def grid(func, nargs, low, high, ntests=10, predicate=None, catch=None):
    '''
    Evaluates a function over a range of argument values.
    
    This can be time consuming, especially if the function to be evaluated is
    particularly intensive: for m tests over n arguments, the function will be 
    evaluated m**n times

    func: The function to be grid searched
    low: The lowest value to test
    high: The highest value to test
    ntests: Number of argument values to test between low and high
    predicate: a function that returns True if the configuration of 
               arguments should be evaluated.
    '''
    if predicate is None:
        predicate = lambda *x: True

    testvals = [np.linspace(low, high, ntests) for x in xrange(nargs)]

    for testargs in itertools.product(*testvals):

        if not predicate(*testargs):
            continue

        try:
            val = func(*testargs)
        except LinAlgError:
            val = np.nan

        yield testargs, val