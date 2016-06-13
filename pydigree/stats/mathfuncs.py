import itertools

import numpy as np
from numpy.linalg import LinAlgError
from scipy.optimize import approx_fprime

def is_positive_definite(mat):
    "Returns true if a matrix is positive definite (all eigvals > 0)"
    return all(np.linalg.eigvals(mat) > 0)


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