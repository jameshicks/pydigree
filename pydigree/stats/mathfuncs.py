"Misc math functions"

import itertools

import numpy as np
from numpy.linalg import LinAlgError


def is_positive_definite(mat):
    """ 
    Evaluates if a matrix is positive definite (all eigvals > 0) 

    :param mat: Matrix to test
    :type mat: matrix

    :returns: positive-definiteness
    :rtype: bool
    """
    return all(np.linalg.eigvals(mat) > 0)


def grid(func, nargs, low, high, ntests=10, predicate=None):
    '''
    Evaluates a function over a range of argument values.
    
    This can be time consuming, especially if the function to be evaluated is
    particularly intensive: for m tests over n arguments, the function will be 
    evaluated m**n times

    :param func: The function to be grid searched
    :param low: The lowest value to test
    :param high: The highest value to test
    :param ntests: Number of argument values to test between low and high
    :param predicate: a function that returns True if the configuration of 
                      arguments should be evaluated.
    
    :type func: callable
    :type predicate: callable
    '''

    if predicate is None:
        predicate = lambda *x: True

    testvals = [np.linspace(low, high, ntests) for x in range(nargs)]

    for testargs in itertools.product(*testvals):

        if not predicate(*testargs):
            continue

        try:
            val = func(*testargs)
        except LinAlgError:
            val = np.nan

        yield testargs, val
