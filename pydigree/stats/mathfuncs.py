import numpy as np
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