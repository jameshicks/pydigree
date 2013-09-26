#!/usr/bin/env python

from math import log as ln
import numpy as np
from numpy.linalg import inv,pinv,det

def restricted_loglikelihood(y,v,x):
    vinv = inv(v)
    p = vinv - vinv * x * pinv(x.T * vinv * x) * x.T * vinv
    return - 0.5 * ( ln(det(v) + ln( x.T * vinv * x) + y.T * p * y )
