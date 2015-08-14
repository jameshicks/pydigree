#!/usr/bin/env python

from itertools import izip_longest
from operator import mul as multiply
from math import log
from pydigree.cyfuncs import runs, runs_gte


def count(val, iter):
    """ Counts how many times a value (val) occurs in an iterable """
    return sum(1 for x in (y for y in iter if y is not None) if val == x)


def table(seq):
    """
    For each unique value in seq, runs count() on it. Returns a dictionary in
    the form of {value1: count, value2: count}.
    """
    seq = [x for x in seq]
    keys = set(seq)
    return dict([(k, seq.count(k)) for k in keys])


def flatten(x):
    """
    Recursively flattens lists.
    From:
    stackoverflow.com/questions/5409224/python-recursively-flatten-a-list
    """
    try:
        it = iter(x)
    except TypeError:
        yield x
    else:
        for i in it:
            for j in flatten(i):
                yield j


def grouper(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks
    Iterable: input iterable
    n: size of the groups (an integer)
    fillvalue: Value to pad the last group with if len(iterable) % n != 0
    
    Returns: a generator
    """
    
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def invert_dict(d):
    """
    Makes the keys the values and the values the keys
    WARNING: No guarantee of dictionary structure if mappings are not unique
    
    Arguments
    ------
    d: Input dictionary

    Returns: a dictionary
    """
    return dict((y, x) for x, y in d.items())

def merge_dicts(*args):
    ''' Merges two dictionaries into one bigger one '''
    odict = {}
    for d in args:
        odict.update(d)
    return odict

# Common stats/math functions
def log_base_change(value, old, new):
    ''' 
    Changes the base of the logarithm used on `value` from `old` to `new`

    Arguments: 
    value: the value to be converted
    old: old base (numeric)
    new: new base (numeric)

    returns: a float
    '''
    return value / log(old, new) 


def product(iter):
    """
    Reduces an iterable by multiplication. Analogous to sum, but with
    multiplication instead of addition.
    """
    # This should really be a python builtiin
    if not iter:
        # sum([]) returns 0, so this will return 1
        return 1
    return reduce(multiply, iter)


def cumsum(iter):
    """
    Cumulative sum:
    >>> pydigree.cumsum([0,1,2,3,4])
    [0, 1, 3, 6, 10]
    
    Arguments
    ------
    iter: the iterable to be cumsum'ed
    
    Returns: a list of the cumulative sums
    """
    if not iter:
        return []
    value = 0
    g = [None] * len(iter)
    for idx, x in enumerate(iter):
        value += x
        g[idx] = value
    return g
