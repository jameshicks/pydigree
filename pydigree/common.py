#!/usr/bin/env python
import numpy as np
from itertools import zip_longest
from operator import mul as multiply
from functools import reduce
from math import log

from pydigree.cydigree.cyfuncs import runs, runs_gte, interleave
from pydigree.cydigree.cyfuncs import all_same_type, is_sorted


def count(val, iterable):
    """ 
    Counts how many times a value (val) occurs in an iterable, excluding `None`

    
    :param val: The value to be counted
    :param iterable: values to be counted over
    :type iterable: iterable  

    :return: the count of values
    :rtype: int
    """
    return sum(1 for x in (y for y in iterable if y is not None) if val == x)


def table(seq):
    """
    For each unique value in seq, runs count() on it. Returns a dictionary in
    the form of {value1: count, value2: count}.

    :param seq: Values to make a table of

    :rtype: dict
    """
    seq = [x for x in seq]
    keys = set(seq)
    return dict([(k, seq.count(k)) for k in keys])

def mode(seq):
    """ 
    Returns the most common value in a sequence 

    :param seq: sequence of values to evaluate
    """

    if not len(seq):
        raise IndexError('Sequence is empty')
    tab = table(seq)
    return sorted(tab.items(), key=lambda x: x[1], reverse=True)[0][0]

def random_choice(iterable):
    ''' Randomly chooses an item from an iterable '''
    itersize = len(iterable)
    randidx = np.random.randint(0, itersize)
    return iterable[randidx]

def flatten(x):
    """
    Recursively flattens lists.
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
    :param iterable: values to evaluate
    :param n: size of the groups
    :param fillvalue: Value to pad the last group with if len(iterable) % n != 0
    
    :type n: integer
    :rtype: a generator
    """
    
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def invert_dict(d):
    """
    Makes the keys the values and the values the keys
    
    .. warning:: 
        No guarantee of dictionary structure if mappings are not unique
    
    :param d: dictionary to be inverted

    :rtype: dict
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

    :returns: log of value in new base  
    :rtype: float
    '''
    return value / log(old, new) 


def product(iter):
    """
    Reduces an iterable by multiplication. Analogous to sum, but with
    multiplication instead of addition.

    :returns: overall product
    :rtype: numeric
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
    

    :param iter: the iterable to be cumsum'ed
    
    Returns: cumulative sums
    :rtype: integer
    """
    if not iter:
        return []
    value = 0
    g = [None] * len(iter)
    for idx, x in enumerate(iter):
        value += x
        g[idx] = value
    return g



