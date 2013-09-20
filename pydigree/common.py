#!/usr/bin/env python

def count(val,iter):
    """ Counts how many times a value (val) occurs in an iterable """
    return sum(1 for x in (y for y in iter if y is not None) if val == x)
def table(seq):
    """
    For each unique value in seq, runs count() on it. Returns a dictionary in
    the form of {value1: count, value2: count}.
    """
    seq = [x for x in seq]
    keys = set(seq)
    return dict([(k,seq.count(k)) for k in keys])
def flatten(x):
    """
    Recursively flattens lists.
    From: http://stackoverflow.com/questions/5409224/python-recursively-flatten-a-list
    """
    try:
        it = iter(x)
    except TypeError:
        yield x
    else:
        for i in it:
            for j in flatten(i):
                yield j

# Common stats functions
def cumsum(iter):
    """
    Cumulative sum:
    >>> pydigree.cumsum([0,1,2,3,4])
    [0, 1, 3, 6, 10]
    """
    value = 0
    g = [None] * len(iter)
    for idx,x in enumerate(iter):
        value += x
        g[idx]=value
    return g

def invcumsum(iter):
    """
    Inverse of cumsum. iter can't be a generator currently

    >>> a = range(5)
    >>> pydigree.invcumsum(pydigree.cumsum(a)) == a
    True
    
    """
    return [x - iter[i-1] if i > 0 else 0 for i,x in enumerate(iter)]

