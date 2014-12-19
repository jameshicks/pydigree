test='test'

from itertools import izip
import numpy as np
cimport numpy as np

cpdef ibs(g1,g2):
 a,b = g1
 c,d = g2
 
 if not (a and b and c and d):
    return None
 if (a == c and b == d) or (a == d and b == c):
    return 2
 elif a == c or a == d or b == c or b == d:
    return 1
 return 0

def get_ibs_states(ind1, ind2, chromosome_index):
    genos1 = izip(*ind1.genotypes[chromosome_index])
    genos2 = izip(*ind2.genotypes[chromosome_index])
    return [ibs(x,y) for x,y in izip(genos1, genos2)]

def runs(sequence, predicate, minlength=2):
    """
    Identifies runs of values for which predicate(value) evaluates True and yields 2-tuples of the start and end (inclusive) indices
    """
    cdef int inrun = False
    cdef int start, stop
    for i,v in enumerate(sequence):
        if not inrun and predicate(v):
            inrun = True
            start = i
        elif inrun and not predicate(v):
            inrun = False
            stop = i - 1
            if stop - start >= minlength:
                yield start, stop

    if predicate(v) and inrun:
        stop = i
        if stop - start >= minlength:
            yield start, stop
    
def set_intervals_to_value(intervals, size, value):
    DTYPE = np.int
    cdef int start = 0
    cdef int stop = 0
    cdef np.ndarray array = np.zeros(size, dtype=DTYPE)
    for start, stop in intervals:
        array[start:(stop+1)] = value
    return array