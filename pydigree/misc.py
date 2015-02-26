from itertools import izip

import numpy as np
from pydigree.cyfuncs import ibs


# Extra genetics functions
def py_ibs(g1, g2, checkmissing=True):
    """
    Returns the number of alleles identical by state between two genotypes
    Arguements: Two tuples
    Returns: an integer
    """
    a, b = g1
    c, d = g2

    # Missing genotypes are coded as 0s
    if not (a and b and c and d):
        return None
    if (a == c and b == d) or (a == d and b == c):
        return 2
    elif a == c or a == d or b == c or b == d:
        return 1
    return 0


def get_ibs_states(ind1, ind2, chromosome_index, missingval=64):
    a, b = ind1.genotypes[chromosome_index]
    c, d = ind2.genotypes[chromosome_index]

    a_eq_c = a == c
    a_eq_d = a == d
    b_eq_c = b == c
    b_eq_d = d == c

    ibs_states = np.zeros(a.shape[0], dtype=np.uint8)
    
    # Catch which genotypes are missing so we can mark them later
    missing = (a == 0) | (b == 0) | (c == 0) | (d == 0) 
    # Any cross-genotype sharing is sufficient to be IBS=1 
    ibs1 = a_eq_c | a_eq_d | b_eq_c | b_eq_d
    # Both alleles are IBS for IBS=2. 
    ibs2 = (a_eq_c & b_eq_d) | (a_eq_d & b_eq_c) 

    ibs_states = np.zeros(a.shape[0], dtype=np.uint8)
    ibs_states[missing] = missingval
    ibs_states[ibs1] = 1
    ibs_states[ibs2] = 2 
    
    return ibs_states


def py_get_ibs_states(ind1, ind2, chromosome_index):
    genos1 = izip(*ind1.genotypes[chromosome_index])
    genos2 = izip(*ind2.genotypes[chromosome_index])
    return [ibs(x,y) for x,y in izip(genos1, genos2)]
