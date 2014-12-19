from itertools import izip

from pydigree.cyfuncs import ibs, get_ibs_states

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

def py_get_ibs_states(ind1, ind2, chromosome_index):
    genos1 = izip(*ind1.genotypes[chromosome_index])
    genos2 = izip(*ind2.genotypes[chromosome_index])
    return [ibs(x,y) for x,y in izip(genos1, genos2)]
