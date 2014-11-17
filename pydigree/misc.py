#!/usr/bin/env python


# Extra genetics functions
def ibs(g1, g2, checkmissing=True):
    """
    Returns the number of alleles identical by state between two genotypes
    Arguements: Two tuples
    Returns: an integer
    """
    a,b,c,d = g1 + g2
    # Missing genotypes are coded as 0s
    if 0 in {a,b,c,d}:
        return None
    if (a == c and b == d) or (a == d and b == c):
        return 2
    elif a == c or a == d:
        return 1
    return 0
