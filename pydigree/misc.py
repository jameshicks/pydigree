#!/usr/bin/env python


# Extra genetics functions
def is_missing_genotype(g):
    return g == tuple([0, 0])


def ibs(g1, g2, checkmissing=True):
    """
    Returns the number of alleles identical by state between two genotypes
    Arguements: Two tuples
    Returns: an integer
    """
    if checkmissing:
        if is_missing_genotype(g1) or is_missing_genotype(g2):
            return None
    
    g1 = set(g1)
    g2 = set(g2)
    if g1 == g2:
        return 2
    if g1 & g2:
        return 1
    return 0
