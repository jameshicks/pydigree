#!/usr/bin/env python


### Extra genetics functions
def is_missing_genotype(g):
    return g == tuple([0, 0])


def ibs(g1, g2):
    """
    Returns the number of alleles identical by state between two genotypes
    Arguements: Two tuples
    Returns: an integer
    """
    g1, g2 = sorted(g1), sorted(g2)
    if is_missing_genotype(g1) or is_missing_genotype(g2):
        return None
    if g1 == g2:
        return 2
    g2s = set(g2)
    if g1[0] in g2s or g1[1] in g2s:
        return 1
    return 0
