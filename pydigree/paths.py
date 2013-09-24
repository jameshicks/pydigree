#!/usr/bin/env python

from common import table

### Functions for finding paths through pedigrees

def common_ancestors(ind1,ind2):
    """ Common ancestors of ind1 and ind2. A quick set operation on Individual.ancestors() """
    return ind1.ancestors() & ind2.ancestors()
def path_downward(start,end,path=[]):
    """
    Returns a list of paths (if they exist) from an ancestor 'start' to a descentdant 'end'
    """
    path = path + [start]
    if start == end: return [path]
    paths = []
    for child in start.children:
        newpath = path_downward(child,end,path=path)
        for p in newpath: paths.append(p)
    return paths
def paths_through_ancestor(ind1,ind2,ancestor):
    """
    Returns a list of paths through the pedigree between ind1 and ind2 that
    pass through a specific ancestor.
    """
    paths = []
    ups = path_downward(ancestor,ind1)
    downs = path_downward(ancestor,ind2)
    for u in ups:
        for d in downs:
            # Often you get paths that go up and then come straight
            # back down, going through ind1. We dont want those and we're
            # not going to return them in the paths.
            if u[1] == d[1]: continue
            # Since you're pathing downward in both partial paths through
            # the common ancestor, the path from ind1 -> ancestor needs to
            # be reversed. Then we skip the first element of the downward
            # path because its the last element of the upward path.
            path = u[::-1] + d[1:]
            # In sort of a generalization of the first check for the up
            # and straight-back-down scenario, valid paths can only go
            # through an individual ONCE. So if any path has an individual
            # more than once, we'll skip that path entirely
            pt = table(path)
            if [pt[k] for k in pt if pt[k] > 1]:
                continue
            paths.append(path)
    return paths
def paths(ind1,ind2):
    """
    Returns a list of all valid paths through the pedigree connecting
    individual 1 and individual 2. A valid path can only go through an
    individual once, or it's not a valid path.

    This function consists of repeated calls to paths_through_ancestor
    for each common ancestor of ind1 and ind2
    """
    paths = path_downward(ind1,ind2) + path_downward(ind2,ind1)
    common = common_ancestors(ind1,ind2)
    for ancestor in common:
        paths.extend(paths_through_ancestor(ind1,ind2,ancestor))
    return paths
def kinship(ind1, ind2):
    """
    Returns the Malecot kinship coefficient for ind1 and ind2, calculated by path counting.
    The kinship coefficient is the probability that a randomly selected pair of alleles
    (one from each individual) is IBD.
    
    This quantity is calculated by finding the common ancestors of ind1 and ind2.
    For each ancestor, the paths between ind1 and ind2 are found.
    
    The kinship coefficient is calculated as the sum for every ancestor of the sum
    of ((1/2) ** N) * (1+ F) for each path, where N is the number of individuals in the path
    and F is the inbreeding coefficient for that ancestor.
    """
    # The most common edge case is when both individuals are founders.
    # In that scenario, the inds are unrelated by definition (or else
    # they wouldn't be founders, right?). So we'll return 0 right away.
    if ind1.is_founder() and ind2.is_founder():
        return 0
    # The number of meioses is the sum of one minus the lengths of the paths
    #  connecting the two individuals. The kinship coefficient is then
    #  .5 ** len(path)
    def kin(pathlength,ancF): return (0.5 ** pathlength) * (1+ancF)
    partialkin = []
    partialkin.append( sum( kin(len(p), ind1.inbreeding()) for p in path_downward(ind1,ind2)) )
    partialkin.append( sum( kin(len(p), ind2.inbreeding()) for p in path_downward(ind2,ind1)) )
    commonancestors = common_ancestors(ind1,ind2)
    for ancestor in commonancestors:
        ancf = ancestor.inbreeding()
        k = sum( kin(len(p), ancf) \
                 for p in paths_through_ancestor(ind1,ind2,ancestor) )
        partialkin.append(k)
    return sum(partialkin)
