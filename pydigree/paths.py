"Functions for finding paths through pedigrees and genealogies"

from pydigree.common import table


def common_ancestors(ind1, ind2):
    """
    Common ancestors of ind1 and ind2.
    
    Recursively searches ancestors for both individuals, and
    then performs a set intersection on on each set of ancestors

    :param ind1: the first individual
    :param ind2: the second individual
    :type ind1: Individual
    :type ind2: Individual

    :returns: Common ancestors
    :rtype: set
    """
    return ind1.ancestors() & ind2.ancestors()


def path_downward(start, end, path=None):
    """
    Returns a list of paths (if they exist) from an ancestor (start)
    to a descentdant (end).

    :param start: The individual at the start of the path
    :param end: The individual to be found
    :path: a path to append to
    
    :returns: A list of individuals, in path order
    """
    if path is None:
        path = []

    path = path + [start]
    if start == end:
        return [path]
    identified_paths = []
    for child in start.children:
        newpath = path_downward(child, end, path=path)
        for p in newpath:
            identified_paths.append(p)
    return identified_paths


def paths_through_ancestor(ind1, ind2, ancestor):
    """
    Finds all paths through a genealogy between ind1 and ind2
    that pass through a specific ancestor.

    :param ind1: the first endpoint
    :param ind2: the other endpoint
    :param ancestor: the ancestor to pass through
    :type ind1: Individual
    :type ind2: Individual
    :type ancestor: Individual

    :returns: identified paths
    :rtype: list of lists of Individuals 
    """
    identified_paths = []
    ups = path_downward(ancestor, ind1)
    downs = path_downward(ancestor, ind2)
    for u in ups:
        for d in downs:
            # Often you get paths that go up and then come straight
            # back down, going through ind1. We dont want those and we're
            # not going to return them in the paths.
            if u[1] == d[1]:
                continue
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
            identified_paths.append(path)
    return identified_paths


def paths(ind1, ind2):
    """
    Returns a list of all valid paths through the pedigree connecting
    individual 1 and individual 2. A valid path can only go through an
    individual once, or it's not a valid path.

    This function consists of repeated calls to paths_through_ancestor
    for each common ancestor of ind1 and ind2

    For pedigrees, where typically many kinship coefficients must be
    calculated simulateously, kinships are calculated by this function.
    See notes on pedigree.kinship for more information.

    :param ind1: the first individual
    :param ind2: the second individual
    :type ind1: Individual
    :type ind2: Individual

    :returns: identified paths
    :rtype: list of lists of Individuals
    """
    identified_paths = path_downward(ind1, ind2) + path_downward(ind2, ind1)
    common = common_ancestors(ind1, ind2)
    for ancestor in common:
        identified_paths.extend(paths_through_ancestor(ind1, ind2, ancestor))
    return identified_paths


def kinship(ind1, ind2):
    """
    Returns the Malecot kinship coefficient for ind1 and ind2, calculated
    by path counting. The kinship coefficient is the probability that a
    randomly selected pair of alleles (one from each individual) is IBD.

    This quantity is calculated by finding the common ancestors of ind1
    and ind2. For each ancestor, the paths between ind1 and ind2 are found.

    The kinship coefficient is calculated as the sum for every ancestor of
    the sum of ((1/2) ** N) * (1+ F) for each path, where N is the number of
    individuals in the path and F is the inbreeding coefficient for that
    ancestor.

    :param ind1: the first individual
    :param ind2: the second individual
    :type ind1: Individual
    :type ind2: Individual
    :returns: Malecot's coefficient of coancestry
    :rtype: float
    """
    # In paths.fraternity, it asks for kinships between parents.
    # If the individual is a founder, they won't send individual
    # objects here, but instead Nonetypes. We'll return 0 for that
    # scenario.
    if ind1 is None or ind2 is None:
        return 0

    # Kinship with yourself (or a clone of you) is 1/2
    if ind1 == ind2:
        return 0.5

    # The most common edge case is when both individuals are founders.
    # In that scenario, the inds are unrelated by definition (or else
    # they wouldn't be founders, right?). So we'll return 0 right away.
    if ind1.is_founder() and ind2.is_founder():
        return 0

    # The number of meioses is the sum of one minus the lengths of the paths
    #  connecting the two individuals. The kinship coefficient is then
    #  .5 ** len(path)
    def kin(pathlength, ancF):
        "Quickly calculate the kinship from the path length"
        return (0.5 ** pathlength) * (1 + ancF)

    partialkin = []
    partialkin.append(
        sum(kin(len(p), ind1.inbreeding()) for p in path_downward(ind1, ind2)))
    partialkin.append(
        sum(kin(len(p), ind2.inbreeding()) for p in path_downward(ind2, ind1)))
    commonancestors = common_ancestors(ind1, ind2)
    for ancestor in commonancestors:
        ancf = ancestor.inbreeding()
        k = sum(kin(len(p), ancf)
                for p in paths_through_ancestor(ind1, ind2, ancestor))
        partialkin.append(k)
    return sum(partialkin)


def fraternity(ind1, ind2):
    """
    The coefficient of fraternity is the probability that both alleles in
    a pair of individuals are IBD.
    This is equal to: k(x_m,y_m) * k(x_f,y_f) + k(x_m,y_f) * k(x_f,y_m),
    where x_m represents the mother of x, etc. and k(x,y) represents the
    Malecot kinship between x and y.

    :param ind1: the first individual
    :param ind2: the second individual
    :type ind1: Individual
    :type ind2: Individual
    :returns: coefficient of fraternity
    :rtype: float
    """
    if ind1 is None or ind2 is None:
        return 0
    if ind1.is_founder() and ind2.is_founder():
        return 0
    a = kinship(ind1.mother, ind2.mother) * kinship(ind1.father, ind2.father)
    b = kinship(ind1.mother, ind2.father) * kinship(ind1.father, ind2.mother)
    return a + b
