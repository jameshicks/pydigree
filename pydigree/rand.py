"Convenience functions for randomness"

import numpy as np

def choice(seq):
    """
    Randomly chooses an item from a sequence. 
    Probabilities are uniform.

    :param seq: choices
    :returns: randomly chosen item
    """
    l = len(seq)
    idx = np.random.randint(0, l)
    return seq[idx]

def sample_with_replacement(seq, n):
    """
    Choose a sample of n items with replacement from a sequence c

    :param seq: sequence to choose from
    :param n: sample size
    :returns: randomly chosen items
    :rtype: list
    """
    l = len(seq)
    return [seq[np.random.randint(0,l)] for x in range(n)]