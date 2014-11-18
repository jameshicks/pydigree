from itertools import izip, combinations

import numpy as np

from pydigree.common import runs
from pydigree.misc import ibs, get_ibs_states

from pydigree import Population, PedigreeCollection

def sgs_pedigrees(pc, phaseknown=False):
    shared = {}
    for pedigree in pedigrees:
        shared[pedigree] = sgs_population(pedigree)
    return shared

def sgs_population(pop, seed_size=500, phaseknown=False):
    shared = {}
    for ind1, ind2 in combinations(pop, 2):
        if not (ind1.has_genotypes() and ind2.has_genotypes()):
            continue
        pair = frozenset({ind1, ind2})
        shared[pair] = []
        for chridx, chromosome in enumerate(ind1.chromosomes):
            shares = make_intervals(_sgs_unphased(ind1, ind2, chridx, seed_size=seed_size))
            shared[pair].append(list(shares))
    return shared

def _sgs_unphased(ind1, ind2, chromosome_idx, seed_size=200):
    ''' Returns IBD states for each marker along a chromosome '''
    
    identical = get_ibs_states(ind1, ind2, chromosome_idx)

    ibd_states = np.zeros(ind1.chromosomes[chromosome_idx].nmark())

    # First get the segments that are IBD=1
    ibd1 = _process_segments(identical, min_seg=seed_size, predicate=lambda x: x > 0 or x is None)
    for start, stop in ibd1:
        ibd_states[start:(stop+1)] = 1

    # Then the segments that are IBD=2
    ibd2 = _process_segments(identical, min_seg=seed_size, predicate=lambda x: x in {2, None})
    for start, stop in ibd2:
        ibd_states[start:(stop+1)] = 2
    return ibd_states

def _process_segments(identical, min_seg=100, predicate=lambda x: x):    
    # IBD segments are long runs of identical genotypes
    ibd = runs(identical, predicate, minlength=min_seg)
    
    # Genotype errors are things that happen. If theres a small gap between
    # two IBD segments, we'll chalk that up to a genotyping error and join
    # them together.
    ibd = join_gaps(ibd, max_gap=2)

    return ibd

# Support functions

def join_gaps(seq, max_gap=1):
    seq = list(seq)

    if not seq:
        return
    elif len(seq) == 1:
        yield seq[0]
        return

    iseq = iter(seq)
    # Get the first item
    prev_start, prev_stop = iseq.next()
    for start, stop in iseq:
        if start - prev_stop > max_gap:
            yield prev_start, prev_stop
            prev_start, prev_stop = start, stop
        else:
            prev_stop = stop
    yield prev_start, stop 

def make_intervals(ibdarray):
    ibdarray = ibdarray.copy()
    # Get the intervals that are IBD=2 and remove them from the array
    ibd2_tracts = [x for x in runs(ibdarray, lambda x: x==2)]
    for start, stop in ibd2_tracts:
        ibdarray[start:(stop+1)] -= 1
    # Now get the remaining IBD=1 tracts and remove them from the array
    ibd1_tracts = [x for x in runs(ibdarray, lambda x: x > 0)]
    for start, stop in ibd1_tracts:
        ibdarray[start:(stop+1)] -= 1

    return ibd1_tracts + ibd2_tracts

def intervals_to_array(intervals, nmark):
    array = np.zeros(nmark)
    for start, stop in intervals:
        array[start:(stop+1)] += 1
    return array
