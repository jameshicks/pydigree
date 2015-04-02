from itertools import izip, combinations

import numpy as np

from pydigree.common import runs
from pydigree.misc import ibs, get_ibs_states
from pydigree.cyfuncs import set_intervals_to_value, runs_gte_uint8


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
            shares = make_intervals(sgs_unphased(ind1, ind2, chridx, seed_size=seed_size))
            shared[pair].append(list(shares))
    return shared


def sgs_unphased(ind1, ind2, chromosome_idx, seed_size=255,
                 min_length=1, length_unit='mb'):
    ''' Returns IBD states for each marker along a chromosome '''
    
    chromosome = ind1.chromosomes[chromosome_idx] 
    identical = get_ibs_states(ind1, ind2, chromosome_idx)
    nmark = chromosome.nmark()

    # First get the segments that are IBD=1
    ibd1 = list(_process_segments(identical, min_seg=seed_size, min_val=1, chromobj=chromosome))
    ibd1 = set_intervals_to_value(ibd1, nmark, 1)

    # Then the segments that are IBD=2
    ibd2 = list(_process_segments(identical, min_seg=seed_size, min_val=2, chromobj=chromosome))
    ibd2 = set_intervals_to_value(ibd2, nmark, 2)
    return np.maximum(ibd1, ibd2)


def _process_segments(identical, min_seg=100, min_val=1, chromobj=None):    
    # IBD segments are long runs of identical genotypes
    ibd = runs_gte_uint8(identical, min_val, minlength=min_seg)

    if chromobj:
        ibd = filter_segments(chromobj, ibd)
    
    # Genotype errors are things that happen. If theres a small gap between
    # two IBD segments, we'll chalk that up to a genotyping error and join
    # them together.
    ibd = join_gaps(ibd, max_gap=2)

    return ibd


def filter_segments(chromosome, intervals,  min_size=1.0, min_density=100, size_unit='mb'):
    if size_unit == 'mb':
        locations = chromosome.physical_map
        min_size *= 1e6
        min_density /= 1e6

    elif size_unit == 'cm':
        locations = chromosome.genetic_map
    else:
        raise ValueError('Invalid size unit: {}'.format(size_unit))
    
    def meets_criteria(seg):
        start, stop = seg
        nmarkers = stop - start
        size = locations[stop] - locations[start]
        density = nmarkers / float(size)
        return size >= min_size and density >= min_density
    
    return [seg for seg in intervals if meets_criteria(seg)]


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
