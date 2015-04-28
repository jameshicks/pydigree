from itertools import izip, combinations

import numpy as np

from pydigree.common import runs
from pydigree.misc import ibs, get_ibs_states
from pydigree.cyfuncs import set_intervals_to_value, runs_gte_uint8


from pydigree import Population, PedigreeCollection

class Segment(object):
    def __init__(self, ind1, ind2, chromobj, startidx, stopidx):
        self.ind1 = ind1
        self.ind2 = ind2
        self.chromosome = chromobj
        self.start = startidx
        self.stop = stopidx

    @property
    def physical_location(self):
        return self.chromosome.physical_map[self.start], self.chromosome.physical_map[self.stop]
    
    @property
    def genetic_location(self):
        return self.chromosome.genetic_map[self.start], self.chromosome.genetic_map[self.stop]
    
    @property
    def physical_size(self):
        start, stop = self.physical_location
        return stop - start
    
    @property
    def genetic_size(self):
        start, stop = self.genetic_location
        return stop - start

    @property
    def nmark(self):
        return self.stop - self.start

def sgs_pedigrees(pc, phaseknown=False):
    shared = {}
    for pedigree in pedigrees:
        shared[pedigree] = sgs_population(pedigree)
    return shared

def sgs_population(pop, seed_size=500, phaseknown=False, min_length=1, size_unit='mb', min_density=100):
    shared = {}
    for ind1, ind2 in combinations(pop.individuals, 2):
        if not (ind1.has_genotypes() and ind2.has_genotypes()):
            continue
        pair = frozenset({ind1, ind2})
        shared[pair] = []
        for chridx, chromosome in enumerate(ind1.chromosomes):
            shares = sgs_unphased(ind1, ind2, chridx, seed_size=seed_size, 
                                  min_length=min_length, size_unit=size_unit, min_density=min_density)
            shared[pair].append(shares)
    return shared


def sgs_unphased(ind1, ind2, chromosome_idx, seed_size=255,
                 min_length=1, size_unit='mb', min_density=100, array=False):
    ''' Returns IBD states for each marker along a chromosome '''
    
    chromosome = ind1.chromosomes[chromosome_idx] 
    identical = get_ibs_states(ind1, ind2, chromosome_idx)
    nmark = chromosome.nmark()

    # First get the segments that are IBD=1
    ibd1_segs = list(_process_segments(identical, min_seg=seed_size,
                                  min_val=1, chromobj=chromosome, 
                                  min_length=min_length, size_unit=size_unit, min_density=min_density))
    ibd1 = set_intervals_to_value(ibd1_segs, nmark, 1)

    # Then the segments that are IBD=2
    ibd2_segs = list(_process_segments(identical, min_seg=seed_size,
                                  min_val=2, chromobj=chromosome, 
                                  min_length=min_length, size_unit=size_unit, min_density=min_density))
    ibd2 = set_intervals_to_value(ibd2_segs, nmark, 2)
    ibd = np.maximum(ibd1, ibd2)
    if array:
        return ibd
    
    segs = [Segment(ind1, ind2, chromosome, start, stop) for start, stop in make_intervals(ibd)]
    return segs

def _process_segments(identical, min_seg=100, min_val=1, chromobj=None, 
                      min_density=100, size_unit='mb', min_length=1):    
    # IBD segments are long runs of identical genotypes
    ibd = runs_gte_uint8(identical, min_val, minlength=min_seg)

    if chromobj:
        ibd = filter_segments(chromobj, ibd, min_length=min_length,
                              size_unit=size_unit, min_density=min_density)
    
    # Genotype errors are things that happen. If theres a small gap between
    # two IBD segments, we'll chalk that up to a genotyping error and join
    # them together.
    ibd = join_gaps(ibd, max_gap=2)

    return ibd


def filter_segments(chromosome, intervals,  min_length=1.0, min_density=100, size_unit='mb'):
    size_unit = size_unit.lower()
    if size_unit == 'mb':
        locations = chromosome.physical_map
        min_length *= 1e6
        min_density /= 1e6
    elif size_unit == 'kb':
        locations = chromosome.physical_map
        min_length *= 1000
        min_density /= 1e3
    elif size_unit == 'cm':
        locations = chromosome.genetic_map
    else:
        raise ValueError('Invalid size unit: {}'.format(size_unit))
    
    def meets_criteria(seg):
        start, stop = seg
        nmarkers = stop - start
        size = locations[stop] - locations[start]
        density = nmarkers / float(size)
        return size >= min_length and density >= min_density
    
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
