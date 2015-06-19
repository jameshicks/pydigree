from itertools import izip, combinations, product

import numpy as np

from pydigree.common import runs
from pydigree.ibs import ibs, get_ibs_states
from pydigree.cyfuncs import set_intervals_to_value, runs_gte_uint8


from pydigree import Population, PedigreeCollection


class Segment(object):
    __slots__ = ['ind1', 'ind2', 'chromosome', 'start', 'stop', '_chridx']

    def __init__(self, ind1, ind2, chromobj, startidx, stopidx):
        self.ind1 = ind1
        self.ind2 = ind2
        self.chromosome = chromobj
        self._chridx = ind1.chromosomes.index(chromobj)
        self.start = startidx
        self.stop = stopidx

    @property
    def physical_location(self):
        ''' Locations of the start and end of the segment in base pairs '''
        return self.chromosome.physical_map[self.start], self.chromosome.physical_map[self.stop]

    @property
    def genetic_location(self):
        ''' Locations of the start and end of the segment in centimorgans '''
        return self.chromosome.genetic_map[self.start], self.chromosome.genetic_map[self.stop]

    @property
    def physical_size(self):
        ''' The size of the segmend in base pairs '''
        start, stop = self.physical_location
        return stop - start

    @property
    def genetic_size(self):
        ''' The size of the segment in centimorgans '''
        start, stop = self.genetic_location
        return stop - start

    @property
    def nmark(self):
        ''' The number of markers in the segment '''
        return self.stop - self.start

    @property
    def missing(self):
        chridx = self._chridx
        miss1 = ind1.genotypes[chridx][
            0].missing | ind1.genotypes[chridx][1].missing
        miss2 = ind2.genotypes[chridx][
            0].missing | ind2.genotypes[chridx][1].missing
        miss = miss1 | miss2
        return miss[self.start:(self.stop + 1)]

    @property
    def missing_rate(self):
        ''' The number of missing genotypes in the segment '''
        return self.missing.sum() / float(self.nmark)


def sgs_pedigrees(pc, phaseknown=False):
    ''' Performs within-pedigree SGS for each pedigree in a pedigree collection '''
    shared = {}
    for pedigree in pedigrees:
        shared[pedigree] = sgs_population(pedigree)
    return shared


def sgs_population(pop, seed_size=500, phaseknown=False, min_length=1,
                   size_unit='mb', min_density=100, maxmiss=0.25):
    ''' Performs SGS between all individuals in a population or pedigree '''
    shared = {}
    for ind1, ind2 in combinations(pop.individuals, 2):
        if not (ind1.has_genotypes() and ind2.has_genotypes()):
            continue
        pair = frozenset({ind1, ind2})
        shared[pair] = []
        for chridx, chromosome in enumerate(ind1.chromosomes):
            shares = sgs_unphased(ind1, ind2, chridx, seed_size=seed_size,
                                  min_length=min_length, size_unit=size_unit, min_density=min_density, maxmiss=0.25)
            shared[pair].append(shares)
    return shared

def sgs_autozygous(ind, chromosome_idx, seed_size=500, phaseknown=False, min_length=1,
                   size_unit='mb', min_density=100, maxmiss=0.25):
    chromosome = ind.chromosomes[chromosome_idx] 
    hapa, hapb = ind.genotypes[chromosome_idx]
    homozygous = (hapa == hapb).astype(np.uint8)
    autozygous_segs = list(_process_segments(homozygous, min_seg=seed_size,
                                       min_val=1, chromobj=chromosome,
                                       min_length=min_length, size_unit=size_unit,
                                       min_density=min_density, maxmiss=0.25))
    return [Segment(ind, ind, chromosome, start, stop) for start, stop in autozygous_segs]

def sgs_unphased(ind1, ind2, chromosome_idx, seed_size=255,
                 min_length=1, size_unit='mb', min_density=100,
                 maxmiss=0.25, array=False):
    ''' Returns IBD states for each marker along a chromosome '''

    chromosome = ind1.chromosomes[chromosome_idx]
    identical = get_ibs_states(ind1, ind2, chromosome_idx)
    nmark = chromosome.nmark()

    # First get the segments that are IBD=1
    ibd1_segs = list(_process_segments(identical, min_seg=seed_size,
                                       min_val=1, chromobj=chromosome,
                                       min_length=min_length, size_unit=size_unit,
                                       min_density=min_density, maxmiss=0.25))
    ibd1 = set_intervals_to_value(ibd1_segs, nmark, 1)

    # Then the segments that are IBD=2
    ibd2_segs = list(_process_segments(identical, min_seg=seed_size,
                                       min_val=2, chromobj=chromosome,
                                       min_length=min_length, size_unit=size_unit,
                                       min_density=min_density, maxmiss=0.25))
    ibd2 = set_intervals_to_value(ibd2_segs, nmark, 2)
    ibd = np.maximum(ibd1, ibd2, dtype=np.uint8)
    if array:
        return ibd

    segs = make_intervals(ibd)
    segs = [Segment(ind1, ind2, chromosome, start, stop)
            for start, stop in segs]
    return segs


def _process_segments(identical, min_seg=100, min_val=1, chromobj=None,
                      min_density=100, size_unit='mb', min_length=1, maxmiss=0.25):
    # IBD segments are long runs of identical genotypes
    ibd = runs_gte_uint8(identical, min_val, minlength=min_seg)

    if not ibd:
        return ibd

    # Genotype errors are things that happen. If theres a small gap between
    # two IBD segments, we'll chalk that up to a genotyping error and join
    # them together.
    ibd = join_gaps(ibd, max_gap=2)

    if chromobj:
        ibd = filter_segments(chromobj, ibd, identical, min_length=min_length,
                              size_unit=size_unit, min_density=min_density, maxmiss=0.25)

    return ibd


def filter_segments(chromosome, intervals, identical, min_length=1.0,
                    min_density=100, size_unit='mb', maxmiss=0.25):
    ''' Perform quality control filtering on SGS results '''
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

    missing = identical > 2

    def meets_criteria(seg):
        start, stop = seg
        nmarkers = stop - start
        size = locations[stop] - locations[start]
        density = nmarkers / float(size)
        missrate = missing[start:stop].sum() / float(nmarkers)
        return size >= min_length and density >= min_density and missrate <= maxmiss

    return [seg for seg in intervals if meets_criteria(seg)]


def ibd_state(shared, ind1, ind2, locus):
    ''' Gets the IBD state between two individuals at a locus '''
    pair = frozenset([ind1, ind2])
    if pair not in shared:
        return 0
    ibd = sum(1 for start, stop in shared[pair] if start <= locus <= (stop+1))
    return ibd


def ibd_matrix(shared, individuals, locus):
    nind = len(individuals)
    mat = []
    for ind1 in individuals:
        row = [ibd_state(shared, ind1, ind2, locus) for ind2 in individuals]
        mat.append(row)
    mat = np.matrix(mat, dtype=np.uint8)
    return mat

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
    ibdarray = np.array(ibdarray, dtype=np.uint8)
    ibdarray = ibdarray.copy()
    # Get the intervals that are IBD=2 and remove them from the array
    ibd2_tracts = [x for x in runs_gte_uint8(ibdarray, 2)]
    for start, stop in ibd2_tracts:
        ibdarray[start:(stop + 1)] -= 1
    # Now get the remaining IBD=1 tracts and remove them from the array
    ibd1_tracts = [x for x in runs_gte_uint8(ibdarray, 1)]
    for start, stop in ibd1_tracts:
        ibdarray[start:(stop + 1)] -= 1

    return ibd1_tracts + ibd2_tracts


def intervals_to_array(intervals, nmark):
    array = np.zeros(nmark)
    for start, stop in intervals:
        array[start:(stop + 1)] += 1
    return array
