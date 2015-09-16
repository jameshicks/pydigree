from itertools import izip, combinations, product, chain, imap
from functools import partial
import multiprocessing

import numpy as np
from scipy.sparse import lil_matrix

from pydigree.common import runs
from pydigree.ibs import ibs, get_ibs_states
from pydigree.cyfuncs import set_intervals_to_value, runs_gte_uint8

from pydigree import Population, PedigreeCollection
from pydigree import Individual, ChromosomeTemplate


class SGSAnalysis(object):

    def __init__(self, pairs=None):
        if pairs:
            self.pairs = pairs
        else:
            self.pairs = {}

    def __getitem__(self, item):
        if type(item) is not frozenset:
            item = frozenset(item)
        return self.pairs[item]

    def __setitem__(self, k, v):
        self.pairs[k] = v

    def __contains__(self, item):
        item = frozenset(item)
        return item in self.pairs

    def merge(self, other):
        ''' Merge two SGSAnalysis objects '''
        self.pairs.update(other.pairs)

    def ibd_state(self, ind1, ind2, locus, location_type='index', onlywithin=False):
        ''' 
        Gets the IBD state between two individuals at a locus

        locus: a 2-tuple in the form (chromosome, position)
        location_type: the type of location units specified. Valid entries are 
        'index' (the index of the markers), 'physical' (the positions in bp),
        and 'genetic' (location in cM)

        Returns 0 in the case of ind1 and ind2 not having any identified
        segments (or not being in the analysis at all)
        '''
        pair = frozenset([ind1, ind2])
        if onlywithin and (ind1.population.label != ind2.population.label):
            return 0
        if pair not in self.pairs:
            return 0
        else:
            return self.pairs[pair].ibd_state(locus, location_type)

    def ibd_matrix(self, individuals, locus, location_type='index', onlywithin=False):
        '''
        Creates an IBD matrix for a set of individuals at a locus

        individuals: A list of individuals to form the matrix for
        locus: a 2-tuple in the form (chromosome, position)
        location_type: the type of location units specified. Valid entries are 
        'index' (the index of the markers), 'physical' (the positions in bp),
        and 'genetic' (location in cM)
        '''
        nind = len(individuals)
        mat = []
        for ind1 in individuals:
            row = [self.ibd_state(ind1, ind2, locus, location_type, onlywithin=onlywithin)
                   for ind2 in individuals]
            mat.append(row)
        mat = np.matrix(mat) / 2.0
        mat = lil_matrix(mat)
        mat.setdiag(1)

        return mat

    def chromwide_ibd(self, chromidx, size=None, onlywithin=False, onlybetween=False):
        if onlybetween and onlywithin:
            raise ValueError(
                'Cannot set onlywithin and onlywithin simultaneously')
        if size is None:
            size = max([x.stop for x in self.segments])

        arry = np.zeros(size, dtype=np.uint)
        for seg in self.segments:
            if onlywithin and seg.ind1.full_label[0] != seg.ind2.full_label[0]:
                continue
            if onlybetween and seg.ind1.full_label[0] == seg.ind2.full_label[0]:
                continue
            arry[seg.start:seg.stop] += 1
        return arry

    @property
    def segments(self):
        ''' Returns an iterable of all the segments in the Analysis '''
        return chain.from_iterable([x.segments for x in self.pairs.values()])

    @property
    def individuals(self):
        """ Returns a set of all individuals present in the analysis """
        inds = reduce(frozenset.union, self.pairs.keys())
        return inds

    def update_segment_references(self, pedigrees):
        '''
        Replace individual labels in the Segment objects of SGSAnalyses read 
        from text with references to the actual individual object 
        '''
        pedindlabs = frozenset([x.full_label for x in pedigrees.individuals])
        sgsindlabs = frozenset(self.individuals)

        if sgsindlabs - pedindlabs:
            raise ValueError('Not all individuals present in SGS are present'
                             'in the given pedigrees')

        chroms = {chrom.label: chrom for chrom in pedigrees.chromosomes}
        chromindices = {chrom: i for i, chrom in enumerate(pedigrees.chromosomes)}
        for segment in self.segments:
            try:
                segment.ind1 = pedigrees[segment.ind1]
                segment.ind2 = pedigrees[segment.ind2]

                if type(segment.chromosome) is str:
                    segment.chromosome = chroms[segment.chromosome]
                    segment._chridx = chromindices[segment.chromosome]

                pstart, pstop = segment.physical_location
                segment.start = segment.chromosome.physical_map.index(pstart)
                segment.stop = segment.chromosome.physical_map.index(pstop)

            except KeyError:
                if segment.ind1 not in pedindlabs:
                    raise ValueError('{} not in pedigree'.format(segment.ind1))
                elif segment.ind2 not in pedindlabs:
                    raise ValueError('{} not in pedigree'.format(segment.ind2))
                else:
                    raise Exception('Unknown error')
            except ValueError:
                raise ValueError('Positions not in chromosome data')

        def pairlookup(pair):
            newpair = frozenset({pedigrees[ind] for ind in pair})
            return newpair
        self.pairs = {pairlookup(k): v for k, v in self.pairs.items()}


class SGS(object):

    def __init__(self, ind1, ind2, segments=None):
        self.ind1 = ind1
        self.ind2 = ind2
        self.segments = segments if segments is not None else []

    def __getitem__(self, idx):
        return self.segments[idx]

    def __iter__(self):
        for x in self.segments:
            yield x

    def append(self, value):
        self.segments.append(value)

    def extend(self, value):
        self.segments.extend(value)

    def ibd_state(self, locus, location_type='index'):
        '''
        locus: a 2-tuple in the form (chromosome, position)
        location_type: the type of location units specified. Valid entries are 
        'index' (the index of the markers), 'physical' (the positions in bp),
        and 'genetic' (location in cM)
        '''
        chrom, pos = locus

        if location_type == 'index':
            ibd = sum(1 for x in self.segments if x._chridx == chrom and
                      x.start <= pos <= (x.stop+1))
        elif location_type == 'physical':
            segs = [(x._chridx, x.physical_start, x.physical_stop)
                    for x in self.segments]
            ibd = sum(1 for segchrom, start, stop in segs
                      if segchrom == chrom and (start <= pos <= stop))
        elif location_type == 'genetic':
            segs = [(x._chridx, x.genetic_start, x.genetic_stop)
                    for x in self.segments]
            ibd = sum(1 for segchrom, start, stop in segs
                      if segchrom == chrom and (start <= pos <= stop))
        return ibd


class Segment(object):
    __slots__ = ['ind1', 'ind2', 'chromosome', 'start', 'stop',
                 '_chridx', 'physical_location', 'genetic_location']

    def __init__(self, ind1, ind2, chromosome, startidx, stopidx,
                 physical_location=None, genetic_location=None):
        self.ind1 = ind1
        self.ind2 = ind2
        self.chromosome = chromosome

        if isinstance(ind1, Individual):
            self._chridx = ind1.chromosomes.index(chromosome)

        self.start = startidx
        self.stop = stopidx

        # Get positions from chomosome if possible
        if isinstance(self.chromosome, ChromosomeTemplate):
            physical_start = self.chromosome.physical_map[self.start]
            physical_stop = self.chromosome.physical_map[self.stop]
            self.physical_location = physical_start, physical_stop

        # If not that, try directly
        elif physical_location is not None:
            pstart, pstop = physical_location
            physical_start = int(pstart)
            physical_stop = int(pstop)
            self.physical_location = physical_start, physical_stop

        # Give up and go with None
        else:
            self.physical_location = None

        # Do the same thing with genetic positions
        if isinstance(self.chromosome, ChromosomeTemplate):
            genetic_start = self.chromosome.genetic_map[self.start]
            genetic_stop = self.chromosome.genetic_map[self.stop]
            self.genetic_location = genetic_start, genetic_stop

        elif genetic_location is not None:
            gstart, gstop = genetic_location
            genetic_start = float(gstart)
            genetic_stop = float(gstop)
            self.genetic_location = genetic_start, genetic_stop

        else:
            self.genetic_location = None

    @property
    def marker_labels(self):
        startlab = self.chromosome.labels[self.start]
        stoplab = self.chromosome.labels[self.stop]
        return startlab, stoplab

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
        miss1 = (ind1.genotypes[chridx][0].missing |
                 ind1.genotypes[chridx][1].missing)
        miss2 = (ind2.genotypes[chridx][0].missing |
                 ind2.genotypes[chridx][1].missing)
        miss = miss1 | miss2
        return miss[self.start:(self.stop + 1)]

    @property
    def missing_rate(self):
        ''' The number of missing genotypes in the segment '''
        return self.missing.sum() / float(self.nmark)


def sgs_pedigrees(pc, phaseknown=False):
    '''
    Performs within-pedigree SGS for each pedigree in a pedigree collection
    '''
    shared = SGSAnalysis()
    for pedigree in pedigrees:
        shared[pedigree] = sgs_population(pedigree)
    return shared


def _pair_sgs(pair, seed_size=500, phaseknown=False,
              min_length=1, size_unit='mb',
              min_density=100, maxmiss=0.25,
              onlywithin=False):
    ind1, ind2 = pair
    if not (ind1.has_genotypes() and ind2.has_genotypes()):
        return None

    if onlywithin and (ind1.full_label[0] != ind2.full_label[0]):
        return None

    if ind1 == ind2:
        results = SGS(ind1, ind2)
        for chridx, chromosome in enumerate(ind1.chromosomes):
            shares = sgs_autozygous(ind1, chridx,
                                    seed_size=seed_size,
                                    min_length=min_length,
                                    size_unit=size_unit,
                                    min_density=min_density,
                                    maxmiss=maxmiss)
            results.extend(shares)

    else:
        results = SGS(ind1, ind2)
        for chridx, chromosome in enumerate(ind1.chromosomes):
            shares = sgs_unphased(ind1, ind2, chridx,
                                  seed_size=seed_size,
                                  min_length=min_length,
                                  size_unit=size_unit,
                                  min_density=min_density,
                                  maxmiss=maxmiss)
            results.extend(shares)
    return results


def sgs_population(pop, seed_size=500, phaseknown=False,
                   min_length=1, size_unit='mb',
                   min_density=100, maxmiss=0.25,
                   onlywithin=False, njobs=1):
    ''' Performs SGS between all individuals in a population or pedigree '''
    pair_sgs = partial(_pair_sgs, seed_size=seed_size,
                         min_length=min_length,
                         size_unit=size_unit,
                         min_density=min_density,
                         maxmiss=maxmiss)
    shared = SGSAnalysis()
    pairs = [x for x in combinations(pop.individuals, 2)]

    njobs = int(njobs)
    if njobs == 1:
        res = imap(pair_sgs, pairs)
    elif njobs > 1:
        pool = multiprocessing.Pool(processes=njobs)
        res = pool.imap(pair_sgs, pairs, chunksize=1000)
    else:
        raise ValueError('Bad value for njobs: {}'.format(njobs))

    for result in res:
        pair = frozenset([result.ind1, result.ind2])
        shared[pair] = result
    return shared


def sgs_autozygous(ind, chromosome_idx, seed_size=500,
                   phaseknown=False, min_length=1, size_unit='mb',
                   min_density=100, maxmiss=0.25):
    chromosome = ind.chromosomes[chromosome_idx]
    hapa, hapb = ind.genotypes[chromosome_idx]
    homozygous = (hapa == hapb).astype(np.uint8)
    autozygous_segs = list(_process_segments(homozygous,
                                             min_seg=seed_size,
                                             min_val=1,
                                             chromobj=chromosome,
                                             min_length=min_length,
                                             size_unit=size_unit,
                                             min_density=min_density,
                                             maxmiss=maxmiss))
    return [Segment(ind, ind, chromosome, start, stop)
            for start, stop in autozygous_segs]


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
                                       min_length=min_length,
                                       size_unit=size_unit,
                                       min_density=min_density,
                                       maxmiss=maxmiss))
    ibd1 = set_intervals_to_value(ibd1_segs, nmark, 1)

    # Then the segments that are IBD=2
    ibd2_segs = list(_process_segments(identical, min_seg=seed_size,
                                       min_val=2, chromobj=chromosome,
                                       min_length=min_length,
                                       size_unit=size_unit,
                                       min_density=min_density,
                                       maxmiss=maxmiss))
    ibd2 = set_intervals_to_value(ibd2_segs, nmark, 2)
    ibd = np.maximum(ibd1, ibd2, dtype=np.uint8)
    if array:
        return ibd

    segs = make_intervals(ibd)
    segs = [Segment(ind1, ind2, chromosome, start, stop)
            for start, stop in segs]
    return segs


def _process_segments(identical, min_seg=100, min_val=1, chromobj=None,
                      min_density=100, size_unit='mb',
                      min_length=1, maxmiss=0.25):
    # IBD segments are long runs of identical genotypes
    ibd = runs_gte_uint8(identical, min_val, minlength=min_seg)

    if not ibd:
        return ibd

    # Genotype errors are things that happen. If theres a small gap between
    # two IBD segments, we'll chalk that up to a genotyping error and join
    # them together.
    ibd = join_gaps(ibd, max_gap=2)

    if chromobj:
        ibd = filter_segments(chromobj, ibd, identical,
                              min_length=min_length,
                              size_unit=size_unit,
                              min_density=min_density,
                              maxmiss=maxmiss)

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
        return (size >= min_length and
                density >= min_density and
                missrate <= maxmiss)

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
