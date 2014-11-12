from itertools import izip
from pydigree.common import runs
from pydigree import Population, PedigreeCollection

def sgs_pedigrees(pc, phaseknown=False):
    shared = {}
    for pedigree in pedigrees:
        shared[pedigree] = sgs_population(pedigree)
    return shared

def sgs_population(pop, phaseknown=False):
    shared = {}
    for ind1, ind2 in itertools.combinations(pop, 2):
        for chridx,chromosome in enumerate(ind1.population.chromosomes):
            shared[pair] = _sgs_unphased(ind1, ind2, chridx)
                

def _sgs_unphased(ind1, ind2, chromosome_idx):
    ''' Returns IBD states for each marker along a chromosome '''
    min_seg = 100
    
    genos1 = izip(ind1.genotypes[chromosome_idx])
    genos2 = izip(ind2.genotypes[chromosome_idx])
    identical = [ibs(g1,g2)  for g1, g2 in izip(genos1, genos2)]

    ibd_states = np.zeros(ind1.population.chromosomes[chromosome_idx].nmark())

    # First get the segments that are IBD=1
    ibd1 = _process_segments(identical, predicate=lambda x: x > 0 or x is None})
    for start, stop in ibd1:
        ibd_states[start:(stop+1)] = 1

    # Then the segments that are IBD=2
    ibd2 = _process_segments(identical, predicate=lambda x: x in {2, None})
    for start, stop in ibd2:
        ibd_states[start:(stop+1)] = 2
    
    return ibd_states

def _process_segments(identical, predicate=lambda x: x)    
    # IBD segments are long runs of identical genotypes
    ibd = runs(identical, lambda x: x, minlength=min_seg)
    
    # Genotype errors are things that happen. If theres a small gap between
    # two IBD segments, we'll chalk that up to a genotyping error and join
    # them together.
    ibd = join_gaps(ibd, max_gap=2)

    return ibd

# Support functions

def join_gaps(seq, max_gap=1):
    seq = iter(seq)
    # Get the first item
    prev_start, prev_stop = seq.next()
    for start, stop in seq:
        if start - prev_stop > max_gap:
            yield prev_start, prev_stop
            prev_start, prev_stop = start, stop
        else:
            prev_stop = stop
    yield prev_start, stop 

