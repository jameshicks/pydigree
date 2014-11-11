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
        if phaseknown:
            segments = [for haplotype_pair in product((0,1), repeat=2):
                

def _sgs_unphased(ind1, ind2, chromosome_idx):
    min_seg = 100
    
    genos1 = izip(ind1.genotypes[chromosome_idx])
    genos2 = izip(ind2.genotypes[chromosome_idx])
    identical = [ibs(g1,g2)  for g1, g2 in izip(genos1, genos2)]

    return _process_segments(identical)


def _process_segments(identical)    
    # IBD segments are long runs of identical genotypes
    ibd = runs(identical, lambda x: x, minlength=min_seg)
    
    # Genotype errors are things that happen. If theres a small gap between
    # two IBD segments, we'll chalk that up to a genotyping error and join
    # them together.
    ibd = join_gaps(ibd, max_gap=2)

    return ibd


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


            
        
