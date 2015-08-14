from pydigree.sgs import SGSAnalysis, SGS, Segment
from pydigree.io import smartopen as open

def write_sgs(data, filename):
    with open(filename, 'w') as o:
        for pair, intervals in data.items():
            for i, chromintervals in enumerate(intervals):
                if not chromintervals:
                    continue
                ind1, ind2 = pair
                l = [ind1.pedigree.label, ind1.label,
                     ind2.pedigree.label, ind2.label, i] + chromintervals
                o.write('\t'.join(l) + '\n')


class GermlineRecord(object):

	" A class for working with records in GERMLINE formatted files"
    def __init__(self, payload):
        l = payload.strip().split()
        self.ind1 = l[0:2]
        self.ind2 = l[2:4]
        self.chromosome = l[4]

        self.unit = l[11]

        # Pick the function to convert starts and stops into the appropriate
        # numerical type: int for physical positions (bp), float for genetic
        # locations (cm)
        numberfy = float if self.unit.lower() == 'cm' else int

        self.start, self.stop = [int(x) for x in l[5:7]]

    @property
    def pair(self):
        return frozenset([self.ind1, self.ind2])

    @property
    def location(self):
        return (self.start, self.stop)

    @property
    def bp_locations(self):
        return self.unit.lower() == 'mb'


def read_germline(filename):
    '''
    Reads a GERMLINE formatted SGS filename into an SGSAnalysis object

    GERMLINE files are text files with the format:

    	0) Family ID 1
        1) Individual ID 1
        2) Family ID 2
        3) Individual ID 2
        4) Chromosome
        5) Segment start (bp/cM)
        6) Segment end (bp/cM)
        7) Segment start (SNP)
        8) Segment end (SNP)
        9) Total SNPs in segment
        10) Genetic length of segment
        11) Units for genetic length (cM or MB)
        12) Mismatching SNPs in segment
        13) 1 if Individual 1 is homozygous in match; 0 otherwise
        14) 1 if Individual 2 is homozygous in match; 0 otherwise

    This function only uses 0-6 and 11.
    '''
    analysis = SGSAnalysis()
    with open(filename) as f:
        for line in f:
            rec = GermlineRecord(line)

            if rec.pair not in analysis:
                analysis[rec.pair] = SGS()

            phys_loc = (rec.location if rec.bp_locations else None)
            gen_loc = (rec.location if not rec.bp_locations else None)
            seg = Segment(rec.ind1, rec.ind2, rec.chromosome, None, None,
                          physical_location=phys_loc,
                          genetic_location=gen_loc)
            
            analysis[rec.pair].append(seg)
    return analysis
