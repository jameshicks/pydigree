"Utilities for file I/O with SGS data"

from pydigree.sgs import SGSAnalysis, SGS, Segment
from pydigree.io import smartopen


def write_sgs(data, filename):
    """
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
    """

    with smartopen(filename, 'w') as o:
        for segment in data.segments:
            oline = []

            ind1 = segment.ind1.full_label
            ind2 = segment.ind2.full_label
            oline.extend(ind1)
            oline.extend(ind2)

            chrom = [segment.chromosome.label]
            physical = segment.physical_location
            labs = segment.marker_labels
            nmark = [segment.nmark]
            psize = [segment.physical_size / 1e6]  # Megabases, not basepairs
            oline.extend(chrom)
            oline.extend(physical)
            oline.extend(labs)
            oline.extend(nmark)
            oline.extend(psize)
            unit = ['MB']
            # Extra info GERMLINE gives you like mismatch rate
            misc = 'X', 'X', 'X'
            oline.extend(unit)
            oline.extend(misc)

            oline = '\t'.join([str(x) for x in oline])

            o.write(oline)
            o.write('\n')


class GermlineRecord(object):

    " A class for working with records in GERMLINE formatted files"
    def __init__(self, line):
        """
        Create the record from a line
        """
        l = line.strip().split()
        self.ind1 = tuple(l[0:2])
        self.ind2 = tuple(l[2:4])
        self.chromosome = l[4]

        self.unit = l[11]

        # Pick the function to convert starts and stops into the appropriate
        # numerical type: int for physical positions (bp), float for genetic
        # locations (cm)
        numberfy = float if self.unit.lower() == 'cm' else int

        self.start, self.stop = [numberfy(x) for x in l[5:7]]

    @property
    def pair(self):
        "Returns the individuals for the segment"
        return frozenset([self.ind1, self.ind2])

    @property
    def location(self):
        "Returns the location of the segment"
        return (self.start, self.stop)

    @property
    def bp_locations(self):
        """
        Are the segments provided in physical or genetic positions?

        :returns: True if physical, False if centimorgans
        :rtype: bool
        """
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
        10) Length of segment
        11) Units for genetic length (cM or MB)
        12) Mismatching SNPs in segment
        13) 1 if Individual 1 is homozygous in match; 0 otherwise
        14) 1 if Individual 2 is homozygous in match; 0 otherwise

    This function only uses 0-6.
    '''
    analysis = SGSAnalysis()
    with smartopen(filename) as f:
        for line in f:
            rec = GermlineRecord(line)

            if rec.pair not in analysis:
                analysis[rec.pair] = SGS(rec.ind1, rec.ind2)

            phys_loc = (rec.location if rec.bp_locations else None)
            seg = Segment(rec.ind1, rec.ind2, rec.chromosome, None, None,
                          physical_location=phys_loc)
            
            analysis[rec.pair].append(seg)
    return analysis
