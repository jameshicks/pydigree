"Utilities for reading BEAGLE formatted genotype data"

from pydigree.population import Population
from pydigree.individual import Individual
from pydigree.genotypes import ChromosomeTemplate
from pydigree.io import smartopen 
from pydigree.io.base import genotypes_from_sequential_alleles as gt_from_seq
from pydigree.exceptions import FileFormatError
from pydigree.common import grouper


class BeagleMarkerRecord(object):
    "A class corresponding to a marker in a BEAGLE marker file"
    __slots__ = ['label', 'pos', 'alleles']

    def __init__(self, line):
        """
        Create the marker record

        :param line: record line
        :type line: string
        """
        name, pos, alleles = line.strip().split(None, 2)
        alleles = alleles.split()
        self.label = name
        self.pos = int(pos)
        self.alleles = alleles

    @property
    def reference(self):
        """
        The reference allele for the marker

        :returns: reference allele
        :rtype: string
        """
        return self.alleles[0]

    @property
    def alternates(self):
        """
        The alternate alleles for the marker
        
        :returns: Each alternate allele
        :rtype: iterable
        """
        return self.alleles[1:]


class BeagleGenotypeRecord(object):
    "A class corresponding to a record in a BEAGLE genotype file"
    __slots__ = ['identifier', 'label', 'data']

    def __init__(self, line):
        """
        Create the genotype record object.

        :param line: record line
        :type line: string
        """
        l = line.strip().split()
        self.identifier = l[0]
        self.label = l[1]
        self.data = l[2:]

    def is_phenotype_record(self):
        "Returns true if record corresponds to a phenotype field"
        return self.identifier in 'ACT'


def read_beagle_markerfile(filename, label=None):
    """ 
    Reads marker locations from a BEAGLE formatted file
    
    :param filename: The file to be read
    :param label: An optional label to give the chromosome, since the BEAGLE
        format does not require it
    
    :type filename: string

    :rtype: ChromosomeTemplate
    """
    with smartopen(filename) as f:
        chrom = ChromosomeTemplate(label=label)

        last_pos = -1
        for line in f:
            rec = BeagleMarkerRecord(line)

            if rec.pos < 0:
                raise FileFormatError(
                    'Bad position for genotype: {}'.format(rec.pos))
            elif rec.pos <= last_pos:
                raise FileFormatError('Makers in file out of order')

            chrom.add_genotype(None, map_position=None, label=rec.label, 
                               bp=rec.pos, reference=rec.reference, 
                               alternates=rec.alternates)
            last_pos = rec.pos

    return chrom


def read_beagle_genotypefile(filename, pop, missingcode='0'):
    '''
    Reads BEAGLE formatted genotype files
    
    Arguments

    :param filename: Filename of BEAGLE genotype file
    :param pop: the population to add these individuals to
    :param missingcode: The value that indicates a missing genotype
    
    :type missingcode: string
    :rtype: void
    '''
    with smartopen(filename) as f:
        for line in f:
            rec = BeagleGenotypeRecord(line)

            if rec.identifier == 'I':
                inds = [Individual(pop, label) for label in rec.data[::2]]
            elif rec.is_phenotype_record:
                for ind, pheno_status in zip(inds, rec.data[::2]):
                    if rec.identifier == 'A':
                        pheno_status = pheno_status == '2'
                    else:
                        try:
                            pheno_status = float(pheno_status)
                        except ValueError:
                            pass
                    ind.phenotypes[rec.label] = pheno_status
            else:
                # We've reached the genotypes, and we're skipping out
                break
        f.seek(0)
        gtrows = [list(grouper(BeagleGenotypeRecord(x).data, 2))
                  for x in f if x.startswith('M')]
        genotypes = zip(*gtrows)
        for ind, sequentialalleles in zip(inds, genotypes):
            ind.genotypes = gt_from_seq(ind.chromosomes,
                                        sequentialalleles,
                                        missing_code=missingcode)


def read_beagle(genofile, markerfile):
    '''
    Reads BEAGLE formatted genotype data

    :param genofile: Filename containing genotype information for individuals
    :param markerfile: Filename containing marker location and allele 
        information corresponding to genofile

    :type genofile: string
    :type markerfile: string

    :rtype: Population
    '''
    pop = Population()
    chrom = read_beagle_markerfile(markerfile)
    chrom.finalize()
    pop.chromosomes.add_chromosome(chrom)

    read_beagle_genotypefile(genofile, pop)

    return pop
