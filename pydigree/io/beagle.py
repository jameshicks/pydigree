import numpy as np



from pydigree.individual import Individual
from pydigree.genotypes import ChromosomeTemplate, Alleles
from pydigree.io import smartopen as open
from pydigree.exceptions import FileFormatError
from pydigree.common import grouper


class BeagleMarkerRecord(object):
    __slots__ = ['label', 'pos', 'alleles']

    def __init__(self, line):
        name, pos, alleles = line.strip().split(None, 2)
        alleles = alleles.split()
        self.label = name
        self.pos = int(pos)
        self.alleles = alleles

    @property
    def reference(self):
        return alleles[0]

    @property
    def alternates(self):
        return allels[1:]

class BeagleGenotypeRecord(object):
    __slots__ = ['identifier', 'label', 'data']

    def __init__(self, line):
        l = line.strip().split()
        self.identifier = l[0]
        self.label = l[1]
        self.data = l[2:]

    def is_phenotype_record(self):
        return self.identifier in 'ACT'


def read_beagle_markerfile(filename, label=None):
    ''' 
    Reads marker locations from a beagle format file
    
    Arguments
    -----
    filename: The file to be read
    label: An optional label to give the chromosome, since the BEAGLE
        format does not require it

    Returns: a ChromosomeTemplate object
    '''
    with open(filename) as f:
        chrom = ChromosomeTemplate(label=label)

        last_pos = -1
        for line in f:
            rec = BeagleMarkerRecord(line)

            if rec.pos < 0:
                raise FileFormatError(
                    'Bad position for genotype: {}'.format(rec.pos))
            elif rec.pos <= last_pos:
                raise FileFormatError('Makers in file out of order')

            chrom.add_genotype(None, cm=None, label=rec.label, bp=rec.pos,
                               reference=rec.reference, alternates=rec.alternates)
            last_pos = rec.pos

    return chrom


def read_beagle_genotypefile(filename, pop, missingcode='0'):
    '''
    Reads BEAGLE formatted genotype files
    
    Arguments
    ------
    filename: Filename of BEAGLE genotype file
    pop: the population to add these individuals to
    missingcode: The value that indicates a missing genotype

    Returns: Nothing
    '''
    with open(filename) as f:
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
            ind.genotypes = genotypes_from_sequential_alleles(ind.chromosomes,
                                                              sequentialalleles,
                                                              missingcode=missingcode)


def read_beagle(genofile, markerfile):
    '''
    Reads BEAGLE formatted genotype data

    Arguments
    ------
    genofile: File containing genotype information for individuals
    markerfile: File containing marker location and allele information

    Returns: a Population object
    '''
    pop = Population()
    pop.chromosomes.append(read_beagle_markerfile(markerfile))

    read_beagle_genotypefile(genofile, pop)

    return pop
