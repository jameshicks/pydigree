import numpy as np


class GeneticEffect(object):

    '''
    GeneticEffect is a class for objects that relate loci to phenotypes.

    Effects are parameterized with a and k (see also Lynch & Walsh p. 62), 
    illustrated by the following ascii art diagram:

    effect:      0                        (1+k)a            2a 
                 |---------------------------|--------------|
    genotype:  A1/A1                       A1/A2          A2/A2 
    '''

    def __init__(self, locus, a, k=0, chromosomes=None):
        self.locus = locus
        self.chromosomes = chromosomes
        self.a = a
        self.k = k

    def genotypic_value(self, individual):
        gt = individual.get_genotype(self.locus)
        eff = self.a * gt.count(2)  # Additive effect
        eff += (1+self.k) * self.a if gt[0] != gt[1] else 0  # Dominance effect
        return eff

    @property
    def alpha(self):
        'Returns the average effect of allelic substitution for the locus'
        if not self.chromosomes:
            raise ValueError('Chromosomes not specified')

        chridx, locidx = self.locus
        majfreq = 1.0 - self.chromosomes[chridx].frequencies[locidx]
        return self.a * (1.0 + self.k * (majfreq - (1 - majfreq)))

    @property
    def locus_additive_variance(self):
        # See Lynch & Walsh p. 69
        if not self.chromosomes:
            raise ValueError('Chromosomes not specified')

        chridx, locidx = self.locus
        majfreq = 1.0 - self.chromosomes[chridx].frequencies[locidx]

        return 2 * majfreq * (1 - majfreq) * self.alpha ** 2

    @property
    def locus_dominance_variance(self):
        # See Lynch & Walsh p. 69
        if not self.chromosomes:
            raise ValueError('Chromosomes not specified')

        chridx, locidx = self.locus
        majfreq = 1.0 - self.chromosomes[chridx].frequencies[locidx]
        return (2 * majfreq * (1 - majfreq) * self.a * self.k) ** 2


class Architecture(object):

    """
    Architecture is a static class that relates genotypes to phenotypes.

    Main Effects
    ------
    For main (i.e non-epistatic) effects, you supply a tuple of
    (chromosome, position) to indicate location. 

    Then you would add the effect to the trait architecture with:
    t.add_effect_liability(location, a, k)

    Predicted phenotypes for individual objects can be given 
    by t.predict_phenotype(individual)
    """

    def __init__(self, name, traittype, chromosomes=None):
        self.name = name
        self.chromosomes = chromosomes
        self.effects = []
        self.noise = None
        if traittype not in ['quantitative', 'dichotomous']:
            raise ValueError('Not a valid trait type!')
        else:
            self.traittype = traittype
        self.liability_threshold = None

    def __str__(self):
        return "Trait {} ({}): {} main effects".format(self.name,
                                                       self.traittype,
                                                       len(self.effects))


    def set_liability_threshold(self, threshold):
        if self.traittype != 'dichotomous':
            raise ValueError('Thresholds only for dichotomous traits')
        else:
            self.liability_threshold = threshold

    def add_effect(self, locus, a=0, k=0, effect=None):
        """
        Add a main effect

        Arguments
        ------
        location: a chromosome, index tuple
        a: The additive effect of each minor allele
        k: The dominance effect of at the locus, where k is the 
           deviation from additivity (default 0)
        effect: a GeneticEffect object if you don't want to supply a and k
        """
        chrom, marker = location

        if effect is None:
            eff = GeneticEffect(locus, a, k, chromosomes=self.chromosomes)
        self.effects.append(eff)

    def add_noise(self, mean=0, sd=1):
        self.noise = (mean, sd)

    @property
    def additive_genetic_variance(self):
        return sum(x.locus_additive_variance for x in self.effects)

    def predict_phenotype(self, individual):
        phenotype = [eff.genotypic_value(individual) for eff in self.effects]
        phenotype = sum(phenotype)

        if self.noise:
            mu, sigma = self.noise
            phenotype += np.random.normal(mu, sigma)

        if self.traittype == 'dichotomous':
            if self.liability_threshold is None:
                raise ValueError('No liability threshold set')
            return 1 if phenotype >= self.liability_threshold else 0
        return phenotype

    @staticmethod
    def from_file(filename):
        with open(filename) as f:
            type, name = f.readline().strip().split()
            trait = Architecture(type, name)
            for line in f:
                l = line.strip().split()
                if len(l) != 5:
                    # TODO: implement epistatic effects in file
                    raise NotImplementedError(
                        'Epistatic effects not yet implemented')
                chrom, loc, allele_a, allele_b, a, k = line.strip().split()
                locus = chrom, loc
                trait.add_effect(locus, a, k)
