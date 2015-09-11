import numpy as np
from pydigree.common import rescale_variable


class QuantitativeGeneticEffect(object):

    '''
    QuantitativeGeneticEffect is a class for objects that relate loci to 
    phenotypes.

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
        N = gt.count(2)  # Number of minor alleles
        if N == 0:
            return 0
        elif N == 1:
            return self.a * (1 + self.k)
        elif N == 2:
            return 2 * self.a
        else:
            raise ValueError('Bad genotype: {}'.format(gt))

    @property
    def expected_genotypic_value(self):
        chridx, locidx = self.locus
        maf = self.chromosomes[chridx].frequencies[locidx]

        mu_g = 0
        mu_g += (1-maf)**2 * 0  # Genotypic value for major homozygote
        mu_g += 2 * maf * (1-maf) * self.a * (1 + self.k)  # Heterozygote value
        mu_g += (maf ** 2) * 2 * self.a  # Genotypic value for minor homozygote

        return mu_g

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
    t.add_effect(location, a, k)

    When h2 is specified, phenotypes have an appropriate amount of random
    normal noise added to them so that the heritability of the trait in a 
    population in Hardy-Weinberg equilibrium (infinitely large, randomly mating, 
    no migration/selection, etc) equals h2. 

    When trait_mean and trait_sd are specified, the output phenotype is shifted
    and rescaled so that it is distribution N(trait_mean, trait_sd). In the 
    presence of rescaling, the additive values of the genotype effects are 
    not particularly meaningful, since they are not rescaled internally. The
    effects then should be considered as relative weights of the effects of
    different loci, and not directly predictive.   

    Predicted phenotypes for individual objects can be given 
    by t.predict_phenotype(individual)
    """

    def __init__(self, name, traittype, h2=None, chromosomes=None,
                 trait_mean=None, trait_sd=None):
        self.name = name
        self.chromosomes = chromosomes
        self.effects = []
        self.h2 = h2

        self.trait_mean = float(trait_mean)
        self.trait_sd = float(trait_sd)

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

    def add_effect(self, locus, a=None, k=0, obj=None):
        """
        Add a main effect

        Arguments
        ------
        locus:  a chromosome, index tuple
        a:                 The additive effect of each minor allele
        k:                 The dominance effect of at the locus, where k is the 
                           deviation from additivity (default 0)
        phenotypic_effect: The additive effect of an allele on the scale of the
                           output trait
        obj:               a QuantitativeGeneticEffect object if you don't 
                           want to supply a and k
        """
        chrom, marker = locus

        if obj is None:
            if a 
            obj = QuantitativeGeneticEffect(locus,
                                            a,
                                            k,
                                            chromosomes=self.chromosomes)
        self.effects.append(obj)

    @property
    def unscaled_expected_genotypic_value(self):
        return sum(x.expected_genotypic_value for x in self.effects)

    @property
    def unscaled_additive_genetic_variance(self):
        return sum(x.locus_additive_variance for x in self.effects)

    @property
    def unscaled_environmental_variance(self):
        if self.h2 is None:
            raise ValueError('Trait heritability not set!')

        # h2 = V_a / (V_a + V_e)
        # A little algebra gives us V_e =  V_a/h2 - V_a
        add = self.unscaled_additive_genetic_variance
        return (add / self.h2) - add

    @property
    def unscaled_total_variance(self):
        return (self.unscaled_additive_genetic_variance +
                self.unscaled_environmental_variance)

    @property
    def rescaling(self):
        " Returns True if phenotypes are rescaled to a new distribution "
        return self.trait_mean is not None and self.trait_sd is not None

    def rescale_phenotype(self, phenotype, old_mean, old_sd):
        """
        Takes a phenotype value distributed N(old_mean, old_sd) and rescales 
        to N(self.trait_mean, self.trait_sd) 
        """
        new_mean = self.trait_mean
        new_sd = self.trait_sd
        newphen = new_mean + (phenotype - old_mean) * (new_sd / float(old_sd))
        newphen = rescale_variable(phenotype,
                                   old_mean, old_sd,
                                   new_mean, new_sd)
        return newphen

    def _rescale_to_genotype_effect(self, phenotypic_effect):
        """
        Rescales an effect on a phenotypic value to the scale
        genotypic values are on, if phenotypic rescaling is turned on.
        """
        if not self.rescaling:
            return effect

        unscaled_sd = np.sqrt(self.unscaled_total_variance)
        genotype_effect = phenotypic_effect * unscaled_sd / self.trait_sd

    def predict_phenotype(self, individual):
        phenotype = [eff.genotypic_value(individual) for eff in self.effects]
        phenotype = sum(phenotype)

        if self.h2:
            enviro = self.unscaled_environmental_variance
            phenotype += np.random.normal(0, np.sqrt(enviro))

        if self.rescaling:
            # P = G + E
            #
            # The phenotype is the sum of the genotypic and environmental
            # values. So the expected phenotype value is the sum of the
            # expected genotype value and the expected environmental value.
            # Since in these simulations E[environmental_value] = 0, we don't
            # have to worry about it.
            old_mean = self.unscaled_expected_genotypic_value

            phenotype = self.rescale_phenotype(phenotype,
                                               old_mean,
                                               np.sqrt(self.total_variance))

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
