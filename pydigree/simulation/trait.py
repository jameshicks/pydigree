"Synthetic quantitative genetic traits"
import numpy as np

from pydigree.io import smartopen
from pydigree.genotypes import ChromosomeTemplate


class QuantitativeGeneticEffect(object):

    '''
    QuantitativeGeneticEffect is a class for objects that relate loci to 
    phenotypes.

    :ivar a: the additive effect of an allele
    :type a: numeric

    :ivar k: the deviation from additivity for the heterozygote
    :type k: numeric

    :ivar chromosomes: Chromosomes to get frequencies from
    :type chromosomes: ChromosomeSet
    '''

    def __init__(self, locus, a, k=0, chromosomes=None):
        """
        Create an effect.

        :param locus: the chromosome and marker index of the locus 
            causing the effect 
        :param a: additive effect
        :param k: dominance effect
        :param chromosomes: chromosomes for the population the effect acts on
        :type locus: tuple
        :type a: float
        :type k: float
        :type chromosomes: ChromosomeSet
        """
        self.locus = locus
        self.chromosomes = chromosomes
        self.a = a
        self.k = k

    def genotypic_value(self, individual):
        """
        The value that this individual gets from this genotype:

        ============= ======
        Minor Alleles Effect
        ============= ======
        0             0
        1             (1+k)a
        2             2a
        ============= ======

        :param individual:
        :type individual: Individual
        """
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
        """
        The expected genotypic value at the locus. This is the mean of the 
        three genotype effects, weighted by their frequency.

        :rtype: float
        """
        chridx, locidx = self.locus
        
        q = self.chromosomes[chridx].frequencies[locidx]
        p = 1 - q

        mu_g = 0

        # Genotypic value for major homozygote. Since we're counting minor
        # alleles, this always evaluates to 0, and we'll skip calculating it.
        # 
        # mu_g += p ** 2 * 0  
        
        # Heterozygote genotypic value
        mu_g += 2 * q * p * self.a * (1 + self.k)  
        
        # Genotyping value for minor homozygote
        mu_g += (q ** 2) * 2 * self.a  # Genotypic value for minor homozygote

        return mu_g

    @property
    def alpha(self):
        r"""
        Returns the average effect of allelic substitution for the locus,
        :math:`\alpha = a(1 + k(p - q))`.

        "[alpha] represents the average change in genotypic value that results 
        when a [minor allele] is randomly substituted for a major allele"

        Lynch & Walsh, Genetics and Analysis of Quantitative Traits, p. 68 

        :rtype: float
        """
        if not self.chromosomes:
            raise ValueError('Chromosomes not specified')

        chridx, locidx = self.locus
        majfreq = 1.0 - self.chromosomes[chridx].frequencies[locidx]
        return self.a * (1.0 + self.k * (majfreq - (1 - majfreq)))

    @property
    def locus_additive_variance(self):
        """
        Returns the additive variance due to the locus, based on 
        the frequencies in the chromosomeset

        :math:`\sigma^2_a = 2pq \alpha ^2`

        :returns: :math:`\sigma ^{2}_{A_{x}}`
        :rtype: float
        """
        # See Lynch & Walsh p. 69
        if not self.chromosomes:
            raise ValueError('Chromosomes not specified')

        chridx, locidx = self.locus
        majfreq = 1.0 - self.chromosomes[chridx].frequencies[locidx]

        return 2 * majfreq * (1 - majfreq) * self.alpha ** 2

    @property
    def locus_dominance_variance(self):
        """
        Returns the dominance variance due to the locus, based on 
        the frequencies in the chromosomeset

        :math:`\sigma^2_d = (2pqak)^2`

        :returns: :math:`\sigma ^{2}_{D_{x}}`
        :rtype: float
        """
        # See Lynch & Walsh p. 69
        if not self.chromosomes:
            raise ValueError('Chromosomes not specified')

        chridx, locidx = self.locus
        majfreq = 1.0 - self.chromosomes[chridx].frequencies[locidx]
        return (2 * majfreq * (1 - majfreq) * self.a * self.k) ** 2


class QuantitativeTrait(object):

    """
    QuantitativeTrait is a static class that relates genotypes to phenotypes.

    For main (i.e non-epistatic) effects, you supply a tuple of
    (chromosome, position) to indicate location. 

    Then you would add the effect to the trait architecture with:
    t.add_effect(location, a, k)

    When h2 is specified, phenotypes have an appropriate amount of random
    normal noise added to them so that the heritability of the trait in a 
    population in Hardy-Weinberg equilibrium (infinitely large, randomly 
    mating, no migration/selection, etc) equals h2. 

    The trait can be forced to have a mean at a desired location by computing 
    trait_mean = intercept + expected_genotypic_value. The genotypic value can
    be directly calculated from the specified effects, and the intercept is 
    automatically calculated by the formula.

    Predicted phenotypes for individual objects can be given 
    by t.predict_phenotype(individual)
    """

    def __init__(self, name, traittype, h2=1.0, mean=0, chromosomes=None):
        self.name = name
        self.chromosomes = chromosomes
        self.effects = []
        self.h2 = h2
        self.mean = mean
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
        """
        Sets value above which an individual is considered affected

        :param threshold: liability threshold
        :type threshold: numeric
        """
        if self.traittype != 'dichotomous':
            raise ValueError('Thresholds only for dichotomous traits')
        else:
            self.liability_threshold = threshold

    def add_effect(self, locus, a=0, k=0, effect=None):
        """
        Add a main genetic effect for a locus.

        :param locus:  a chromosome, index tuple
        :param a: The additive effect of each minor allele
        :param k: The dominance effect of at the locus, where k is the 
                deviation from additivity (default 0)
        :param effect: a QuantitativeGeneticEffect object if a and k
            are not supplied
        """

        if effect is None:
            eff = QuantitativeGeneticEffect(locus,
                                            a,
                                            k,
                                            chromosomes=self.chromosomes)
        self.effects.append(eff)

    @property
    def expected_genotypic_value(self):
        """
        The expected genotypic value for an individual in the population.
        
        This is the sum of the expected genotypic values for each of the effect
        loci

        :rtype: float
        """
        return sum(x.expected_genotypic_value for x in self.effects)

    @property
    def intercept(self):
        """
        The mean phenotypic value in the population
        """
        # E[P_environment] == 0, so we can ignore it
        return self.mean - (self.expected_genotypic_value if self.h2 < 1.0
                            else 0)

    @property
    def additive_genetic_variance(self):
        """
        Total additive variance.

        The sum of the additive variance of each effect.

        :rtype: float
        """
        return sum(x.locus_additive_variance for x in self.effects)

    @property
    def environmental_variance(self):
        """
        When h2 is fixed, this returns the corresponding environmental 
        variance in the population that must be added to result in the 
        specified h2 value.

        :rtype: float
        """
        if self.h2 is None:
            raise ValueError('Trait heritability not set!')

        # h2 = V_a / (V_a + V_e)
        # A little algebra gives us V_e =  V_a/h2 - V_a
        add = self.additive_genetic_variance
        return (add / self.h2) - add

    @property
    def total_variance(self):
        """
        The total variance of the phenotype

        :rtype: float
        """
        return self.additive_genetic_variance + self.environmental_variance

    def rescale(self, mean, sd):
        """
        Rescale trait to have the distribution N(mean, sd). Modifies genotype
        effects so that they sum to the appropriate variance.

        :param mean: the desired mean
        :param sd: the desired sd
        :rtype: void   
        """

        old_add_variance = self.additive_genetic_variance
        new_add_variance = self.h2 * (sd**2)

        # Scaling works easiest when done as standard deviations
        scaling_factor = np.sqrt(new_add_variance) / np.sqrt(old_add_variance)
        for effect in self.effects:
            effect.a *= scaling_factor

        self.mean = mean

    def predict_phenotype(self, individual):
        """
        Generates a predicted phenotype for an individual based off their 
        genotype

        :param individual: subject to have a phenotype prediected


        :returns: Trait value if quantitative, affectation status if dichotomous
        :rtype: double (quantitative); 0 or 1 (dichotomous)
        """

        phenotype = [eff.genotypic_value(individual) for eff in self.effects]
        phenotype = self.intercept + sum(phenotype)

        # Random variates can't be generated for variance <= 0,
        # so we have to check if there even is an environmental
        # component to the trait (i.e. h2 < 1) for the trait we're
        # simulating before we try to add it to the phenotype
        if self.h2 < 1.0:
            enviro = self.environmental_variance
            phenotype += np.random.normal(0, np.sqrt(enviro))

        if self.traittype == 'dichotomous':
            if self.liability_threshold is None:
                raise ValueError('No liability threshold set')
            return 1 if phenotype >= self.liability_threshold else 0
        return phenotype

    def add_dummy_polygene_chromosomes(self, population, nloc,
                                       mean=0,
                                       sd=1,
                                       freqs=None,
                                       polylabel='Polygene'):
        """
        Creates many independently segregating chromosomes that 
        additively influence the trait.

        :param population: The population to add the chromosomes to
        :param nloc:       The number of dummy chromosomes to create
        :param mean:       Mean locus additive effect
        :param sd:         Standard deviation of locus additive effect
        :param freqs: frequencies of each polygene
        :param polylabel:  Label to add give chromosome

        :type nloc: integer
        :type mean: float
        :type sd: float
        :type freqs: sequence of floats
        :type polylabel: string

        :rtype: void
        """
        if freqs is None:
            freqs = np.zeros(nloc, dtype=np.float) + 0.5

        if sd == 0:
            effects = [mean] * nloc
        else:
            effects = np.random.normal(mean, sd, nloc)

        for i, effect in enumerate(effects):
            # Create the chromosome
            lab = '{}{}'.format(polylabel, i)
            c = ChromosomeTemplate(label=lab)
            c.add_genotype(freqs[i], 0)
            population.add_chromosome(c)

            # Add the effect
            locus = i, 0
            self.add_effect(locus, a=effect, k=0)

    @staticmethod
    def from_file(filename):
        """
        Reads a trait from a file

        :param filename: path to file
        :type filename: string

        :rtype: QuantitativeTrait
        """
        with smartopen(filename) as f:
            trait_type, name = f.readline().strip().split()
            trait = QuantitativeTrait(trait_type, name)
            for line in f:
                l = line.strip().split()
                
                if len(l) != 5:
                    # TODO: implement epistatic effects in file
                    raise NotImplementedError(
                        'Epistatic effects not yet implemented')
                chrom, loc, _, _, a, k = line.strip().split()
                locus = chrom, loc
                trait.add_effect(locus, a, k)
        
        return trait
