import numpy as np

class Architecture(object):

    """
    Architecture is a static class that relates genotypes to phenotypes.

    Main Effects
    ------
    For main (i.e non-epistatic) effects, you supply a tuple of
    (chromosome, position) to indicate location. You also supply a dictionary
    of genotype/effect pairs. For example, to add a fully dominant effect of
    allele A  that adds 1 to the trait value at chromosome 4, position 50,
    you would prepare a dictionary in the form:

    effects = {('A','A'): 1, ('A','a'): 1, ('a','a'): 0}
    location = (4,50)

    Then you would add the effect to the trait architecture with:
    t.add_effect_liability(location,effects)

    For convenience, any genotype not found in the effects dictionary has no
    effect on the trait. For example a fully reccessive effect of allele 'a'
    on the trait, giving an effect of 1 could be a dictionary that looks like
    this: {('a','a'): 1}


    """

    def __init__(self, name, type, chromosomes=None):
        self.name = name
        self.chromosomes = chromosomes
        self.effects = {}
        self.noise = None
        if type not in ['quantitative', 'dichotomous']:
            raise ValueError('Not a valid trait type!')
        else:
            self.traittype = type
        self.liability_threshold = None

    def __str__(self):
        return "Trait {} ({}): {} main effects".format(self.name,
                                                       self.traittype,
                                                       len(self.effects))

    def __transformdict(self, d):
        d2 = {}
        for k in d:
            d2[frozenset(k)] = d[k]
        return d2

    def _geteffect(self, individual, location):
        g = frozenset(individual.get_genotype(location))
        try:
            return self.effects[location][g]
        except KeyError:
            return 0

    def set_liability_threshold(self, threshold):
        if self.traittype != 'dichotomous':
            raise ValueError('Thresholds only for dichotomous traits')
        else:
            self.liability_threshold = threshold

    def add_effect(self, location, effects):
        """
        Add a main effect

        Arguments
        ------
        location: a chromosome, index tuple
        effects: a dictionary describing the effect for each genotype. eg:
        Dominance: {(1,1): 1, (0,1): 1, (0,0): 0 }
        Additive: {(1,1): 2, (0,1): 1, (0,0): 0 }

        or

        a numeric variable that will be made into an additive effect
        """
        chrom, marker = location
        if not isinstance(effects, dict):
            try:
                effect = float(effects)
                effects = {(1,1): 0, (1,2): effect, (2,2): 2 * effect}
            except ValueError:
                raise ValueError('Unparsable effect: {}'.format(effects))
        effects = self.__transformdict(effects)
        try:
            self.effects[location].update(effects)
        except KeyError:
            self.effects[(chrom, marker)] = effects

    def add_noise(self, mean=0, sd=1):
        self.noise = (mean, sd)

    def predict_phenotype(self, individual):
        liability = [self._geteffect(individual, x) for x in self.effects]
        liabilitysum = sum(liability)

        if self.noise:
            mu, sigma = self.noise
            liabilitysum += np.random.normal(mu, sigma)

        if self.traittype == 'dichotomous':
            if self.liability_threshold is None:
                raise ValueError('No liability threshold set')
            return 1 if liabilitysum >= self.liability_threshold else 0
        return liabilitysum

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
                chr, pos, allele_a, allele_b, effect = line.strip().split()
                locus = chr, pos
                eff = {(allele_a, allele_b): effect}
                trait.add_effect(locus, eff)
