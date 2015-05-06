#!/usr/bin/env python


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

    Epistatic Effects
    ------
    Pydigree can also simulate epistatic effects. That is, combinations of
    genotypes that add effects non additively. For these the locations are
    given as a list of (chr, pos) tuples (as seen in the main effects section).
    The dictionary of effects is given as the key,value pair of a list of
    genotypes and the effect of that genotype on the trait.

    The list of genotypes reflects the order of the positions. That is, if you
    supply positions 1,2,3, the genotypes should look like:
    [(genotype at position 1), (genotype at position 2),
    (genotype at position 3)]

    Like with the main effects, any genotype combination not in the effect
    dictionary is assumed to have 0 effect on the trait.
    """

    def __init__(self, name, type):
        self.name = name
        self.effects = {}
        self.epistatic_effects = {}
        if type not in ['quantitative', 'dichotomous']:
            raise ValueError('Not a valid trait type!')
        else:
            self.traittype = type
        self.liability_threshold = None

    def __str__(self):
        return "Trait %s (%s): %s main effects, %s epistatic effects" \
            % (self.name, self.traittype,
               len(self.effects), len(self.epistatic_effects))

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

    def _getepistaticeffects(self, individual, locations):
        locations = tuple(tuple(loc) for loc in locations)
        g = tuple(frozenset(individual.get_genotype(chr, pos)
                            for chr, pos in locations))
        if g not in self.epistatic_effects[locations]:
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
        """
        chrom, marker = location
        effects = self.__transformdict(effects)
        try:
            self.effects[location].update(effects)
        except KeyError:
            self.effects[(chrom, marker)] = effects

    def add_epistatic_effect(self, locations, effects):
        locations = tuple(tuple(x) for x in locations)
        effects = self.__transformdict(effects)
        self.epistatic_effects[locations] = effects

    def predict_phenotype(self, individual):
        liability = [self._geteffect(individual, x) for x in self.effects]
        liability += [self._getepistaticeffect(individual, x)
                      for x in self.epistatic_effects]
        liabilitysum = sum(liability)
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
