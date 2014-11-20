#!/usr/bin/env python

# Packages we're going to use
import random
import math

# Abstract base class for population and pedigree
from collections import MutableMapping

# Other pydigree objects
from pydigree.individual import Individual
from pydigree.chromosome import Chromosome
from pydigree.common import *
from pydigree.recombination import recombine

# C extension functions
from pydigree._pydigree import sample_with_replacement, random_pairs
from pydigree._pydigree import choice_with_probs

missing_genotype = (0, 0)
def is_missing_genotype(g):
    return g == missing_genotype


# Population growth models.
# These are easy enough to supply your own if you wanted
def exponential_growth(p, r, t):
    """
    Models exponential growth over discrete generations.
    p: initial population
    r: growth rate
    t: number of generations
    """
    return p * math.exp(r * t)


def logistic_growth(p, r, k, t):
    """
    Models logistic growth over discrete generations.
    p: initial population
    r: growth rate
    k: final population
    t: number of generations
    """
    return (p * k) / (p + (k - p) * math.exp(-r * t))


# Classes
class Population(MutableMapping):
    # Methods for mapping types

    def __init__(self, intial_pop_size=0, name=None):
        self.chromosomes = []
        self.pool = None
        self.population = {}
        self.n0 = intial_pop_size
        self.name = name

    # Things I have to implement for the ABC
    #
    #
    def __iter__(self):
        return iter(self.population.values())

    def __getitem__(self, key):
        return self.population[key]

    def __contains__(self, item):
        return item in self.population.values()

    def __len__(self):
        return len(self.population)

    def __setitem__(self, key, value):
        self.population[key] = value

    def __delitem__(self, key):
        del self.population[key]

    def size(self):
        """ Returns the number of individuals in the population. """
        return len(self.population)

    def remove_ancestry(self):
        """ Makes every individual in the population a founder """
        for x in self:
            x.remove_ancestry()

    # Adding and removing people
    #
    #
    def register_individual(self, ind):
        if ind.id in self.population:
            raise ValueError('ID %s already in population!' % ind.id)
        self.population[ind.id] = ind

    def remove_individual(self, ind):
        del self[ind.id]

    def random_identifier(self):
        idx = 0
        while True:
            id = str(idx)
            if self.name:
                id = name + id
            if id not in self.population.keys():
                yield id
            idx += 1
    def individuals(self):
        return [x for x in self]
    def males(self):
        """ Returns list of males in population """
        return [x for x in self if x.sex == 0]

    def females(self):
        """ Returns list of females in population """
        return [x for x in self if x.sex == 1]

    def founders(self):
        """ Returns a list of founders in population """
        return [x for x in self if x.is_founder()]

    def nonfounders(self):
        """ Returns a list of founders in population """
        return [x for x in self if not x.is_founder()]

    # Chromosome functions
    #
    #
    def add_chromosome(self, chrom):
        """ Adds a chromosome to the population """
        self.chromosomes.append(chrom)

    def chromosome_count(self):
        """ Returns the number of chromosomes in this population """
        return len(self.chromosomes)

    # Random mating
    #
    #
    def mate(self, ind1, ind2, id):
        """
        Creates an individual as the child of two specificied individual
        objects and randomly chooses a sex.
        """
        child = Individual(self, id, ind1, ind2, random.choice([0, 1]))
        return child

    def random_mating_generation(self, gensize):

        def rand_mate(pop, name=None):
            name = str(name) if not self.name else self.name + str(name)
            return self.mate(random.choice(self.males()),
                             random.choice(self.females()),
                             name)

        newpop = {}
        for i in xrange(gensize):
            newind = rand_mate(self, name=i)
            newpop[newind.id] = newind
        self.population = newpop

    def founder_individual(self, register=True, sex=None):
        if sex is not None:
            sex = sex.lower()
        sexd = {'m': 0, 'f': 1, None: random.choice([0, 1])}
        i = Individual(self, self.size(), None, None, sexd[sex])
        if register:
            self.register_individual(i)
        return i

    # Genotype Functions
    #
    #
    def get_founder_genotypes(self):
        if not self.pool:
            raise ValueError('Nothing in chromosome pool')
        g = []
        for i, x in enumerate(self.pool):
            g.append([random.choice(x), random.choice(x)])
        return g

    def get_genotypes(self):
        for x in self:
            x.get_genotypes()

    def get_linkage_equilibrium_genotypes(self):
        return [[c.linkageequilibrium_chromosome(),
                 c.linkageequilibrium_chromosome()]
                for c in self.chromosomes]

    def clear_genotypes(self):
        """ Clears genotypes for every one in population """
        for x in self:
            x.clear_genotypes()

    # Frequency functions
    #
    def missingness(self, location):
        """
        Returns the percentage of individuals in the population missing a
        genotype at location.
        """
        genotypes = [x.get_genotype(location) for x in self]
        tab = table(genotypes)
        return tab[missing_genotype] / float(len(genotypes))

    def alleles(self, location, constraint=None):
        """
        Returns the set of available alleles at a locus in the population.

        The argument constraint is a function that acts on an individual. If
        constraint(individual) can evaluate as True that person is included
        """
        if not constraint:
            constraint = lambda x: True
        gen = (x for x in self if constraint(x))
        alleles = reduce(set.union, (set(x.get_genotype(location))
                                     for x in gen if x.has_genotypes())) - {0}
        return alleles

    def allele_list(self, location, constraint=None):
        """
        The list of available alleles at a location in this population

        The argument constraint is a function that acts on an individual. If
        constraint(individual) can evaluate as True that person is included
        """
        if not constraint:
            constraint = lambda x: True
        gen = (x for x in self if constraint(x))
        alleles = list(flatten(x.get_genotype(location)
                               for x in gen if x.has_genotypes()))
        return [x for x in alleles if x != 0]

    def allele_frequency(self, location, allele, constraint=None):
        """
        Returns the frequency (as a percentage) of an allele in this population

        The argument constraint is a function that acts on an individual. If
        constraint(individual) can evaluate as True that person is included
        """
        alleles = self.allele_list(location, constraint=constraint)
        freqtab = table(alleles)
        if allele not in freqtab:
            return 0
        return freqtab[allele] / float(len(alleles))

    def major_allele(self, location, constraint=None):
        """
        Returns the major allele in this population at a locus.

        Arguments
        -----
        location: the position to find the major allele at
        constraint: a constraint (see population.alleles)

        Returns
        -----
        the major allele at a locus (type depends on how genotypes are stored)
        """
        alleles = self.allele_list(location, constraint=constraint)
        freqtab = table(alleles)
        return sorted(freqtab.keys(), key=lambda x: freqtab[x])[0]

    def haplotype_frequency(self, *args):
        matches = {x for x in self if
                   all(x.has_allele(loc, allele) for loc, allele in args)}
        return len(matches) / (2.0 * self.size())

    def ld(self, locusa, locusb, allelea=None, alleleb=None, method='r2'):
        """
        Returns a measure of linkage disquilibirum between alleles at two
        loci, for methods in {D, Dprime, r2}. If alleles aren't specified,
        the major alleles at each locus are chosen.

        LD metrics are functions of D, the measure of deviance from
        independence in sampling alleles. D = x11 - p1q1, where x11 is the
        frequency of allele 1 at locus p occuring with allele 1 at locus q.

        Lewontin's D' is a normalization of D, since D is dependent on
        frequency.

        D' = D / Dmax, where:
        Dmax = min(p1*q1, p2*q2) if D > 0
        Dmax = min(p2*q1,p1*q2) if D < 0

        r2 is a third metric, equal to
        D**2 / (p1 * p2 * q1 * q2)
        r2 ranges from 0 to 1

        Alleles are in linkage equilibrium if any of these metrics are 0.

        Arguments
        -----
        locusa:  The first locus
        locusb:  The second locus
        allelea: The allele at the first locus
                 (the major allele if not specified)
        alleleb: The allele at the second locus
                 (the major allele if not specified)
        method:  The metric of linkage disequilibrium
                 (r2 [default], Dprime, or D)

        Returns: a double
        """
        method = method.lower()
        if method not in set(['r2', 'dprime', 'd']):
            raise ValueError("LD metric must be one of r2, Dprime or D")
        if not allelea:
            allelea = self.major_allele(locusa)
        if not alleleb:
            alleleb = self.major_allele(locusb)
        pop = [x for x in self if
               (not is_missing_genotype(x.get_genotype(locusa))) and
               (not is_missing_genotype(x.get_genotype(locusb)))]
        # Allele 1
        p1 = self.allele_frequency(locusa, allelea)
        p2 = 1 - p1
        # Allele 2
        q1 = self.allele_frequency(locusa, allelea)
        q2 = 1 - q1
        x = sum(1 for x in pop
                if x.has_allele(locusa, allelea) and
                x.has_allele(locusb, alleleb))
        D = x - p1 * q1
        # All of these metrics are normalizations of D by dividing it by
        # some expression. If D==0 all of them will equal 0.
        if D == 0:
            return 0
        elif method == 'd':
            return D
        elif method == 'dprime':
            if D > 0:
                Dmax = min(p1 * q1, p2 * q2)
            else:
                Dmax = min(p2 * q1, p1 * q2)
            return D / Dmax
        elif method == 'r2':
            return (D ** 2) / (p1 * p2 * q1 * q2)

    # Phenotype functions
    #
    #
    def predict_phenotype(self, trait):
        """
        Shortcut function to call predict_phenotype(trait) for every individual
        in the population.

        Returns: nothing
        """
        for x in self:
            x.predict_phenotype(trait)
