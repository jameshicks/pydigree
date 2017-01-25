import pandas as pd
from itertools import chain
from operator import add
from functools import reduce

from pydigree.common import table


class IndividualContainer(object):

    def apply_inplace(self, func):
        ''' 
        Calls a function on each individual object in the collection

        :param func: A function that takes an :class:`Individual` as a parameter
        :type func: callable

        :rtype: void
        '''

        for ind in self.individuals:
            func(ind)

    def apply(self, func):
        '''
        Calls a function on each individual object in the collection and yields its value

        :param func: A function that takes an :class:`Individual` as a parameter
        :type func: callable

        :rtype: generator
        '''
        for ind in self.individuals:
            yield func(ind)

    # Individual filtering functions
    def males(self):
        """ Returns list of males in population """
        return [x for x in self.individuals if x.sex == 0]

    def females(self):
        """ Returns list of females in population """
        return [x for x in self.individuals if x.sex == 1]

    def founders(self):
        """ Returns a list of founders in population """
        return [x for x in self.individuals if x.is_founder()]

    def nonfounders(self):
        """ Returns a list of founders in population """
        return [x for x in self.individuals if not x.is_founder()]

    def sex_ratio(self):
        """
        Returns the sex ratio in the population, defined as n_male / n_female

        :returns: sex ratio
        :rtype: float
        """
        counts = [0, 0]
        for ind in self.individuals:
            counts[ind.sex == 0] += 1

        return counts[0] / counts[1]

    # Phenotype functions
    #
    #
    def predict_phenotype(self, trait):
        """
        Shortcut function to call predict_phenotype(trait) for every individual
        in the population.

        Returns: nothing
        """
        for x in self.individuals:
            x.predict_phenotype(trait)

    def delete_phenotype(self, trait):
        for ind in self.individuals:
            ind.delete_phenotype(trait)

    def phenotypes(self):
        """ Returns the available phenotypes for analysis """
        return set(reduce(add, [list(x.phenotypes.keys()) for x in
                                self.individuals]))

    def phenotype_dataframe(self, onlyphenotyped=True):
        '''
        Returns a pandas dataframe of the phenotypes available 
        for the individuals present in the collection.

        Arguments:
        onlyphenotyped: Only include individuals who have phenotypes

        Returns: A pandas dataframe
        '''
        records = [x.phenotypes.to_series() for x in self.individuals]
        df = pd.DataFrame.from_records(records)

        if onlyphenotyped:
            df.dropna(how='all', inplace=True)

        return df

    def genotype_as_phenotype(self, locus, minor_allele, label):
        """ 
        Dispatches a genotype_as_phenotype to each individual in the pedigree.
        
        See docstring for Individual.genotype_as_phenotype for more
        details
        """
        for ind in self.individuals:
            ind.genotype_as_phenotype(locus, minor_allele, label)

    # Genotype generation
    def clear_genotypes(self):
        "Removes all the genotypes for individuals"

        for x in self.individuals:
            x.clear_genotypes()

    def get_founder_genotypes(self):
        "Have founder individuals request genotypes"

        for x in self.founders():
            x.get_founder_genotypes()

    def get_genotypes(self):
        "Have individuals request genotypes"
        for x in self.individuals:
            x.get_genotypes()

    # Genotype frequency functions
    def genotype_missingness(self, location):
        """
        Returns the percentage of individuals in the population missing a
        genotype at location.

        Returns: A float
        """
        missing_genotype = (0,0)
        genotypes = [x.get_genotype(location) for x in self.individuals]
        tab = table(genotypes)
        return tab[missing_genotype] / float(len(genotypes))

    def alleles(self, location, constraint=None):
        """
        Returns the set of available alleles at a locus in the population.

        The argument constraint is a function that acts on an individual. If
        constraint(individual) can evaluate as True that person is included

        Returns: a set of the available genotypes at a location
        """
        if not constraint:
            constraint = lambda x: True
        gen = (x for x in self.individuals if constraint(x))
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
        gen = (x for x in self.individuals if constraint(x))
        alleles = list(chain.from_iterable(x.get_genotype(location)
                               for x in gen if x.has_genotypes()))
        return [x for x in alleles if x != 0]

    def allele_frequency(self, location, allele, constraint=None):
        """
        Returns the frequency (as a percentage) of an allele in this population

        :param location: the locus to be evaluated
        :param allele: the allele to be counted
        :param constraint: Function that acts on an individual. If
            constraint(individual) can evaluate as True that person is included
        
        :type constraint: callable
        :rtype: float 
        """
        alleles = self.allele_list(location, constraint=constraint)
        freqtab = table(alleles)
        if allele not in freqtab:
            return 0
        return freqtab[allele] / float(len(alleles))

    def major_allele(self, location, constraint=None):
        """
        Returns the major (most common) allele in this population at a locus.

        
        :param location: the position to find the major allele at
        :param constraint: a constraint (see population.alleles)

        :returns: the major allele at a locus 
        """
        alleles = self.allele_list(location, constraint=constraint)
        freqtab = table(alleles)
        # Reverse sort the table by the count of alleles, and return the first
        # item's first item (i.e. the allele label)
        return sorted(freqtab.items(), key=lambda x: x[1], reverse=True)[0][0]
