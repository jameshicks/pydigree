from itertools import chain
from operator import add

from pydigree.common import table

class IndividualContainer(object):

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
        return set(reduce(add, [x.phenotypes.keys() for x in
                                self.individuals]))

    def phenotype_dataframe(self, onlyphenotyped=True):
        records = [x._phenotypes_to_series() for x in self.individuals]
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
        for x in self.individuals:
            x.clear_genotypes()

    def get_founder_genotypes(self):
        for x in self.founders():
            x.get_founder_genotypes()

    def get_genotypes(self):
        for x in self.individuals:
            x.get_genotypes()

    # Genotype frequency functions
    def genotype_missingness(self, location):
        """
        Returns the percentage of individuals in the population missing a
        genotype at location.
        """
        genotypes = [x.get_genotype(location) for x in self.individuals]
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
        # Reverse sort the table by the count of alleles, and return the first
        # item's first item (i.e. the allele label)
        return sorted(freqtab.items(), key=lambda x: x[1], reverse=True)[0][0]
