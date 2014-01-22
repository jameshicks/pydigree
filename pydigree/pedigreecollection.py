#!/usr/bin/env python

from common import *
from operator import add
from pedigree import Pedigree
from scipy.sparse import block_diag
from collections import MutableMapping


class PedigreeCollection(MutableMapping):

    def __init__(self):
        self.pedigrees = {}

    ### Things I have to implement for the ABC
    ###
    def __iter__(self):
        return (x for x in self.pedigrees.values())

    def __getitem__(self, key):
        return self.pedigrees[key]

    def __contains__(self, item):
        return item in self.pedigrees.values()

    def __len__(self):
        return len(self.pedigrees)

    def __setitem__(self, key, value):
        # All pedigrees have to have the same population
        if not all(ped.population is value.population for ped in self):
            raise ValueError("Pop for %s doesn't match collection" % value)

        self.pedigrees[key] = value

    def __delitem__(self, key):
        del self.pedigrees[key]

    def keys(self):
        return self.pedigrees.keys()

    def individuals(self):
        inds = []
        for pedigree in sorted(self, key=lambda x: x.label):
            inds.extend(sorted((x for x in pedigree), key=lambda x: x.id))
        return inds

    def phenotypes(self):
        """ Returns the available phenotypes for analysis """
        return set(reduce(add, [x.phenotypes.keys() for x in
                                self.individuals()]))

    def population(self):
        """ Returns the population all the pedigrees belong to. """
        k = self.pedigrees.keys()[0]
        return self.pedigres[k]

    def chromosomes(self):
        return self.population.chromosomes

    def add_chromosome(self, chrom):
        for x in self:
            x.add_chromosome(chrom)

    def clear_genotypes(self):
        for x in self:
            x.clear_genotypes()

    def get_founder_genotypes(self):
        for x in self:
            x.get_founder_genotypes()

    def get_genotypes(self):
        for x in self:
            x.get_genotypes()

    ### Matrix functions
    ###
    def additive_relationship_matrix(self):
        """
        Returns a block diagonal matrix of additive relationships
        for each pedigree.

        See notes on Pedigree.additive_relationship_matrix
        """
        return block_diag([x.additive_relationship_matrix() for x in
                           sorted(self, key=lambda x: x.label)], format='bsr')

    def dominance_relationship_matrix(self):
        """
        Returns a block diagonal matrix of dominance relationships
        for each pedigree.

        See notes on Pedigree.dominance_relationship_matrix
        """
        return block_diag([x.dominance_relationship_matrix() for x in
                           sorted(self, key=lambda x: x.label)], format='bsr')

    def mitochondrial_relationship_matrix(self):
        """
        Returns a block diagonal matrix of mitochondrial relationships
        for each pedigree.

        See notes on Pedigree.mitochondrial_relationship_matrix
        """
        return block_diag([x.mitochondrial_relationship_matrix() for x in
                           sorted(self, key=lambda x: x.label)], format='bsr')
