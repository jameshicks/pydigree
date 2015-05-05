#!/usr/bin/env python

from common import *
from operator import add
from collections import MutableMapping

from scipy.sparse import block_diag

from pydigree.pedigree import Pedigree


class PedigreeCollection(MutableMapping):

    def __init__(self):
        self.pedigrees = {}

    ### Things I have to implement for the ABC
    ###
    def __iter__(self):
        return (x for x in self.pedigrees.values())

    def __getitem__(self, key):
        if isinstance(key, tuple) or isinstance(key, list):
            return self.pedigrees[key[0]][key[1]]
        return self.pedigrees[key]

    def __contains__(self, item):
        return item in self.pedigrees.values()

    def __len__(self):
        return len(self.pedigrees)

    def __setitem__(self, key, value):
        self.pedigrees[key] = value

    def __delitem__(self, key):
        del self.pedigrees[key]

    def keys(self):
        return self.pedigrees.keys()

    @property
    def individuals(self):
        inds = []
        for pedigree in sorted(self, key=lambda x: x.label):
            inds.extend(sorted((x for x in pedigree), key=lambda x: x.id))
        return inds

    def _getindividual(self, label):
        for x in self.individuals:
            if x.id == label: return x
        raise KeyError('Individual not in collection')

    def phenotypes(self):
        """ Returns the available phenotypes for analysis """
        return set(reduce(add, [x.phenotypes.keys() for x in
                                self.individuals]))
    
    @property
    def chromosomes(self):
        k = self.pedigrees.keys()[0]
        return self.pedigrees[k].chromosomes

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

    def merge(self, pop):
        for ped in self:
            ped.merge(pop)

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
