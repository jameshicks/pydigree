#!/usr/bin/env python

from common import *
from operator import add
from collections import MutableMapping

from scipy.sparse import block_diag

from pydigree.pedigree import Pedigree


class PedigreeCollection(MutableMapping):

    def __init__(self):
        self.pedigrees = {}

    # Things I have to implement for the ABC
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

    def add_pedigree(self, ped):
        if not isinstance(ped, Pedigree):
            raise ValueError('{} not of type Pedigree')
        elif ped.label in self.keys():
            raise ValueError('A pedigree labeled {} already in collection'.format(ped.label))
        else:
            self[ped.label] = ped

    @property
    def individuals(self):
        '''
        Returns a list of the individuals represented by all pedigrees, 
        sorted by pedigree label, id label 
        '''
        inds = []
        for pedigree in sorted(self, key=lambda x: x.label):
            inds.extend(sorted((x for x in pedigree), key=lambda x: x.label))
        return inds

    def founders(self):
        ''' Returns a list of founder individuals across all pedigrees '''
        return [x for x in self.individuals if x.is_founder()]

    def nonfounders(self):
        ''' Returns a list of founder individuals across all pedigrees '''
        return [x for x in self.individuals if not x.is_founder()]
    
    def _getindividual(self, label):
        for x in self.individuals:
            if x.label == label:
                return x
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
        for x in self.individuals:
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

    # Matrix functions
    ###
    def additive_relationship_matrix(self, ids=None):
        """
        Returns a block diagonal matrix of additive relationships
        for each pedigree.

        See notes on Pedigree.additive_relationship_matrix
        """
        mats = [x.additive_relationship_matrix(ids) for x in
                sorted(self, key=lambda x: x.label)]
        mats = [x for x in mats if x.size > 0]
        return block_diag(mats, format='bsr')

    def dominance_relationship_matrix(self, ids=None):
        """
        Returns a block diagonal matrix of dominance relationships
        for each pedigree.

        See notes on Pedigree.dominance_relationship_matrix
        """
        mats = [x.dominance_relationship_matrix(ids) for x in
                sorted(self, key=lambda x: x.label)]
        mats = [x for x in mats if x.size > 0]
        return block_diag(mats, format='bsr')

    def mitochondrial_relationship_matrix(self, ids=None):
        """
        Returns a block diagonal matrix of mitochondrial relationships
        for each pedigree.

        See notes on Pedigree.mitochondrial_relationship_matrix
        """
        mats = [x.mitochondrial_relationship_matrix(ids) for x in
                sorted(self, key=lambda x: x.label)]
        return block_diag(mats, format='bsr')
