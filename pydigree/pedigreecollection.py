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
    def __iter__(self): return (x for x in self.pedigrees.values())
    def __getitem__(self,key): return self.pedigrees[key]
    def __contains__(self,item): return item in self.pedigrees.values()
    def __len__(self): return len(self.pedigrees)
    def __setitem__(self,key,value): self.pedigrees[key] = value
    def __delitem__(self,key): del self.pedigrees[key]
    def keys(self): return self.pedigrees.keys()
    def individuals(self):
        inds = []
        for pedigree in sorted(self,key=lambda x:x.label):
            inds.extend(sorted((x for x in pedigree),key=lambda x:x.id))
        return inds
    def phenotypes(self):
        """ Returns the available phenotypes for analysis """
        return set(reduce(add,[x.phenotypes.keys() for x in self.individuals()]))
    ### Matrix functions
    ###
    def additive_relationship_matrix(self):
        """
        Returns a block diagonal matrix of additive relationships for each pedigree.
        See notes on Pedigree.additive_relationship_matrix
        """
        return block_diag([x.additive_relationship_matrix() for x in self],format='bsr')
    def dominance_relationship_matrix(self):
        """
        Returns a block diagonal matrix of dominance relationships for each pedigree.
        See notes on Pedigree.dominance_relationship_matrix
        """
        return block_diag([x.dominance_relationship_matrix() for x in self],format='bsr')
    def mitochondrial_relationship_matrix(self):
        """
        Returns a block diagonal matrix of mitochondrial relationships for each pedigree.
        See notes on Pedigree.mitochondrial_relationship_matrix
        """
        return block_diag([x.mitochondrial_relationship_matrix() for x in self],format='bsr')
    
