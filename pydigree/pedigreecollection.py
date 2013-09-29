#!/usr/bin/env python

from pedigree import Pedigree
from scipy.linalg import block_diag
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
    ### Matrix functions
    ###
    def additive_relationship_matrix(self):
        """
        Returns a block diagonal matrix of additive relationships for each pedigree.
        See notes on Pedigree.additive_relationship_matrix
        """
        return block_diag(x.additive_relationship_matrix() for x in self)
    def dominance_relationship_matrix(self):
        """
        Returns a block diagonal matrix of dominance relationships for each pedigree.
        See notes on Pedigree.dominance_relationship_matrix
        """
        return block_diag(x.dominance_relationship_matrix() for x in self)
    def mitochondrial_relationship_matrix(self):
        """
        Returns a block diagonal matrix of mitochondrial relationships for each pedigree.
        See notes on Pedigree.mitochondrial_relationship_matrix
        """
        return block_diag(x.mitochondrial_relationship_matrix() for x in self)
    
