#!/usr/bin/env python

import itertools
import random
from recombination import recombine
from paths import *

def is_missing_genotype(g): g == (0,0)

class Individual(object):
    def __init__(self,population,id,father,mother,sex):
        # Every individual is observed within a population with certain genotypes
        #  available. This makes recombination book-keeping easier. 
        self.population = population

        if id: self.id = id
        elif id is None and self.population:
            self.id = self.population.random_identifier().next()
        else: self.id = None
        
        self.father = father
        self.mother = mother
        self.sex = sex # 0:M 1:F
        self.pedigree = None
        self.genotypes = None
        self.observed_genos = True
        self.phenotypes = {}
        self.children = []
        if isinstance(self.father,Individual): self.father.register_child(self)
        if isinstance(self.mother,Individual): self.mother.register_child(self)
    def __str__(self):
        if self.is_founder(): return 'Individual %s (FOUNDER)' % self.id
        return 'Individual %s (F:%s,M:%s)' % (self.id, self.father.id,self.mother.id)

    def register_child(self,child): self.children.append(child)
    def register_with_parents(self):
        if self.is_founder(): return
        self.father.register_child(self)
        self.mother.register_child(self)

    ### Functions about genotypes
    ###
    def has_genotypes(self): return self.genotypes is not None
    def __fail_on_observed_genos(self):
        if self.observed_genos:
            raise ValueError('Individual has observed genotypes')        
    def get_genotypes(self,linkeq=False):
        self.__fail_on_observed_genos()
        if self.has_genotypes(): return 
        if self.is_founder() and self.population is None:
            raise ValueError('Founder individual %s has no assigned population!' % (self.id)) 
        if self.is_founder():
            if linkeq:
                self.genotypes = self.population.get_linkage_equilibrium_genotypes()
            else:
                self.genotypes = self.population.get_founder_genotypes()
        else:
            self.genotypes = self.fertilize(self.father.gamete(),self.mother.gamete())
    def get_genotype(self,location):
        """
        Returns a tuple of alleles, in the format (paternal allele, maternal allele)
        """
        self.__fail_on_observed_genos()
        if not self.has_genotypes():
            raise ValueError('Individual has no genotypes!')
        chr,pos = location
        return tuple([x[pos] for x in self.genotypes[chr]])
    def has_allele(self,location,allele):
        """ Returns true if individual has the specified allele at location """
        g = self.get_genotype(location)
        if is_missing_genotype(g): return None
        else: return allele in g
    def label_genotypes(self):
        """
        Gives the individual label genotypes. When these genotypes are transmitted on to the next
        generation, you can see where each allele in the next generation came from. Useful for
        gene dropping simulations, probably not much else.
        """
        self.__fail_on_observed_genos()
        g = []
        for x in self.population.chromosomes:
            g.append( [ ['%sP' % self.id ] * len(x.genetic_map), ['%sM' % self.id] * len(x.genetic_map) ]  )
        self.genotypes = g
    def clear_genotypes(self):
        """ Removes genotypes """
        self.observed_genos = False
        self.genotypes = None
    ### Functions about ancestry and family
    ###
    def is_founder(self):
        """ Returns true if individual is a founder """
        return self.father is None and self.mother is None
    def ancestors(self):
        """ Recursively searches for ancestors. """
        if self.is_founder(): return set()
        return set([self.father,self.mother]) | self.father.ancestors() | self.mother.ancestors()
    def descendants(self):
        """ Recursively searches for descendants. """
        return set(self.children + list(flatten([x.descendants() for x in self.children])))
    def siblings(self, include_halfsibs=False):
        """
        Returns this individuals sibliings.
        If include halfsibs is set, half sibs will also be included
        """
        if include_halfsibs:
            return set(self.father.children) | set(self.mother.children)
        else:
            return set(self.father.children) & set(self.mother.children)
    def matriline(self):
        """
        Returns a label by recursively searching for the individual's mother's mother's
        mother's etc. until it reachesa founder mother, in which case it returns that
        ancestor's id.

        Useful reference:
        Courtenay et al. 'Mitochondrial haplogroup X is associated with successful aging in the Amish.'
        Human Genetics (2012). 131(2):201-8. doi: 10.1007/s00439-011-1060-3. Epub 2011 Jul 13.
        """
        if self.is_founder(): return self.id
        else: return self.mother.matriline()
    def remove_ancestry(self):
        """
        Removes ancestry: makes a person a founder. Cannot be used on an individual in a
        pedigree, because the pedigree structure is already set.
        """
        if self.pedigree: raise ValueError('Individual in a pedigree!')
        self.father = None
        self.mother = None
    def inbreeding(self):
        """
        Returns the inbreeding coefficient (F) for the individual.
        """
        # Two edge cases where inbreedings must be 0
        if self.is_founder(): return 0.0
        if self.father.is_founder() or self.mother.is_founder(): return 0.0
        return kinship(self.father,self.mother)
    ### Functions for breeding
    ###
    def gamete(self):
        if not self.genotypes: self.get_genotypes()
        def randbit(): random.choice([True,False])
        g = []
        for i,chr in enumerate(self.genotypes):
            chr1,chr2 = chr
            if randbit():
                chrom = recombine(chr1,chr2,self.population.chromosomes[i].genetic_map)
            else:
                chrom = recombine(chr2,chr1,self.population.chromosomes[i].genetic_map)
            g.append(chrom)
        return g
    def fertilize(self,father,mother):
        return [ [x,y] for x,y in itertools.izip(father,mother) ]
    ### Phenotype functions
    def predict_phenotype(self,trait):
        self.phenotypes[trait.name] = trait.predict_phenotype(self)
    def clear_phenotypes(self):
        """ Removes phenotypes """
        self.phenotypes = {}
    
