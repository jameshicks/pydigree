#!/usr/bin/env python
import itertools
import random
from recombination import recombine
from paths import *


class Individual():
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
    def get_genotypes(self,linkeq=False):
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
        if not self.has_genotypes(): raise ValueError('Individual has no genotypes!')
        chr,pos = location
        return tuple([x[pos] for x in self.genotypes[chr]])
    def label_genotypes(self):
        g = []
        for x in self.population.chromosomes:
            g.append( [ ['%sP' % self.id ] * len(x.genetic_map), ['%sM' % self.id] * len(x.genetic_map) ]  )
        self.genotypes = g
    def clear_genotypes(self): self.genotypes = None
    def clear_phenotypes(self): self.phenotypes = {}
    ### Functions about ancestry
    ###
    def is_founder(self):
        return self.father is None and self.mother is None
    def ancestors(self):
        """ Recursively searches for ancestors. """
        if self.is_founder(): return set()
        return set([self.father,self.mother]) | self.father.ancestors() | self.mother.ancestors()
    def descendants(self):
        """ Recursively searches for descendants. """
        return set(self.children + list(flatten([x.descendants() for x in self.children])))
    def siblings(self): return set(self.father.children) & set(self.mother.children)
    def remove_ancestry(self):
        """
        Removes ancestry: makes a person a founder. Cannot be used on an individual in a
        pedigree, because the pedigree structure is already set.
        """
        if self.pedigree: raise ValueError('Individual in a pedigree!')
        self.father = None
        self.mother = None
    def inbreeding(self):
        if self.is_founder(): return 0.0
        return kinship(self.father,self.mother)
    ### Functions for breeding
    ###
    def gamete(self):
        if not self.genotypes: self.get_genotypes()
        def randbit(): random.choice([True,False])
        g = []
        for i,chr in enumerate(self.genotypes):
            ### FIXME: _pydigree.recombine should work on sequences
            chr1,chr2 = chr
            if randbit():
                chrom = recombine(chr1,chr2,self.population.chromosomes[i].genetic_map)
            else:
                chrom = recombine(chr2,chr1,self.population.chromosomes[i].genetic_map)
            g.append(chrom)
        return g
    def fertilize(self,father,mother):
        return [ [x,y] for x,y in itertools.izip(father,mother) ]
    def pprint_genos(self):
        for i,c in enumerate(self.genotypes):
            print 'Chromosome %s' % i 
            for g in zip(*c):
                print '\t'.join([str(x) for x in g])
    def predict_phenotype(self,trait):
        self.phenotypes[trait.name] = trait.predict_phenotype(self)

    
