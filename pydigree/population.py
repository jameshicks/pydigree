#!/usr/bin/env python

import itertools
import random
import math
import numpy as np
from collections import MutableMapping

from individual import Individual
from chromosome import Chromosome
from paths import kinship
from common import *

from pydigree._pydigree import sample_with_replacement,random_pairs
from pydigree._pydigree import choice_with_probs
from recombination import recombine

### Population growth models.
# These are easy enough to supply your own if you wanted
def exponential_growth(p,r,t):
    """
    Models exponential growth over discrete generations.
    p: initial population
    r: growth rate
    t: number of generations
    """
    return p*math.exp(r*t)
def logistic_growth(p,r,k,t):
    """
    Models logistic growth over discrete generations.
    p: initial population
    r: growth rate
    k: final population
    t: number of generations
    """
    return (p*k) / ( p + (k-p) * math.exp(-r*t) )
                

# Classes
class Population(MutableMapping):
    ### Methods for mapping types
    def __init__(self,intial_pop_size=0,name=None):
        self.chromosomes = []
        self.pool = []
        self.population = {}
        self.n0 = intial_pop_size
        self.name = name

    ### Things I have to implement for the ABC
    ###
    ###
    def __iter__(self): return iter(self.population.values())
    def __getitem__(self,key): return self.population[key]
    def __contains__(self,item): return item in self.population.values()
    def __len__(self): return len(self.population)
    def __setitem__(self,key,value): self.population[key] = value
    def __delitem__(self,key): del self.population[key]

    def size(self): return len(self.population)
    def remove_ancestry(self):
        for x in self: x.remove_ancestry()

    ### Adding and removing people
    ###
    ###
    def register_individual(self,ind):
        if ind.id in self.population:
            raise ValueError('ID %s already in population!' % ind.id)
        self.population[ind.id] = ind
    def random_identifier(self):
        idx = 0
        while True:
            id = str(idx)
            if self.name: id = name + id
            if id not in self.population.keys():
                yield id
            idx+1
    ### Chromosome functions
    ###
    ###
    def add_chromosome(self,chrom): self.chromosomes.append(chrom)
    def chromosome_count(self): return len(self.chromosomes)
    
    ### Chromosome pool functions
    ###
    ###
    def chrom_pool_size(self): return len(self.pool[0])
    def initialize_pool(self):
        self.pool = [None]*self.chromosome_count()
        for i,q in enumerate(self.chromosomes):
            self.pool[i]=[q.linkageequilibrium_chromosome() for x in range(2*self.n0)]
    def iterate_pool(self,gensize):
        """
        Iterate pool simulates a generation of random mating
        between chromosomes instead of individuals. The pool of
        population chromosomes then contains chromosomes from the
        new generation. 

        """
        gensize = int(gensize)
        for i,c in enumerate(self.chromosomes):
            # Chromosomes have a 1/2 chance of being recombined with another
            def choose_chrom(pool,map):
                q,w = random.choice(pool),random.choice(pool)
                return recombine(q,w,map)
            np = [choose_chrom(self.pool[i],c.genetic_map) for x in range(gensize)]
            self.pool[i] = np

    ### Random mating
    ###
    ###
    def mate(self,ind1,ind2,id):
        child = Individual(self,id,ind1,ind2,random.choice([0,1]))
        return child
    def random_mating_generation(self,gensize):
        def rand_mate(pop,name=None):
            name = str(name) if not self.name else self.name + str(name)
            return self.mate(random.choice(self.males()), random.choice(self.females()),name)
        newpop = {}
        for i in xrange(gensize):
            newind = rand_mate(self,name=i)
            newpop[newind.id] = newind
        self.population = newpop

    def _founder_ind(self):
        return Individual(self,self.size(),None,None,random.choice([0,1]))
        
    def add_founder_individual(self):
        i = self._founder_ind()
        self.register_individual(i)
    ### Genotype/Allele Functions
    ###
    ###
    def get_founder_genotypes(self):
        g = []
        for i,x in enumerate(self.pool):
            g.append( [random.choice(x),random.choice(x)] )
        return g
    def get_genotypes(self):
        for x in self: x.get_genotypes()
    def get_linkage_equilibrium_genotypes(self):
        return [ [c.linkage_equilibrium_chromosome(),c.linkage_equilibrium_chromosome()]
                for c in self.chromosomes]
    def clear_genotypes(self):
        for x in self: x.clear_genotypes()
    def alleles(self,location):
        """ The list of available alleles at a location in this population """
        return list(flatten(x.get_genotype(location) for x in self if x.has_genotypes()))
    def allele_frequency(self,location,allele):
        """ Returns the frequency (as a percentage) of an allele in this population """
        alleles = self.alleles(location)
        freqtab = table(alleles)
        if allele not in freqtab: return 0
        return freqtab[allele] / float(len(alleles))
    def males(self):
        """ Returns list of males in population """
        return [x for x in self if x.sex == 0]
    def females(self):
        """ Returns list of females in population """
        return [x for x in self if x.sex == 1]

    ### Phenotype functions
    ###
    ###
    def predict_phenotype(self,trait):
        """
        Shortcut function to call predict_phenotype(trait) for every individual
        in the population.
        """
        for x in self: x.predict_phenotype(trait)

class Pedigree(Population):
    def __init__(self):
        Population.__init__(self)
        self.kinmat = {}
    def bit_size(self):
        t = table([x.is_founder() for x in self])
        return 2 * t[False] - t[True]
    def kinship(self,id1,id2):
        """
        Get the Malecot coefficient of coancestry for two individuals in the pedigree.
        (See notes in pydigree.paths.kinship). For pedigree objects, results are stored
        to reduce the calculation time for kinship matrices.

        This is a convenience wrapper for paths.kinship, which takes pedigree objects as
        arguments. This function takes id labels and looks them up in the pedigree, and
        calls paths.kinship on those pedigree objects. 
        """
        pair = frozenset([id1,id2])
        if pair not in self.kinmat:
            k = kinship(self[id1],self[id2])
            self.kinmat[pair] = k
            return k
        else: return self.kinmat[pair]
    def inbreeding(self,id):
        """
        Like Pedigree.kinship, this is a convenience function for getting inbreeding
        coefficients for individuals in pedigrees by their id label. As inbreeding
        coefficients are the kinship coefficient of the parents, this function calls
        Pedigree.kinship to check for stored values.
        """
        ind = self[id]
        return self.kinship(ind.father.id,ind.mother.id)
    def simulate_ibd_states(self):
        """
        Simulate IBD patterns by gene dropping: Everyone's genotypes reflect the
        founder that they received the genotype from.  
        """
        for x in self:
            if x.is_founder(): x.label_genotypes()
        for x in self:
            if x.is_founder(): continue
            x.get_genotypes()
    def makeA(self):
        """
        Calculates an additive relationship matrix (the A matrix) for quantiatitive genetics.
        A_ij = 2 * kinship(i,j) if i != j. (See the notes on function 'kinship')
        A_ij = 1 + inbreeding(i) if i == j (inbreeding(i) is equivalent to kinship(i.father,i.mother))

        Important: the rows/columns are sorted on ids. If you're not sure about this, try
        sorted(x.id for x in ped) to see the ordering.
        """
        inds = sorted(x.id for x in self)
        mat = []
        for a in ids:
            row = []
            for b in ids:
                if a == b: row.append(1 + self.inbreeding(a))
                else: row.append(2 * self.kinship(a,b))
            mat.append(row)
        return np.array(mat)
