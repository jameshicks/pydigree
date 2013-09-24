#!/usr/bin/env python

# Packages we're going to use
import itertools
import random
import math
import numpy as np

# Abstract base class for population and pedigree
from collections import MutableMapping

# Other pydigree objects
from individual import Individual
from chromosome import Chromosome
from paths import kinship, fraternity
from common import *
from misc import is_missing_genotype
from recombination import recombine

# C extension functions
from pydigree._pydigree import sample_with_replacement,random_pairs
from pydigree._pydigree import choice_with_probs


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
    def remove_individual(self,ind):
        del self[ind.id]
    def random_identifier(self):
        idx = 0
        while True:
            id = str(idx)
            if self.name: id = name + id
            if id not in self.population.keys():
                yield id
            idx+1
    def males(self):
        """ Returns list of males in population """
        return [x for x in self if x.sex == 0]
    def females(self):
        """ Returns list of females in population """
        return [x for x in self if x.sex == 1]
    ### Chromosome functions
    ###
    ###
    def add_chromosome(self,chrom): self.chromosomes.append(chrom)
    def chromosome_count(self): return len(self.chromosomes)
    
    ### Chromosome pool functions
    ###
    ###
    def chrom_pool_size(self):
        """ Returns the size of the pool of available chromosomes"
        return len(self.pool[0])
    def initialize_pool(self):
        """ Initializes a pool of chromosomes for simulation """
        self.pool = [None]*self.chromosome_count()
        for i,q in enumerate(self.chromosomes):
            self.pool[i]=[q.linkageequilibrium_chromosome() for x in xrange(2*self.n0)]
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
            np = [choose_chrom(self.pool[i],c.genetic_map) for x in xrange(gensize)]
            self.pool[i] = np
    ### Random mating
    ###
    ###
    def mate(self,ind1,ind2,id):
        """
        Creates an individual as the child of two specificied individual objects
        randomly chooses a sex. 
        """
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
    ### Genotype Functions
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
        """ Clears genotypes for every one in population """
        for x in self: x.clear_genotypes()
    ### Frequency functions
    ### 
    def alleles(self,location, constraint=None):
        """
        The list of available alleles at a location in this population

        The argument constraint is a function that acts on an individual. If
        constraint(individual) can evaluate as True that person is included
        """
        if not constraint:
            constraint = lambda x: True
        gen = (x for x in self if constraint(x))
        alleles = list(flatten(x.get_genotype(location) for x in gen if x.has_genotypes()))
        return [x for x in alleles if x != 0]
    def allele_frequency(self,location,allele,constraint=None):
        """
        Returns the frequency (as a percentage) of an allele in this population

        The argument constraint is a function that acts on an individual. If
        constraint(individual) can evaluate as True that person is included
        """
        alleles = self.alleles(location,constraint=constraint)
        freqtab = table(alleles)
        if allele not in freqtab: return 0
        return freqtab[allele] / float(len(alleles))
    def major_allele(self,location,constraint=None):
        alleles = self.alleles(location,constraint=constraint)
        freqtab = table(alleles)
        return sorted(freqtab.keys(), key=lambda x:freqtab[x])[0]
    def ld(self,locusa,locusb,allelea=None,alleleb=None, method='r2'):
        """
        Returns a measure of linkage disquilibirum between alleles at two loci,
        for methods in {D, Dprime, r2}. If alleles aren't specified, the
        major alleles at each locus are chosen.
        
        LD metrics are functions of D, the measure of deviance from independence in
        sampling alleles. D = x11 - p1q1, where x11 is the frequency of allele 1 at
        locus p occuring with allele 1 at locus q.
        
        Lewontin's D' is a normalization of D, since D is dependent on frequency.
        D' = D / Dmax, where:
        Dmax = min(p1*q1, p2*q2) if D > 0
        Dmax = min(p2*q1,p1*q2) if D < 0
        
        r2 is a third metric, equal to
        D**2 / (p1 * p2 * q1 * q2)
        r2 ranges from 0 to 1
        
        Alleles are in linkage equilibrium if any of these metrics are 0. 
        """
        method = method.lower()
        if method not in set(['r2','dprime','d']):
            raise ValueError("LD metric must be one of r2,Dprime or D") 
        if not allelea: allelea = self.major_allele(locusa)
        if not alleleb: alleleb = self.major_allele(locusb)
        pop = [x for x in self if \
               (not is_missing_genotype(x.get_genotype(locusa))) and \
               (not is_missing_genotype(x.get_genotype(locusb)))]
        # Allele 1
        p1 = self.allele_frequency(locusa,allelea)
        p2 = 1 - p1
        # Allele 2
        q1 = self.allele_frequency(locusa,allelea)
        q2 = 1 - q1
        x = sum(1 for x in pop \
                if x.has_allele(locusa,allelea) and x.has_allele(locusb,alleleb))
        D = x - pa*qa
        # All of these metrics are normalizations of D by dividing it by
        # some expression. If D==0 all of them will equal 0.
        if D == 0: return 0 
        elif method == 'd': return D
        elif method == 'dprime':
            if D > 0: Dmax = min(p1*q1,p2*q2)
            else: Dmax = min(p2*q1,p1*q2)
            return D / Dmax
        elif method == 'r2':
            return (D**2) / (p1 * p2 * q1 * q2)
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
        self.fratmat = {}
    def __prepare_nonfounder_contraint(self,con):
        if not con: return lambda x: x.is_founder()
        else: return lambda x: x.is_founder() and con(x)
    def bit_size(self):
        """
        Returns the bit size of the pedigree. The bitsize is defined as 2*n-f
        where n is the number of nonfounders and f is the number of founders.

        This represents the number of bits it takes to represent the inheritance
        vector in the Lander-Green algorithm. 
        """
        t = table([x.is_founder() for x in self])
        return 2 * t[False] - t[True]
    ### Frequency
    def alleles(self,location, constraint=None, nonfounders=False):
        """
        Like Population.alleles, except constrained to founder individuals

        If nonfounders is True, it's just a call to Population.alleles.
        """
        if nonfounders:
            return Population.alleles(self,location,constraint)
        con = self.__prepare_nonfounder_constraint(constraint)
        return Population.alleles(self,location,con)
    def allele_frequency(self,location,allele,constraint=None, nonfounders=False):
        """
        Like Population.alleles, except constrained to founder individuals.
        If nonfounders is True, it's just a call to Population.alleles.        
        """
        if nonfounders: return Population.allele_frequency(self,location,allele,constraint)
        constraint = self.__prepare_nonfounder_constraint(constraint)
        return Population.allele_frequency(self,location,allele,constraint=constraint)
    ### Relationships
    ###
    def kinship(self,id1,id2):
        """
        Get the Malecot coefficient of coancestry for two individuals in the pedigree.
        (See notes in pydigree.paths.kinship). For pedigree objects, results are stored
        to reduce the calculation time for kinship matrices.

        This is a convenience wrapper for paths.kinship, which takes pedigree objects as
        arguments. This function takes id labels and looks them up in the pedigree, and
        calls paths.kinship on those individual objects. 
        """
        pair = frozenset([id1,id2])
        if pair not in self.kinmat:
            k = kinship(self[id1],self[id2])
            self.kinmat[pair] = k
            return k
        else: return self.kinmat[pair]
    def fraternity(self,id1,id2):
        """
        Like Pedigree.kinship, this is a convenience function for getting fraternity
        coefficients for two pedigree memebers by their ID label. This is a wrapper
        for paths.fraternity
        """
        pair = frozenset([id1,id2])
        if pair not in self.fratmat:
            f = fraternity(self[id1],self[id2])
            self.fratmat[pair] = f
            return f
        else: return self.fratmat[pair]
    def inbreeding(self,id):
        """
        Like Pedigree.kinship, this is a convenience function for getting inbreeding
        coefficients for individuals in pedigrees by their id label. As inbreeding
        coefficients are the kinship coefficient of the parents, this function calls
        Pedigree.kinship to check for stored values.
        """
        ind = self[id]
        return self.kinship(ind.father.id,ind.mother.id)
    def makeA(self):
        """
        Calculates an additive relationship matrix (the A matrix) for quantiatitive genetics.
        A_ij = 2 * kinship(i,j) if i != j. (See the notes on function 'kinship')
        A_ij = 1 + inbreeding(i) if i == j (inbreeding(i) is equivalent to kinship(i.father,i.mother))

        Important: the rows/columns are sorted on ids. If you're not sure about this, try
        sorted(x.id for x in ped) to see the ordering.
        """
        ids = sorted(x.id for x in self)
        mat = []
        for a in ids:
            row = []
            for b in ids:
                if a == b: row.append(1 + self.inbreeding(a))
                else: row.append(2 * self.kinship(a,b))
            mat.append(row)
        return np.array(mat)
    def makeD(self):
        """
        Calculates the dominance genetic relationship matrix (the D matrix) for quantitative genetics.
        D_ij = .25 * fraternity(i,j) if i != j (See notes on function 'fraternity')
        D_ij = 1 if i == j.

        Important: the rows/columns are sorted on ids. If you're not sure about this, try
        sorted(x.id for x in ped) to see the ordering.
        """
        ids = sorted(x.id for x in self)
        mat = []
        for a in ids:
            row = []
            for b in ids:
                if a == b: row.append(1)
                else: row.append(.25 * self.fraternity(a,b))
            mat.append(row)
        return np.array(mat)
    ### Gene dropping
    ###
    def simulate_ibd_states(self):
        """
        Simulate IBD patterns by gene dropping: Everyone's genotypes reflect the
        founder chromosome that they received the genotype from. You can then use
        misc.ibs to determine IBD state. This effectively an infinite-alleles simulation.
        """
        for x in self:
            if x.is_founder(): x.label_genotypes()
        for x in self:
            if x.is_founder(): continue
            x.get_genotypes()
