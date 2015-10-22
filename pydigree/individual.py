#!/usr/bin/env python

import itertools
from itertools import izip

import numpy as np
import pandas as pd

from pydigree.recombination import recombine
from pydigree.paths import *
from pydigree.common import *
from pydigree.genotypes import Alleles, LabelledAlleles
from pydigree.exceptions import IterationError

missing_genotype = (0, 0)


def is_missing_genotype(g):
    return g == missing_genotype


class Individual(object):

    def __init__(self, population, label, father=None, mother=None, sex=None):
        # Every individual is observed within a population with certain
        # genotypes available. This makes recombination book-keeping easier.
        self.population = population
        if id:
            self.label = label
        else:
            self.label = None
        self.father = father
        self.mother = mother
        self.sex = sex  # 0:M 1:F
        self.pedigree = None
        self.genotypes = None
        self.observed_genos = False
        self.phenotypes = {}
        self.attrib = {}
        self.children = []
        if isinstance(self.father, Individual):
            self.father.register_child(self)
        if isinstance(self.mother, Individual):
            self.mother.register_child(self)

    def __str__(self):
        try:
            pop, lab = self.full_label
            if self.is_founder():
                return 'Individual {}:{} (FOUNDER)'.format(pop, lab)
            else:
                return 'Individual {}:{} (F:{},M:{})'.format(pop, lab,
                                                             self.father.label,
                                                             self.mother.label)
        except AttributeError:
            return 'Individual %s (Unlinked)' % self.label

    def __repr__(self):
        return self.__str__()

    def register_child(self, child):
        self.children.append(child)

    def register_with_parents(self):
        if self.is_founder():
            return
        self.father.register_child(self)
        self.mother.register_child(self)

    @property
    def full_label(self):
        ''' Returns a 2-tuple of pedigree label and individual label '''
        return (self.population.label, self.label)

    # Functions about genotypes
    #
    @property
    def chromosomes(self):
        ''' Returns a list of the individuals ChromosomeTemplate objects '''
        if self.pedigree is not None:
            return self.pedigree.chromosomes
        else:
            return self.population.chromosomes

    def _init_genotypes(self, blankchroms=True, dtype=None):
        """ 
        Initializes genotypes so that all genotypes are missing if blankchroms 
        is true, otherwise, just sets what would be the chromosome to None
        """
        if blankchroms:
            self.genotypes = [(chrom.empty_chromosome(dtype=dtype),
                               chrom.empty_chromosome(dtype=dtype))
                              for chrom in self.chromosomes]
        else:
            self.genotypes = [[None, None] for chrom in self.chromosomes]

    def has_genotypes(self):
        """ Returns True if an individual has genotypes """
        return self.genotypes is not None

    def __fail_on_observed_genos(self):
        if self.observed_genos:
            raise ValueError('Individual has observed genotypes')

    def __fail_on_non_genotyped(self):
        if not self.has_genotypes():
            raise ValueError('Individual has no genotypes')

    def get_genotypes(self, linkeq=False):
        """
        Retrieve genotypes from a chromosome pool if present, or else a
        chromosome generated under linkage equilibrium
        """
        self.__fail_on_observed_genos()
        if not self.population.pool:
            linkeq = True
        if self.has_genotypes():
            return
        if self.is_founder() and self.population is None:
            raise ValueError('Founder ind %s has no population!'
                             % (self.label))
        if self.is_founder():
            pop = self.population
            if linkeq:
                self.genotypes = pop.get_linkage_equilibrium_genotypes()
            else:
                self.genotypes = pop.get_founder_genotypes()
        else:
            self.genotypes = Individual.fertilize(self.father.gamete(),
                                                  self.mother.gamete())

    def get_genotype(self, loc, checkhasgeno=True):
        """
        Returns a tuple of alleles, in the format:
        paternal allele, maternal allele
        """
        if checkhasgeno and not self.has_genotypes():
            raise ValueError('Individual has no genotypes!')

        return (self.genotypes[loc[0]][0][loc[1]],
                self.genotypes[loc[0]][1][loc[1]])

    def get_constrained_genotypes(self, constraints, linkeq=True):
        '''
        Gets genotypes from parents (or population if the individual is a
        founder) subject to constraints. Used by the simulation objects
        in pydigree.simulation
        '''
        self.get_genotypes(linkeq=linkeq)
        for constraint in constraints:
            locus, chromatid, allele, method = constraint
            gt = list(self.get_genotype(locus))
            gt[chromatid] = allele
            self.set_genotype(locus, tuple(gt))

    def set_genotype(self, location, genotype):
        """ Manually set a genotype """
        self.__fail_on_non_genotyped()
        self.genotypes[location[0]][0][location[1]] = genotype[0]
        self.genotypes[location[0]][1][location[1]] = genotype[1]

    def _set_genotypes(self, gts):
        self.genotypes = gts

    def update(self, other):
        '''
        Takes another individual object, merges/updates phenotypes with the other
        individual object and REPLACES self's genotypes with other's
        '''
        self.phenotypes.update(other.phenotypes)
        self.genotypes = other.genotypes

    def has_allele(self, location, allele):
        """
        Returns True if individual has the specified allele at location
        """
        g = self.get_genotype(location)
        if is_missing_genotype(g):
            return None
        else:
            return allele in g

    def label_genotypes(self):
        """
        Gives the individual label genotypes. When these genotypes are
        transmitted on to the next generation, you can see where each allele
        in the next generation came from. Useful for gene dropping simulations,
        probably not much else.
        """
        self.__fail_on_observed_genos()

        def labelled_chromatids(i, c):
            a = LabelledAlleles.founder_chromosome(self, i, 0, chromobj=c)
            b = LabelledAlleles.founder_chromosome(self, i, 1, chromobj=c)
            return (a, b)

        g = [labelled_chromatids(i, c) for i, c in enumerate(self.chromosomes)]
        self.genotypes = g

    def delabel_genotypes(self):
        '''
        When an individual has label genotypes, replaces the labels with 
        the ancestral allele corresponding to the label
        '''
        for chromoidx, chromosome in enumerate(self.genotypes):
            for chromaidx, chromatid in enumerate(chromosome):
                newchromatid = chromatid.delabel()
                self.genotypes[chromoidx][chromaidx] = newchromatid

    def clear_genotypes(self):
        """ Removes genotypes """
        self.observed_genos = False
        self.genotypes = None

    # Functions about ancestry and family
    #
    def is_founder(self):
        """ Returns true if individual is a founder """
        return self.father is None and self.mother is None

    def is_marryin_founder(self):
        """
        Returns true if an individual is a marry-in founder
        i.e.: the individual is a founder (depth: 0) and has a child with 
        depth > 1
        """
        if not self.is_founder():
            return False
        return any(x.depth > 1 for x in self.children)

    def parents(self):
        ''' Returns the individual's father and mother in a 2-tuple '''
        return self.father, self.mother

    def ancestors(self):
        """
        Recursively searches for ancestors.

        Returns: A set object with Individual objects for all the ancestors
        of this Individual object
        """
        if self.is_founder():
            return set()
        return set([self.father, self.mother]) | \
            self.father.ancestors() | \
            self.mother.ancestors()

    def descendants(self):
        """
        Recursively searches for descendants.

        Returns: a set of individual objects for all the descendants
        of this individual
        """
        return set(self.children + list(flatten([x.descendants()
                                                 for x in self.children])))

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
        Returns a label by recursively searching for the individual's mother's
        mother's mother's etc. until it reachesa founder mother, in which case
        it returns that ancestor's id.

        Useful reference:
        Courtenay et al. 'Mitochondrial haplogroup X is associated with
        successful aging in the Amish.' Human Genetics (2012). 131(2):201-8.
        doi: 10.1007/s00439-011-1060-3. Epub 2011 Jul 13.
        """
        if self.is_founder():
            return self.label
        else:
            return self.mother.matriline()

    def patriline(self):
        """
        Returns a label by recursively searching for the individual's mother's
        father's father's etc. until it reaches a founder father, in which case
        it returns that ancestor's id.

        Analagous to individual.matriline.
        """
        if self.is_founder():
            return self.label
        else:
            return self.father.patriline()

    @property
    def depth(self):
        """
        Returns the depth of an individual in the pedigree, a rough measure of
        what generation in the pedigree the individual is. Defined as:
        depth = 0 if individual is a founder, else the maximum of the
        depth of each parent

        Arguements: None
        Returns: An integer
        """
        if self.is_founder():
            return 0
        elif 'depth' in self.attrib:
            return self.attrib['depth']
        else:
            d = 1 + max(self.father.depth, self.mother.depth)
            self.attrib['depth'] = d
            return d

    def remove_ancestry(self):
        """
        Removes ancestry: makes a person a founder. Cannot be used on an
        individual in a pedigree, because the pedigree structure is
        already set.
        """
        if self.pedigree:
            raise ValueError('Individual in a pedigree!')
        self.father = None
        self.mother = None

    def inbreeding(self):
        """
        Returns the inbreeding coefficient (F) for the individual.
        """
        # Two edge cases where inbreedings must be 0
        if self.is_founder():
            return 0.0
        if self.father.is_founder() or self.mother.is_founder():
            return 0.0
        if 'inbreed' in self.attrib:
            return self.attrib['inbreed']
        else:
            self.attrib['inbreed'] = kinship(self.father, self.mother)
            return self.attrib['inbreed']

    # Functions for breeding
    #

    def gamete(self):
        """ Provides a set of half-genotypes to use with method fertilize """
        if not self.genotypes:
            self.get_genotypes()

        g = [recombine(chrom[0],
                       chrom[1],
                       self.chromosomes[i].genetic_map)
             for i, chrom in enumerate(self.genotypes)]

        return g

    def constrained_gamete(self, constraints, attempts=1000):
        # Constraints here is a list of ((location, index), alleles) tuples
        # for alleles that the gamete has to have
        for x in xrange(attempts):
            g = self.gamete()
            success = True
            for loc, allele in constraints:
                chr, pos = loc
                if g[chr][pos] != allele:
                    success = False
                    continue
            if success:
                return g
        else:
            raise IterationError('Ran out of constrained gamete attempts')

    @staticmethod
    def fertilize(father, mother):
        """
        Combines a set of half-genotypes (from method gamete) to a full
        set of genotypes
        """
        return [[x, y] for x, y in itertools.izip(father, mother)]

    # Phenotype functions
    def has_phenotype(self, phenotype):
        return (phenotype in self.phenotypes and
                self.phenotypes[phenotype] is not None)

    def predict_phenotype(self, trait):
        """ Predicts phenotypes from a given trait architecture and sets it """
        self.phenotypes[trait.name] = self.predicted_phenotype(trait)

    def predicted_phenotype(self, trait):
        return trait.predict_phenotype(self)

    def delete_phenotype(self, trait):
        "Removes a phenotype from the phenotype dictionary"
        if trait not in self.phenotypes:
            return
        del self.phenotypes[trait]
    
    def clear_phenotypes(self):
        """ Removes phenotypes """
        self.phenotypes = {}

    def genotype_as_phenotype(self, locus, minor_allele, label):
        """
        Creates a phenotype record representing the additive effect of a minor 
        allele at the specified locus. That is, the phenotype created is the 
        number of copies of the minor allele at the locus.

        Arguments:
        locus: the site to count minor alleles
        minor_allele: the allele to count
        label: The name of the phenotype to be added

        Returns: Nothing.
        """
        if not self.has_genotypes():
            self.phenotypes[label] = None
            return
            
        gt = self.get_genotype(locus)
        if is_missing_genotype(gt):
            val = None
        else: 
            val = gt.count(minor_allele)

        self.phenotypes[label] = val


    def _phenotypes_to_series(self):
        return pd.Series(self.phenotypes)
