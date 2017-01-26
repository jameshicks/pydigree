"A class representing individuals"

from pydigree.recombination import recombine
from pydigree.paths import kinship
from pydigree.common import flatten
from pydigree.genotypes import LabelledAlleles
from pydigree.exceptions import IterationError
from pydigree.phenotypes import Phenotypes

# TODO: Move this somewhere more useful
missing_genotype = (0, 0)

def is_missing_genotype(g):
    return g == missing_genotype


class Individual(object):
    '''
    An object for working with the phenotypes and genotypes of an individual
    in a genetic study or simulation
    '''
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
        self.phenotypes = Phenotypes()
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
        ''' 
        Add a child for this individual

        :param child: child object
        :type child: Individual
        '''
        self.children.append(child)

    def register_with_parents(self):
        ''' Inform the parent Individuals that this is their child '''
        if self.is_founder():
            return
        self.father.register_child(self)
        self.mother.register_child(self)

    @property
    def full_label(self):
        ''' 
        :return: the population label and individual label 
        :rtype: 2-tuple
        '''
        return (self.population.label, self.label)

    # Functions about genotypes
    #
    @property
    def chromosomes(self):
        ''' 
        Returns a list of the individuals ChromosomeTemplate objects 
        '''
        if self.pedigree is not None:
            return self.pedigree.chromosomes
        else:
            return self.population.chromosomes

    def _init_genotypes(self, blankchroms=True, dtype=None, sparse=True):
        """ 
        Initializes genotypes so that all genotypes are missing if blankchroms 
        is true, otherwise, just sets what would be the chromosome to None
        """
        gt = []
        if blankchroms:
            gt = [(chrom.empty_chromosome(dtype=dtype, sparse=sparse),
                   chrom.empty_chromosome(dtype=dtype, sparse=sparse))
                   for chrom in self.chromosomes]
        else:
            gt = [[None, None] for chrom in self.chromosomes]

        self.genotypes = gt

    def has_genotypes(self):
        """ 
        Tests if individual has genotypes

        :returns: genotype available
        :rtype: bool 
        """
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
    
        :param linkeq: should the retrieved genotypes be generated
            under linkage equilibrium?
        :type linkeq: bool

        :returns: Nothing

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
        Returns alleles at a position.

        :returns: Genotype tuple
        :rtype: tuple
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
            locus, chromatid, allele, _ = constraint
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
        Takes another individual object, merges/updates phenotypes with the 
        other individual object and replace self's genotypes with other's

        .. warning::
            This does not merge genotypes, it replaces them! Additionally,
            phenotypes in the other data overwrite the existing ones 

        :param other: the data to update with
        :type other: Individual
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
        Gives the individual label genotypes. 

        When label genotypes are transmitted on to the next generation, you 
        can see where each allele in the next generation came from. This is
        sseful for gene dropping simulations and reducing the memory footprint
        of forward-time simulation. 
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

        :returns: A collection of all the ancestors of the individual
        :rtype: set of Individuals
        """
        if self.is_founder():
            return set()
        return set([self.father, self.mother]) | \
            self.father.ancestors() | \
            self.mother.ancestors()

    def descendants(self):
        """
        Recursively searches for descendants.

        :returns: A collection of all the descendants of the individual
        :rtype: set of Individuals
        """
        return set(self.children + list(flatten([x.descendants()
                                                 for x in self.children])))

    def siblings(self, include_halfsibs=False):
        """
        Returns this individuals sibliings.

        :param include_halfsibs: Include half-siblings
        :type include_halfsibs: bool
        
        :returns: A collection of the individual's siblings
        :rtype: set of Individuals
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

        :returns: Label of the matriline founder
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
        
        :returns: Label of patriline founder
        """
        if self.is_founder():
            return self.label
        else:
            return self.father.patriline()

    @property
    def depth(self):
        """
        Returns the depth of an individual in a genealogy, a rough measure of
        what generation in the pedigree the individual is. Defined as:
        depth = 0 if individual is a founder, else the maximum of the
        depth of each parent

        :returns: Indiviual depth
        :rtype: integer
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
        """ 
        Provides a set of half-genotypes to use with method fertilize
        
        :returns: a collection of AlleleContainers
        :rtype: list
        """
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
        for _ in range(attempts):
            g = self.gamete()
            success = True
            for loc, allele in constraints:
                chrom, pos = loc
                if g[chrom][pos] != allele:
                    success = False
                    continue
            if success:
                return g
        raise IterationError('Ran out of constrained gamete attempts')

    @staticmethod
    def fertilize(father, mother):
        """
        Combines a set of half-genotypes (from method gamete) to a full
        set of genotypes
        """
        return [[x, y] for x, y in zip(father, mother)]

    def predict_phenotype(self, trait):
        """ Predicts phenotypes from a given trait architecture and sets it """
        self.phenotypes[trait.name] = self.predicted_phenotype(trait)

    def delete_phenotype(self, trait):
        "Removes a phenotype from the phenotype dictionary"
        self.phenotypes.delete_phenotype(trait)
    

    def genotype_as_phenotype(self, locus, minor_allele, label):
        """
        Creates a phenotype record representing the additive effect of a minor 
        allele at the specified locus. That is, the phenotype created is the 
        number of copies of the minor allele at the locus.

        :param locus: the site to count minor alleles
        :param minor_allele: the allele to count
        :param label: The name of the phenotype to be added

        :returns: void
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

