#!/usr/bin/env/python

import numpy as np

from common import *
from population import Population
from paths import kinship, fraternity


class Pedigree(Population):
    def __init__(self, label=None):
        Population.__init__(self)
        self.label = label
        self.kinmat = {}
        self.fratmat = {}

    def __prepare_nonfounder_contraint(self, con):
        if not con:
            return lambda x: x.is_founder()
        else:
            return lambda x: x.is_founder() and con(x)

    def bit_size(self):
        """
        Returns the bit size of the pedigree. The bitsize is defined as 2*n-f
        where n is the number of nonfounders and f is the number of founders.
        This represents the number of bits it takes to represent the
        inheritance vector in the Lander-Green algorithm.

        Returns: an integer
        """
        t = table([x.is_founder() for x in self])
        return 2 * t[False] - t[True]

    ### Frequency
    def alleles(self, location, constraint=None, nonfounders=False):
        """
        Like Population.alleles, except constrained to founder individuals

        If nonfounders is True, it's just a call to Population.alleles.

        Returns: a list
        """
        if nonfounders:
            return Population.alleles(self, location, constraint)
        con = self.__prepare_nonfounder_constraint(constraint)
        return Population.alleles(self, location, con)

    def allele_frequency(self, location, allele,
                         constraint=None, nonfounders=False):
        """
        Like Population.alleles, except constrained to founder individuals.
        If nonfounders is True, it's just a call to Population.alleles.

        Returns: a double
        """
        if nonfounders:
            return Population.allele_frequency(self,
                                               location, allele, constraint)
        constraint = self.__prepare_nonfounder_constraint(constraint)
        return Population.allele_frequency(self,
                                           location, allele,
                                           constraint=constraint)

    def ld(self):
        raise NotImplementedError('LD not meaningful for pedigrees?')

    ### Relationships
    ###
    def kinship(self, id1, id2):
        """
        Get the Malecot coefficient of coancestry for two individuals in
        the pedigree. (See notes in pydigree.paths.kinship). For pedigree
        objects, results are stored to reduce the calculation time for kinship
        matrices.

        This is a convenience wrapper for paths.kinship, which takes individual
        objects as arguments. This function takes id labels and looks them up
        in the pedigree, and calls paths.kinship on those individual objects.

        Returns: a double
        """
        pair = frozenset([id1, id2])
        if pair not in self.kinmat:
            k = kinship(self[id1], self[id2])
            self.kinmat[pair] = k
            return k
        else:
            return self.kinmat[pair]

    def fraternity(self, id1, id2):
        """
        Like Pedigree.kinship, this is a convenience function for getting
        fraternity coefficients for two pedigree memebers by their ID label.

        This is a wrapper for paths.fraternity

        Returns: a double
        """
        pair = frozenset([id1, id2])
        if pair not in self.fratmat:
            f = fraternity(self[id1], self[id2])
            self.fratmat[pair] = f
            return f
        else:
            return self.fratmat[pair]

    def inbreeding(self, id):
        """
        Like Pedigree.kinship, this is a convenience function for getting
        inbreeding coefficients for individuals in pedigrees by their id
        label. As inbreeding coefficients are the kinship coefficient of
        the parents, this function calls Pedigree.kinship to check for
        stored values.

        Returns: a double
        """
        ind = self[id]
        if ind.is_founder():
            return 0.0
        if ind.father.is_founder() or ind.mother.is_founder():
            return 0.0
        return self.kinship(ind.father.id, ind.mother.id)

    def additive_relationship_matrix(self, ids=None):
        """
        Calculates an additive relationship matrix (the A matrix) for
        quantiatitive genetics.

        A_ij = 2 * kinship(i,j) if i != j.
        (See the notes on function 'kinship')
        A_ij = 1 + inbreeding(i) if i == j
        (inbreeding(i) is equivalent to kinship(i.father,i.mother))

        Arguments
        -----
        ids: IDs of pedigree members to include in the matrix

        Important: if not given, the rows/columns are all individuals in the
        pedigree, sorted by id. If you're not sure about this, try
        sorted(x.id for x in ped) to see the ordering.

        Returns: a numpy matrix
        """
        if not ids:
            ids = sorted(x.id for x in self)
        mat = []
        for a in ids:
            row = []
            for b in ids:
                if a == b:
                    row.append(1 + self.inbreeding(a))
                else:
                    row.append(2 * self.kinship(a, b))
            mat.append(row)
        return np.matrix(mat)

    def dominance_relationship_matrix(self, ids=None):
        """
        Calculates the dominance genetic relationship matrix (the D matrix)
        for quantitative genetics.

        D_ij = fraternity(i,j) if i != j
        D_ij = 1 if i == j

        Arguments
        -----
        ids: IDs of pedigree members to include in the matrix

        Important: if not given, the rows/columns are all individuals in the
        pedigree, sorted by id. If you're not sure about this, try
        sorted(x.id for x in ped) to see the ordering.

        Returns: A numpy matrix
        """
        if not ids:
            ids = sorted(x.id for x in self)
        mat = []
        for a in ids:
            row = []
            for b in ids:
                if a == b:
                    row.append(1)
                else:
                    row.append(.25 * self.fraternity(a, b))
            mat.append(row)
        return np.matrix(mat)

    def mitochondrial_relationship_matrix(self, ids=None):
        """
        Calculates the mitochondrial relationship matrix.
        M_ij = 1 if matriline(i) == matriline(j)

        WARNING: These matrices can often be singular!
        Consider the pedigree:
               F---M
                 |
          ---------------
          |      |      |
          C1     C2     C3

        There are only two matrilines in here. One belonging to whoever the
        mother of the father is (not specified in this pedigree) and one
        transmitted from the mother to all three children. The mitochondrial
        relationship matrix would then look like this:


            / 1 0 0 0 0 \
            | 0 1 1 1 1 |  This matrix is obviously singlular and will not
        M = | 0 1 1 1 1 |  be useful in variance component estimation.
            | 0 1 1 1 1 |
            \ 0 1 1 1 1 /

        Arguments
        -----
        ids: IDs of pedigree members to include in the matrix

        Important: if not given, the rows/columns are all individuals in the
        pedigree, sorted by id. If you're not sure about this, try
        sorted(x.id for x in ped) to see the ordering.

        Returns: A numpy matrix

        Reference:
        Liu et al. "Association Testing of the Mitochondrial Genome Using
        Pedigree Data". Genetic Epidemiology. (2013). 37,3:239-247
        """
        if not ids:
            inds = sorted((x for x in self), key=lambda x: x.id)
        else:
            inds = [self[id] for id in ids]
        mat = []
        for a in inds:
            row = [1 if a.matriline() == b.matriline() else 0
                   for b in ids]
            mat.append(row)
        return np.matrix(mat)

    ### Gene dropping
    ###
    def simulate_ibd_states(self):
        """
        Simulate IBD patterns by gene dropping: Everyone's genotypes reflect
        the founder chromosome that they received the genotype from. You can
        then use misc.ibs to determine IBD state. This effectively an
        infinite-alleles simulation.

        Returns: Nothing
        """
        for x in self:
            if x.is_founder():
                x.label_genotypes()
        for x in self:
            if x.is_founder():
                continue
            x.get_genotypes()
