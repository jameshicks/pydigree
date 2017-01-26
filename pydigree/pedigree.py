"A collection of individuals with fixed relationships"

import numpy as np

from pydigree.paths import fraternity
from pydigree.common import table
from pydigree.population import Population


class Pedigree(Population):
    "A collection of individuals with fixed relationships"
    
    def __init__(self, label=None):
        """
        Create a pedigree.

        :param label: pedigree label
        """
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

        :returns: bit size
        :rtype: pedigree
        """
        t = table([x.is_founder() for x in self.individuals])
        return 2 * t[False] - t[True]
   
    # Relationships
    #
    def kinship(self, id1, id2):
        """
        Get the Malecot coefficient of coancestry for two individuals in
        the pedigree. These are calculated recursively.
        For pedigree objects, results are stored to reduce the calculation
        time for kinship matrices.

        
        :param id1: the label of a individual to be evaluated
        :param id2: the label of a individual to be evaluated

        :returns: Malecot's coefficient of coancestry
        :rtype: float

        Reference:
        Lange. Mathematical and Statistical Methods for Genetic Analysis.
        1997. Springer.
        """
        pair = frozenset([id1, id2])
        if pair in self.kinmat:
            return self.kinmat[pair]
        if id1 is None or id2 is None:
            return 0

        # Since with pedigree objects we're typically working with IDs,
        # I define these functions to get parents for IDs by looking them
        # up in the pedigree
        def fa(lab):
            return (self[lab].father.label
                    if self[lab].father is not None else None)

        def mo(lab):
            return (self[lab].mother.label
                    if self[lab].mother is not None else None)

        # Use tuples here to take advantage of the implicit tuple ordering
        # With depth as the first item, it assures that descendants aren't
        # listed before their ancestors.
        t1 = self[id1].depth, id1
        t2 = self[id2].depth, id2
        if id1 == id2:
            k = (1 + self.kinship(fa(id1), mo(id1))) / 2.0
        elif t1 < t2:
            k = (self.kinship(id1, fa(id2)) + self.kinship(id1, mo(id2))) / 2.0
        else:
            k = (self.kinship(id2, fa(id1)) + self.kinship(id2, mo(id1))) / 2.0
        self.kinmat[pair] = k
        return k

    def fraternity(self, id1, id2):
        """
        Like Pedigree.kinship, this is a convenience function for getting
        fraternity coefficients for two pedigree memebers by their ID label.

        This is a wrapper for paths.fraternity

        :param id1: the label of a individual to be evaluated
        :param id2: the label of a individual to be evaluated

        :returns: coefficient of fraternity
        :rtype: float
        """
        pair = frozenset([id1, id2])
        if pair not in self.fratmat:
            f = f = fraternity(self[id1], self[id2])
            self.fratmat[pair] = f
            return f
        else:
            return self.fratmat[pair]

    def inbreeding(self, indlab):
        """
        Like Pedigree.kinship, this is a convenience function for getting
        inbreeding coefficients for individuals in pedigrees by their id
        label. As inbreeding coefficients are the kinship coefficient of
        the parents, this function calls Pedigree.kinship to check for
        stored values.

        :param id: the label of the individual to be evaluated
    
        :returns: inbreeding coefficient
        :rtype: a double
        """
        ind = self[indlab]
        if ind.is_founder():
            return 0.0
        if ind.father.is_founder() or ind.mother.is_founder():
            return 0.0
        return self.kinship(ind.father.label, ind.mother.label)

    def additive_relationship_matrix(self, ids=None):
        """
        Calculates an additive relationship matrix (the A matrix) for
        quantiatitive genetics.

        A_ij = 2 * kinship(i,j) if i != j.
        (See the notes on function 'kinship')
        A_ij = 1 + inbreeding(i) if i == j
        (inbreeding(i) is equivalent to kinship(i.father,i.mother))


        :param ids: IDs of pedigree members to include in the matrix

        Important: if not given, the rows/columns are all individuals in the
        pedigree, sorted by id. If you're not sure about this, try
        sorted(x.label for x in ped) to see the ordering.

        :returns: additive relationship matrix
        :rtype: matrix
        """
        if not ids:
            ids = sorted(x.label for x in self.individuals)
        else:
            ids = [label for ped, label in ids if ped == self.label and
                   label in self.population.keys()]
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


        :param ids: IDs of pedigree members to include in the matrix

        Important: if not given, the rows/columns are all individuals in the
        pedigree, sorted by id. If you're not sure about this, try
        sorted(x.label for x in ped) to see the ordering.

        :returns: dominance relationship matrix
        :rtype: matrix
        """
        if not ids:
            ids = sorted(x.label for x in self.individuals)
        else:
            ids = [label for ped, label in ids if ped == self.label and
                   label in self.population.keys()]
        mat = []
        for a in ids:
            row = []
            for b in ids:
                if a == b:
                    row.append(1)
                else:
                    row.append(self.fraternity(a, b))
            mat.append(row)
        return np.matrix(mat)

    def mitochondrial_relationship_matrix(self, ids=None):
        """
        Calculates the mitochondrial relationship matrix.
        M_ij = 1 if matriline(i) == matriline(j)

        :param ids: IDs of pedigree members to include in the matrix

        Important: if not given, the rows/columns are all individuals in the
        pedigree, sorted by id. If you're not sure about this, try
        sorted(x.label for x in ped) to see the ordering.

        Returns: A numpy matrix

        Reference:
        Liu et al. "Association Testing of the Mitochondrial Genome Using
        Pedigree Data". Genetic Epidemiology. (2013). 37,3:239-247
        """
        if not ids:
            inds = sorted((x for x in self.individuals), key=lambda x: x.label)
        else:
            inds = [self[id] for id in ids]
        mat = []
        for a in inds:
            row = [1 if a.matriline() == b.matriline() else 0
                   for b in ids]
            mat.append(row)
        return np.matrix(mat)

    # Gene dropping
    #
    def simulate_ibd_states(self, inds=None):
        """
        Simulate IBD patterns by gene dropping: Everyone's genotypes reflect
        the founder chromosome that they received the genotype from. You can
        then use the ibs function to determine IBD state. This effectively an
        infinite-alleles simulation.

        Returns: Nothing
        """
        self.clear_genotypes()
        for x in self.founders():
            x.label_genotypes()
        if inds:
            for x in inds:
                x.get_genotypes()
        else:
            for x in self.nonfounders():
                x.get_genotypes()
