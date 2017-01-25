from scipy.sparse import block_diag

from pydigree.pedigree import Pedigree
from pydigree.individualcontainer import IndividualContainer


class PedigreeCollection(IndividualContainer):

    def __init__(self, peds=None):
        self.container = {}
        if peds:
            for ped in peds:
                self.add_pedigree(ped)

    def __getitem__(self, key):
        if isinstance(key, tuple) or isinstance(key, list):
            return self.container[key[0]][key[1]]
        return self.container[key]

    def __contains__(self, item):
        return item in self.pedigrees

    def __len__(self):
        return len(self.container)

    def __setitem__(self, key, value):
        self.container[key] = value

    def __delitem__(self, key):
        del self.container[key]

    def keys(self):
        return list(self.container.keys())

    def add_pedigree(self, ped):
        if not isinstance(ped, Pedigree):
            raise ValueError('{} not of type Pedigree')
        elif ped.label in list(self.container.keys()):
            raise ValueError(
                'Pedigree label {} already in collection'.format(ped.label))
        else:
            self[ped.label] = ped

    @property
    def individuals(self):
        '''
        Returns a list of the individuals represented by all pedigrees, 
        sorted by pedigree label, id label 
        '''
        inds = []
        for pedigree in sorted(self.pedigrees, key=lambda x: x.label):
            inds.extend(sorted((x for x in pedigree.individuals), key=lambda x: x.label))
        return inds

    @property
    def pedigrees(self):
        '''
        Returns a list of the pedigree objects contained in the collection
        '''
        return list(self.container.values())


    def _getindividual(self, label):
        for x in self.individuals:
            if x.label == label:
                return x
        raise KeyError('Individual not in collection')

    @property
    def chromosomes(self):
        k = list(self.container.keys())[0]
        return self.container[k].chromosomes

    def add_chromosome(self, chrom):
        for x in self.pedigrees:
            x.add_chromosome(chrom)

    def update(self, pop):
        for ped in self.pedigrees:
            ped.chromosomes = pop.chromosomes
            try:
                ped.update(pop[ped.label])
            except KeyError:
                continue

    # Matrix functions
    ###
    def additive_relationship_matrix(self, ids=None):
        """
        Returns a block diagonal matrix of additive relationships
        for each pedigree.

        See notes on Pedigree.additive_relationship_matrix
        """
        mats = [x.additive_relationship_matrix(ids) for x in
                sorted(self.pedigrees, key=lambda x: x.label)]
        mats = [x for x in mats if x.size > 0]
        return block_diag(mats, format='bsr')

    def dominance_relationship_matrix(self, ids=None):
        """
        Returns a block diagonal matrix of dominance relationships
        for each pedigree.

        See notes on Pedigree.dominance_relationship_matrix
        """
        mats = [x.dominance_relationship_matrix(ids) for x in
                sorted(self.pedigrees, key=lambda x: x.label)]
        mats = [x for x in mats if x.size > 0]
        return block_diag(mats, format='bsr')

    def mitochondrial_relationship_matrix(self, ids=None):
        """
        Returns a block diagonal matrix of mitochondrial relationships
        for each pedigree.

        See notes on Pedigree.mitochondrial_relationship_matrix
        """
        mats = [x.mitochondrial_relationship_matrix(ids) for x in
                sorted(self.pedigrees, key=lambda x: x.label)]
        return block_diag(mats, format='bsr')
