import random

from common import *
from pydigree import write_ped
from pydigree.individual import Individual

# A base class for simulations to inherit from
class Simulation(Object):
    
    def __init__(self, template=pedigrees, replications=1000):
        self.template = template
        self.replications = replications
        self.constraints = {'genotype': [], 'ibd': []}

    def set_trait(self, architecture):
        self.trait = architecture

    def replicate(self):
        raise NotImplementedError("This is a base class don't call me")

    def run(self):
        for x in xrange(self.replicates):
            self.replicate()
            self.write_data(x)

    def write_data(self, replicatenumber):
        filename = '{0}-{1}.ped'.format(self.prefix, replicatenumber)
        write_ped(self.template, filename)

    def predicted_trait_accuracy(self):
        calls = [(ind.predicted_phenotype(trait), ind.phenotypes['affected'])
                 for ind in ped if ind.phenotypes['affected'] is not None]
        # Remember: in python the bools True and False are actually alternate
        # names for the integers 1 and 0, respectively, so you can do
        # arithmetic with them if you so please. Here we sum up all the
        # correct predictions and divide by the number of predictions made.
        return sum(x==y for x,y in calls) / len(calls)

    def read_constraints(self, filename):

        if not self.template:
            raise ValueError()

        with open(filename) as f:
            for line in f:
                line = line.strip()

                if not line or line.startswith('#'):
                    continue

                l = line.split()
                if l[0].lower() == 'genotype':
                    type, ped, id, chr, index, allele, chromatid, method = l
                    locus = (chr, index)
                    ind = self.template[ped][id]
                    self.add_genotype_constraint(ind, locus, allele,
                                                 chromatid, method)
                elif l[0].lower() == 'ibd':
                    type, ped, id, ancestor, chr, index, anc_chromatid = l
                    locus = (chr, index)
                    ind = self.template[ped][id]
                    ancestor = self.template[ped][ancestor]
                    self.add_ibd_constraint(ind, ancestor, locus, anc_chromatid)
                else:
                    raise ValueError('Not a valid constraint (%s)' % l[0])

    def add_genotype_constraint(self, ind, location, allele,
                                chromatid, method='set'):
        if chromatid not in 'PM':
            raise ValueError('Not a valid haplotype. Choose P or M')
        self.constraints['genotype'].append('allele': allele, 'chromatid': chromatid,
                                            'ind': ind, 'location': location)

    def add_ibd_constraint(self, ind, ancestor, location, allele):
        self.constraints['ibd'].append({'anchap': (location[0], location[1], anchap),
                                        'inds': (ind, ancestor)})


class SimulationError(Exception):
    pass
