""" 
A base class for gene dropping simulations to inherit from
"""

from itertools import combinations_with_replacement

from pydigree.ibs import ibs
from pydigree.io.smartopen import smartopen
from pydigree.io.base import write_pedigree
from pydigree.io.plink import write_plink
from pydigree.io.base import write_phenotypes



class GeneDroppingSimulation(object):
    """ 
    A base class for gene dropping simulations to inherit from
    """
    def __init__(self, template=None, label=None, replications=1000, only=None):
        self.template = template
        self.label = label if label is not None else 'unlabeled'
        self.replications = replications
        self.accuracy_threshold = 0.9
        self.constraints = {'genotype': {}, 'ibd': {}}
        self.trait = None
        self.founder_genotype_hooks = []
        self.only = only
        self.trait = None

    def set_trait(self, architecture):
        self.trait = architecture
        if not self.trait.chromosomes:
            self.trait.chromosomes = self.template.chromosomes

    def replicate(self, **kwargs):
        raise NotImplementedError("This is a base class don't call me")

    def get_founder_genotypes(self, linkeq=True):
        geno_constraints = self.constraints['genotype']
        
        for ind in self.template.founders():
            ind.clear_genotypes()
            
            if ind not in geno_constraints:
                ind.get_genotypes(linkeq=linkeq)
            else:
                ind.get_constrained_genotypes(geno_constraints[ind],
                                              linkeq=linkeq)
        self.run_founder_genotype_hooks()

    def run(self, verbose=False, writeibd=False, output_predicate=None, compression=None, output_chromosomes=None):
        #write_map(self.template, '{0}.map'.format(self.label))
        for x in range(self.replications):
            print('Replicate {}/{}'.format(x+1, self.replications))
            self.replicate(
                verbose=verbose, writeibd=writeibd, replicatenumber=x)

            self.write_data(
                x, predicate=output_predicate, compression=compression, output_chromosomes=output_chromosomes)

    def write_data(self, replicatenumber, predicate=None, compression=None,
                   output_chromosomes=None):
        filename = '{0}-{1}'.format(self.label, (replicatenumber + 1))
        write_pedigree(self.template, filename + '.pedigrees')
        write_plink(self.template, filename, predicate=predicate,
                    mapfile=True, compression=compression,
                    output_chromosomes=output_chromosomes)
        write_phenotypes(self.template, filename + '.csv',
                         predicate=predicate)

    def _writeibd(self, replicatenumber):
        # Warning: Don't call this function! If the individuals in the pedigree dont have
        # LABEL genotypes, you're just going to get IBS configurations at each locus, not
        # actual IBD calculations.
        #
        # If you have data you want to identify IBD segments in, check
        # pydigree.sgs
        with smartopen('{0}-{1}.ibd.gz'.format(self.label, replicatenumber + 1), 'w') as of:
            for ped in self.template.pedigrees:
                for ind1, ind2 in combinations_with_replacement(ped.individuals, 2):
                    identical = []
                    for chrom_idx in range(ind1.chromosomes.nchrom()):
                        if ind1 == ind2:
                            genos = zip(*ind1.genotypes[chrom_idx])
                            ibd = [2 * (x == y) for x, y in genos]
                        else:
                            genos1 = zip(*ind1.genotypes[chrom_idx])
                            genos2 = zip(*ind2.genotypes[chrom_idx])
                            ibd = [ibs(g1, g2)
                                   for g1, g2 in zip(genos1, genos2)]
                        identical.extend(ibd)
                    outline = [ped.label, ind1.label, ind2.label] + identical
                    outline = ' '.join([str(x) for x in outline])
                    of.write('{}\n'.format(outline))

    def predicted_trait_accuracy(self, ped):
        calls = [(self.trait.predict_phenotype(ind), 
                  ind.phenotypes['affected'])
                 for ind in ped
                 if ind.phenotypes['affected'] is not None]
        # Remember: in python the bools True and False are actually alternate
        # names for the integers 1 and 0, respectively, so you can do
        # arithmetic with them if you so please. Here we sum up all the
        # correct predictions and divide by the number of predictions made.
        return sum(x == y for x, y in calls) / len(calls)

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
                    _, ped, ind, chrom, index, allele, chromatid, method = l
                    locus = (chrom, index)
                    ind = self.template[ped][ind]
                    self.add_genotype_constraint(ind, locus, allele,
                                                 chromatid, method)
                elif l[0].lower() == 'ibd':
                    _, ped, ind, ancestor, chrom, index, anc_chromatid = l
                    locus = (chrom, index)
                    ind = self.template[ped][ind]
                    ancestor = self.template[ped][ancestor]
                    self.add_ibd_constraint(ind, ancestor,
                                            locus, anc_chromatid)
                else:
                    raise ValueError('Not a valid constraint (%s)' % l[0])

    def add_genotype_constraint(self, ind, location, allele,
                                chromatid, method='set'):
        if not ind.is_founder():
            raise ValueError('Genotype constraints only for founders')
        if chromatid not in 'PM':
            raise ValueError('Not a valid haplotype. Choose P or M')
        chromatid = 1 if chromatid == 'M' else 0
        location = tuple(int(x) for x in location)
        allele = int(allele)
        if ind not in self.constraints['genotype']:
            self.constraints['genotype'][ind] = []
        c = (location, chromatid, allele, method)
        self.constraints['genotype'][ind].append(c)

    def add_ibd_constraint(self, ind, ancestor, location, anchap):
        if anchap not in 'PM':
            raise ValueError('Not a valid haplotype. Choose P or M')
        anchap = 1 if anchap == 'M' else 0
        location = tuple(int(x) for x in location)
        if ind not in self.constraints['ibd']:
            self.constraints['ibd'][ind] = []
        c = (ancestor, location, anchap)
        self.constraints['ibd'][ind].append(c)

    def add_founder_genotype_hook(self, func):
        self.founder_genotype_hooks.append(func)

    def run_founder_genotype_hooks(self):
        for hook in self.founder_genotype_hooks:
            for founder in self.template.founders():
                hook(founder)
