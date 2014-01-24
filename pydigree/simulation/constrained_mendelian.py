from pydigree.common import *
from pydigree.path import paths
from pydigree import Individual

class ConstrainedMendelianSimulation(Simulation, pedigrees):
    def __init__(self):
        for inds in pedigrees.individuals():
            if ind.is_founder():
                continue
            if not ind.father.is_founder() and not ind.mother.is_founder():
                raise ValueError("ConstrainedMendelian only available"
                                 "for outbred pedigrees")
    
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

    def replicate(self):
        self.template.clear_genotypes()
        for ped in self.template:
            for x in ped.founders():
                x.label_genotypes()
            for iconstraint in self.constraints['ibd']:
                chromosome, position, allele = iconstraint['anchap']
                descendant, ancestor = iconstraint['inds']
                descent_path = random.choice(paths(ancestor, descendant))

                for pathindex, ind in enumerate(descent_path):
                    if ind.is_founder():
                        continue
                    fa, mo = ind.father, ind.mo
                    
                    if fa in descent_path:
                        paternal_gamete = fa.constrained_gamete(constraints)
                    else:
                        paternal_gamete = fa.gamete()
                    if mo in descent_path:
                        maternal_gamete = fa.constrained_gamete(constraints)
                    else:
                        maternal_gamete = mo.gamete()

                    genotypes = Individual.fertilize(paternal_gamete,
                                                     maternal_gamete)
                    ind._set_genotypes(genotypes)
            # Get genotypes for everybody else that we're not constraining.
            for ind in ped:
                ind.get_genotypes()

            # Now replace the label genotypes in founders with real ones.
            for founder in ped.founders():
                founder.clear_genotypes()
                founder.get_genotypes()
                # TODO: make sure founder genotype constraints get set
                pass

            # Now replace the label genotypes in the nonfounders with the
            # genotypes of the founders
            for nf in ped.nonfounders():
                for chromoidx, chromosome in enumerate(nf.genotypes()):
                    for chromaidx, chromatid in enumerate(chromosome):
                        for midx, marker in chromatid:
                            founder,which = marker
                            which = 0 if which == 'P' else 1
                            founder_marker = founder[chromoidx][which][midx]
                            nf.genotypes[chromoidx][chromaidx][midx] = founder_marker
