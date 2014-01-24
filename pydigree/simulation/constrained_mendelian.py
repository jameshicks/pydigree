from pydigree.common import *
from simulation import *
from pydigree import paths
from pydigree import Individual

class ConstrainedMendelianSimulation(Simulation):
    def __init__(self, template=None, replications=1000):
        Simulation.__init__(self, template, replications)
        for inds in self.template:
            if ind.is_founder():
                continue
            if not ind.father.is_founder() and not ind.mother.is_founder():
                raise ValueError("ConstrainedMendelian only available"
                                 "for outbred pedigrees")
    
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
            geno_constraints = self.constraints['genotype']
            for ind in ped.founders():
                if ind not in geno_constraints:
                    ind.get_genotypes()
                else:
                    ind.get_constrained_genotypes(geno_constraints,
                                                  linkeq=True)
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
