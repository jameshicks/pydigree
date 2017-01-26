"Genedropping with IBD constraints"
from pydigree.common import random_choice
from pydigree.genotypes import AncestralAllele
from .simulation import GeneDroppingSimulation 
from pydigree.exceptions import SimulationError
from pydigree import paths
from pydigree import Individual
import collections


class ConstrainedMendelianSimulation(GeneDroppingSimulation):
    """
    Performs a gene-dropping simulation constrained to a specific 
    IBD pattern
    """

    def __init__(self, template=None, label=None, replications=1000, only=None):
        GeneDroppingSimulation.__init__(self, template=template, label=label,
                            replications=replications, only=only)
        for ind in self.template.individuals:
            if ind.is_founder():
                continue
            if not (ind.father.is_founder() or ind.mother.is_founder()):
                raise ValueError("ConstrainedMendelian only available"
                                 "for outbred pedigrees")

    def replicate(self, writeibd=False, verbose=False, replicatenumber=0):
        "Creates a replicate from the simulation"
        self.template.clear_genotypes()
  
        for x in self.template.founders():
            x.label_genotypes()
        
        for ind in sorted(self.constraints['ibd'],
                          key=lambda x: x.depth, reverse=True):
            if ind.has_genotypes():
                # If the individual we're looking at has genotypes
                # already, we've seen them earlier while getting
                # genotypes for someone deeper in the pedigree
                continue
            
            constraints = self.constraints['ibd'][ind]

            # TODO: Multiple constraints per individual
            # Right now we're only using the first ([0]) constraint
            constraints = [(x[1], AncestralAllele(x[0], x[2])) 
                           for x in constraints]
            
            location, allele = constraints[0]
            ancestor = allele.ancestor
            descent_path = random_choice(paths(ancestor, ind))

            for path_member in descent_path:
                if path_member.is_founder():
                    continue
                fa, mo = path_member.parents()

                if fa in descent_path:
                    paternal_gamete = fa.constrained_gamete(constraints)
                else:
                    paternal_gamete = fa.gamete()
                if mo in descent_path:
                    maternal_gamete = mo.constrained_gamete(constraints)
                else:
                    maternal_gamete = mo.gamete()

                genotypes = Individual.fertilize(paternal_gamete,
                                                 maternal_gamete)
                path_member._set_genotypes(genotypes)
        
        # Get genotypes for everybody else that we're not constraining.
        for ind in self.template.individuals:
            ind.get_genotypes()

        if writeibd:
            self._writeibd(replicatenumber)

        # Now replace the label genotypes in founders with real ones.
        self.get_founder_genotypes()

        # Now replace the label genotypes in the nonfounders with the
        # genotypes of the founders
        if isinstance(self.only, collections.Callable):
            siminds = [x for x in self.template.nonfounders() if self.only(x)]
        else:
            siminds = self.template.nonfounders()

        for nf in siminds:
            nf.delabel_genotypes()

        # Predict phenotypes
        if self.trait:
            for ind in siminds:
                ind.predict_phenotype(self.trait)

        if verbose:
            for ind in siminds:
                print(ind, ind.get_genotype(location))
