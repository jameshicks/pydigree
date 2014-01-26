from __future__ import division

from architecture import Architecture
from simulation import Simulation, SimulationError

class NaiveGeneDroppingSimulation(Simulation):
    def __init__(self, template=None, replications=1000):
        Simulation.__init__(self, template, replications)
        self.genedrop_attempts = 1000

    def replicate(self):
        self.template.clear_genotypes()

        for ped in self.template:
            for founder in ped.founders():
                if founder in self.constraints['genotype']:
                    founder.get_constrained_genotypes(self.constraints['genotype'][founder],
                                                      linkeq=True)
                else:
                    founder.get_genotypes()

            for attempt in xrange(self.genedrop_attempts):
                for ind in ped:
                    if not ind.is_founder():
                        ind.clear_genotypes()
                    ind.get_genotypes()
                if self.trait:
                    accuracy = self.predicted_trait_accuracy(ped)
                    if accuracy < self.accuracy_threshold:
                        continue
                break
            else:
                raise SimulationError('Ran out of gene dropping attempts!')
            
