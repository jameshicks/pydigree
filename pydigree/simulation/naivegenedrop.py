from __future__ import division

from architecture import Architecture
from simulation import Simulation, SimulationError

class NaiveGeneDroppingSimulation(Simulation):
    def __init__(self, template=None, replications=1000):
        Simulation.__init__(self, template, replications)

    def replicate(self):
        self.template.clear_genotypes()

        for i in peds.individuals():
            if i.is_founder():
                i.get_genotypes(linkeq=self.linkeq)

        for ped in self.template:
            for attempt in xrange(self.genedrop_attempts):
                for ind in ped:
                    if not ind.is_founder():
                        ind.clear_genotypes()
                    ind.get_genotypes()
                if self.trait:
                    accuracy = self.predicted_trait_accuracy()
                    if accuracy >= accuracy_threshold:
                        continue
                break
            else:
                raise SimualtionError('Ran out of gene dropping attempts!')
            
