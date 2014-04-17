from __future__ import division

from pydigree.simulation.architecture import Architecture
from pydigree.simulation.simulation import Simulation, SimulationError


class NaiveGeneDroppingSimulation(Simulation):
    def __init__(self, template=None, replications=1000):
        Simulation.__init__(self, template, replications)
        self.genedrop_attempts = 1000

    def replicate(self, verbose=None):
        for ind in self.template.individuals():
            ind.clear_genotypes()
        for ped in self.template:

            for attempt in xrange(self.genedrop_attempts):

                for ind in ped:
                    ind.clear_genotypes()

                for founder in ped.founders():
                    if founder in self.constraints['genotype']:
                        founder.get_constrained_genotypes(
                            self.constraints['genotype'][founder],
                            linkeq=True)
                    else:
                        founder.get_genotypes()

                for ind in ped:
                    if not ind.is_founder():
                        ind.clear_genotypes()

                for ind in ped:
                    ind.get_genotypes()

                if self.trait:
                    accuracy = self.predicted_trait_accuracy(ped)
                    if accuracy < self.accuracy_threshold:
                        continue
                if verbose:
                    print 'Success (%s%%) after %s attempts' % (accuracy * 100,
                                                                attempt)
                break
            else:
                raise SimulationError('Ran out of gene dropping attempts!')
