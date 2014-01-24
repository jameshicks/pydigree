import random

from common import *
from pydigree.individual import Individual

# A base class for simulations to inherit from
class Simulation(Object):
    
    def __init__(self):
        pass

    def replicate(self):
        raise NotImplementedError("This is a base class don't call me")

    def run(self):
        for x in xrange(self.replicates):
            self.replicate()
            self.write_data(x)

class SimulationError(Exception):
    pass
