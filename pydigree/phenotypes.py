"A phenotype holder"

import pandas as pd

class Phenotypes(object):
    """
    A container for the set of phenotypes and exposures associated with an 
    Individual in a population
    """
    def __init__(self, data=None):
        self.data = dict(data) if data is not None else dict()

    def __contains__(self, key):
        return self.has_phenotype(key)

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, val):
        self.data[key] = val

    def __delitem__(self, key):
        self.clear_phenotype(key)

    def get(self, key, default):
        "Gets a phenotype or default value if not present"
        return self.data.get(key, default)

    def keys(self):
        "Iterate over the available phenotype names"
        return self.data.keys()

    def values(self):
        "Iterate over the available phenotype values"
        return self.data.values()

    def items(self):
        "Iterate over phenotype name, value pairs"
        return self.data.items()

    def has_phenotype(self, key):
        """
        Does the individual have a certain phenotype?

        :param phenotype: what phenotype are we looking for
        :type phenotype: string
    
        :returns: True if phenotype present and not None
        :rtype: bool 
        """
        return key in self.data and self.data[key] != None
    
    def delete_phenotype(self, key):
        "Clears the phenotype from the object"
        try:
            del self.data[key]
        except KeyError:
            pass

    def clear(self):
        "Clears all phenotypes from the object"
        self.data = dict()

    def update(self, other):
        "Updates the current Phenotype object with data from the other"
        if isinstance(other, Phenotypes):
            self.data.update(other.data)
        else:
            self.data.update(other)

    def to_series(self):
        "Returns phenotypes as a pandas Series"
        return pd.Series(self.data)

