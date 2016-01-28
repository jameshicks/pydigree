import bisect


class SortedPairContainer(object):

    def __init__(self, pairs):
        self.container = list(pairs)

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self._getslice(key)
        elif isinstance(key, int):
            return self._getloc(key)
        else:
            raise TypeError("Invalid type for key")

    def __contains__(self, key):
        f = bisect.bisect_left(self.indices, key)
        return self.indices[f] == key

    def __delitem__(self, key):
        if not key in self:
            raise KeyError
        internal_idx = self.getindex(key)

        del self.container[internal_idx]

    @property
    def indices(self):
        return [x[0] for x in self.container]

    @property
    def values(self):
        return [x[1] for x in self.container]

    @property
    def items(self):
        return self.container

    def getindex(self, full_idx, after=False):
        ''' 
        Returns the closest index in the container to the index in the full set
        after: Return the index after the last index found, for slicing 
        '''

        # Locate the leftmost value exactly equal to x
        # From the python standard lib:
        # https://docs.python.org/2/library/bisect.html

        i = bisect.bisect_left(self.indices, full_idx)
        if i != len(self.container) and self.indices[i] == full_idx:
            return i
        raise KeyError

    def _getloc(self, key):
        idx = self.getindex(key)
        return self.container[idx][1]

    def _getslice(self, key):
        ''' Returns a list of pairs within the slice '''
        if key.step is not None:
            raise NotImplementedError

        indices = self.indices
        start = bisect.bisect_left(self.indices, key.start)
        stop = bisect.bisect_left(self.indices, key.stop)
        return self.container[start:stop]
