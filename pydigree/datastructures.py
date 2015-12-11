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

    @property
    def indices(self):
        return [x[0] for x in self.container]

    def getindex(self, full_idx, after=False):
        ''' 
        Returns the closest index in the container to the index in the full set
        after: Return the index after the last index found, for slicing 
        '''
        return bisect.bisect(self.indices, full_idx) + (1 if after else 0)

    def _getloc(self, key):
        idx = bisect.bisect_left(self.indices, key)

        if self.container[idx][0] != key:
            raise KeyError
        else:
            return self.container[idx][1]

    def _getslice(self, key):
        ''' Returns a list of pairs within the slice '''
        if key.step is not None:
            raise NotImplementedError

        start = self.getindex(key.start)
        stop = self.getindex(key.stop)
        print start,stop
        return self.container[start:stop]