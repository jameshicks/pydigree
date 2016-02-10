
class SortedList(object):
    def __init__(self, x, key=None):
        self.container = sorted(x[:]) 
        if key is None: 
            self.key = lambda x: x
        else:
            self.key = key

    def __getitem__(self, key): 
        return self.container[key]

    def __setitem__(self, key, value):
        raise NotImplementedError('Use add and remove methods')

    def __delitem__(self, key):
        del self[key] 

    def __len__(self):
        return len(self.container)

    def __contains__(self, value):
        return value in self.container

    def extend(self, extension):
        if not isinstance(extension, SortedList):
            raise TypeError    
        
        T = self.key
        if T(self.container[-1]) > T(extension[0]):
            raise ValueError('Extension breaks sorting order')
        else:
            self.container.extend(extension)

    def append(self, value):
        T = self.key

        if T(self.container[-1]) > T(value):
            raise ValueError('Value breaks sorting order')
        else:
            self.container.append(value)

    def clear(self):
        self.container = []



    