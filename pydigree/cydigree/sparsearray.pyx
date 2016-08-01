from collections import Sequence

cdef class SparseArray:
    # I imagine at some point I'll have to rewrite this whole thing
    # as a red-black tree or something, but until then we'll just bsearch
    # a sorted list

    cdef readonly list container
    cdef readonly object refcode
    cdef readonly Py_ssize_t size

    def __init__(self, size, refcode=None, initial=None):
        self.container = []
        self.refcode = refcode
        self.size = size

    def __len__(self):
        return self.size

    property values:
        def __get__(self):
            return [x.value for x in self.container]

        def __set__(self, vals):
            cdef Py_ssize_t n = len(vals)
            cdef Py_ssize_t i

            for i in range(n):
                self.container[i].value = vals[i]

    property indices:
        def __get__(self):
            return [x.index for x in self.container]

        def __set__(self, vals):
            cdef Py_ssize_t n = len(vals)
            cdef Py_ssize_t i

            for i in range(n):
                self.container[i].index = vals[i]
    
    cdef Py_ssize_t bsearch(self, Py_ssize_t idx_saught):
        'Return the leftmost internal index for a saught full-index'
        cdef Py_ssize_t high = len(self.container) - 1
        cdef Py_ssize_t low = 0
        cdef Py_ssize_t mid = (high + low) / 2
        cdef SparseArrayElement pivot
        
        while low <= high:
            mid = (high + low) / 2
            pivot = self.container[mid]

            if pivot.index > idx_saught:
                high = mid - 1
            elif pivot.index < idx_saught:
                low = mid + 1
            else:
                return mid

        return low

    cdef inline Py_ssize_t fix_index(self, Py_ssize_t index): 
        if index >= 0:
            return index 
        else:
            return self.size + index

    def __getitem__(self, index):
        if type(index) is slice:
            return self._get_slice(index.start, index.stop)
        elif isinstance(index, Sequence) and type(index[0]) is bool:
            return self._get_bool_indices(index)
        elif isinstance(index, Sequence) and type(index[0]) is int:
            return self._get_multiple_indices(index)
        else:
            return self._get_single_item(index)

    cdef _get_single_item(self, index):
        index = self.fix_index(index)
        if not 0 <= index < self.size:
            raise IndexError('index out of range')     

        if len(self.container) == 0:
            return self.refcode
        
        cdef Py_ssize_t putative_idx = self.bsearch(index)
        cdef SparseArrayElement element = self.container[putative_idx]
        
        if element.index == index:
            return element.value
        else:
            return self.refcode

    cdef _get_slice(self, Py_ssize_t start, Py_ssize_t stop):
        cdef Py_ssize_t slicelen = stop - start
        cdef list nonsparsevals = [x for x in self.container if start <= x.index < stop]
        cdef SparseArray sliced = SparseArray(slicelen, self.refcode)
        sliced.container = nonsparsevals

        return sliced

    cdef _get_bool_indices(self, bools):
        if len(bools) != self.size:
            raise IndexError('Bool index not right length')
        cdef set wantedsites = {i for i,x in enumerate(bools) if x}  
        cdef Py_ssize_t size = len(wantedsites)
        ns = SparseArray(size, self.refcode)
        ns.container = [x for x in self.container if x.index in wantedsites]
        return ns

    cdef _get_multiple_indices(self, indices):
        vals = [self[i] for i in indices]
        cdef SparseArray ns = SparseArray.from_sequence(vals, self.refcode)
        return ns

    def __setitem__(self, index, object value):
        if type(index) is slice:
            self._set_slice(index.start, index.stop, value)

        elif isinstance(index, Sequence):
            if type(index[0]) is bool:
                self._set_bool_indices(index, value)
            elif type(index[0]) is int:
                self._set_multiple_values(index, value)

        else:
            self._set_value(index, value)

    cdef inline _set_value(self, Py_ssize_t index, object value):
        cdef SparseArrayElement newelement
        index = self.fix_index(index)
        if not 0 <= index < self.size:
            raise IndexError('index out of range')  

        if len(self.container) == 0:
            newelement = SparseArrayElement(index, value)
            self.container.append(newelement)
            return 

        cdef Py_ssize_t internal_index = self.bsearch(index)

        if internal_index >= len(self.container):
            newelement = SparseArrayElement(index, value)
            self.container.append(newelement)
            return 

        cdef SparseArrayElement element = self.container[internal_index]
        cdef Py_ssize_t insertion_point
        # There are four scenarios here:
        # 1) There's a non-sparse element at `index` and we have to 
        #    change the value to another non-sparse value
        # 2) There's a non-sparse element at `index` and we need to change
        #    it to the sparse value (i.e. remove the SparseArrayElement)
        # 3) There's sparsity at the index and we have to put something there
        # 4) There's sparsity at the index and we have to leave it sparse

        if index == element.index:
            if value == self.refcode:
                # Scenario 2
                del self.container[internal_index]
            else:
                # Scenario 1
                element.value = value
        # Scenario 3 
        if index != element.index:
            newelement = SparseArrayElement(index, value)
            if element.index < index:
                insertion_point = internal_index + 1
            else:
                insertion_point = internal_index
            self.container.insert(insertion_point, newelement)

        # Scenario 4 (sparse to sparse) doesnt need anything

    cdef _set_slice(self, Py_ssize_t start, Py_ssize_t stop, object value):
        cdef list before = [x for x in self.container if x.index < start]
        cdef list after = [x for x in self.container if x.index >= stop]
        cdef list mid

        if isinstance(value, SparseArray):
            mid = [x for x in value.container if start <= (x.index+start) < stop]       
        elif isinstance(value, Sequence):
            raise NotImplementedError
        else:
            mid = [SparseArrayElement(i, value) for i in range(start, stop)]

        self.container = before + mid + after

    cdef _set_bool_indices(self, bools, vals):
        cdef Py_ssize_t i, fullidx
        cdef object x
        cdef list trueidx = [i for i,x in enumerate(bools) if x]

        if not len(bools) == self.size:
            raise IndexError('Mask not same size as array')        

        self._set_multiple_values(trueidx, vals)

    cdef _set_multiple_values(self, indices, value):
        cdef Py_ssize_t fullidx, i
        cdef Py_ssize_t nidx = len(indices)
        

        if isinstance(value, Sequence):
            if not nidx == len(value):
                raise IndexError('Index and values different sizes')
            
            for i in range(nidx):
                fullidx = indices[i]
                self[fullidx] = value[i]
        else:
            for i in range(nidx):
                fullidx = indices[i]
                self[fullidx] = value

    def __richcmp__(self, other, int op):
        # Op codes
        # <   0
        # ==  2
        # >   4
        # <=  1
        # !=  3
        # >=  5

        if op == 2: # Equality
            return self._eq(other)
        raise NotImplementedError

    def _eq(self, val):
        return self._eqval(val)

    cdef _eqval(self, val):
        sparseval = self.refcode == val
        cdef SparseArray eq = SparseArray(self.size, sparseval)
        if sparseval:
            eq.container = [SparseArrayElement(x.index, False) for 
                            x in self.container if x.value != val]
        else:
            eq.container = [SparseArrayElement(x.index, True) for 
                            x in self.container if x.value == val]
        return eq

    def tolist(self):
        cdef Py_ssize_t i
        cdef set denselocs = {x.index for x in self.container}
        cdef list dense = [(self[i] if i in denselocs else self.refcode) for i in range(self.size)]
        return dense

    @staticmethod
    def from_sequence(seq, refcode): 
        cdef Py_ssize_t n = len(seq)
        cdef SparseArray sa = SparseArray(n, refcode)
        sa.container = [SparseArrayElement(i,x) for i,x in 
                        enumerate(seq) if x != refcode]
        return sa

cdef class SparseArrayElement:
    cdef readonly Py_ssize_t index
    cdef readonly object value

    def __init__(self, index, value):
        self.index = index
        self.value = value

    def __repr__(self):
        return 'SparseArrayElement({},{})'.format(self.index, repr(self.value))

    def __copy__(self):
        return SparseArrayElement(self.index, self.value)

    def __richcmp__(self, other, int op):
        # Op codes
        # <   0
        # ==  2
        # >   4
        # <=  1
        # !=  3
        # >=  5

        if not isinstance(other, SparseArrayElement):
            t = type(other)
            raise ValueError("Can't compare SparseArrayElement and {}".format(t))

        if op == 2: # Equality
            return self.__eq(other)
        elif op == 3: # Not equality
            return not self._eq(other)

    cpdef bint __eq(self, other):
        res = self.index == other.index and self.value == other.value
        return res