from collections import Sequence

from cython.operator cimport dereference as deref, preincrement as inc

from libc.stdint cimport uint32_t, int8_t
from libcpp.map cimport map as stlmap
from libcpp.vector cimport vector


ctypedef uint32_t sparsekey
ctypedef int8_t sparseval

cdef bint sparseval_cmp(sparseval a, sparseval b, int op):
    """
    values for op
    <   0
    ==  2
    >   4
    <=  1
    !=  3
    >=  5
    """

    if op == 0:
        return a < b
    elif op == 2:
        return a == b
    elif op == 4:
        return a > b
    elif op == 1:
        return a <= b
    elif op == 3:
        return a != b
    elif op == 5:
        return a >= b

cdef class SparseArray:
    cdef public stlmap[sparsekey, sparseval] data
    cdef public sparsekey size
    cdef public sparseval ref

    def __cinit__(self, sparsekey array_size, sparseval refcode):
        self.size = array_size
        self.ref = refcode

    def __len__(self):
        return self.size

    def keys(self):
        return [x.first for x in self.data]

    def values(self):
        return [x.second for x in self.data]

    def items(self):
        return [(x.first, x.second) for x in self.data]

    cpdef bint any(self):
        return self.data.size() != 0

    cpdef bint all(self):
        return self.data.size() == self.size

    @staticmethod
    def from_dense(seq, sparseval refcode):
        cdef SparseArray out = SparseArray(len(seq), refcode)

        cdef int k = 0
        cdef sparseval v
        for v in seq:
            if v != refcode:
                out.set_item(k, v)

            k += 1

        return out 

    @staticmethod
    def from_items(seq, sparsekey size, sparseval refcode):
        cdef SparseArray out = SparseArray(size, refcode)

        cdef sparsekey k
        cdef sparseval v
        for k, v in seq:
            out.set_item(k, v)

        return out

    # Value getting
    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.get_slice(index)
        elif isinstance(index, Sequence): 
            if type(index[0]) is bool:
                return self.get_boolidx(index)
            else:
                return self.get_fancy(index)
        else:
            return self.get_item(index)

    cpdef sparseval get_item(self, int k):
        cdef stlmap[sparsekey, sparseval].iterator val = self.data.find(k)

        if val == self.data.end():
            return self.ref

        return deref(val).second

    cdef SparseArray get_fancy(self, indices):
        cdef sparsekey n = len(indices)
        cdef SparseArray out = SparseArray(n, self.ref)

        cdef int i = 0
        for i in range(n):
            out.set_item(i, self.get_item(indices[i]))

        return out

    cdef SparseArray get_boolidx(self, indices):
        cdef int n = len(indices)
        if n != self.size:
            raise ValueError

        cdef int sz = 0
        cdef int i = 0
        for i in range(n):
            if indices[i]: inc(sz)

        cdef SparseArray out = SparseArray(sz, self.ref)
        i = 0
        keepidx = 0
        for i in range(n):
            if indices[i]:
                out.set_item(keepidx, self.get_item(i))
                inc(keepidx)

        return out

    cpdef SparseArray get_slice(self, slice):
        cdef int start = slice.start 
        cdef int stop = slice.stop
        cdef SparseArray out = SparseArray(stop - start, self.ref) 

        for k, v in self.data:
            if start <= k < stop:
                out.set_item(k-start, v)
            if k >= stop:
                break
        return out 


    # Item Setting
    # 

    def __setitem__(self, index, object value):
        if isinstance(index, slice):
            self.set_slice(index, value)
        elif isinstance(index, Sequence):
            if type(index[0]) is bool:
                self.set_boolidx(index, value)
            else:
                self.set_fancy(index, value)
        else:
            self.set_item(index, value)

    cpdef void set_item(self, int k, int v):
        if v == self.ref:
            self.clear(k)
        else:
            self.data[k] = v

    cdef void set_slice_to_sparray(self, sparsekey start, sparsekey stop, SparseArray arr):
        if arr.ref != self.ref:
            self.set_slice(slice(start, stop), arr.ref)
        
        cdef sparsekey k
        cdef sparseval v
        for k,v in arr.data:
            self.set_item(start+k, v)


    cdef void set_slice(self, slice, object values):
        cdef sparsekey start = slice.start
        cdef sparsekey stop = slice.stop

        if isinstance(values, SparseArray):
            self.set_slice_to_sparray(start, stop, values)
            return

        self.clear_range(start, stop)

        cdef int nvals = 0
        cdef int i = 0
        cdef sparseval value

        if isinstance(values, Sequence):
            nvals = len(values)
            i = 0
            
            for i in range(nvals):
                value = values[i]
                if value != self.ref:
                    self.set_item(start + i, value)

        else:
            value = values 
            for i in range(start, stop):
                self.set_item(i, value)

    cdef void set_boolidx(self, indices, values):

        cdef int n = len(indices)
        cdef int i = 0
        cdef bint included = False
        
        cdef sparseval value = 0
        cdef int included_count = 0
        if isinstance(values, Sequence):
            
            for i in range(n):
                included = indices[i]
                
                if included:
                    value = values[included_count] 
                    self.set_item(i, value)
                    inc(included_count)

        else:
            value = values
            for i in range(n):
                included = indices[i]
                if included:
                    self.set_item(i, value)


    cdef void set_fancy(self, indices, values):
        cdef int i = 0
        cdef int n = len(indices)
        
        cdef sparseval value = 0
        if isinstance(values, Sequence):
            for i in range(n):
                value = values[i]
                self.set_item(indices[i], values[i])

        else: 
            value = values
            for i in range(n):
                self.set_item(indices[i], value)

    cpdef void clear(self, sparsekey k):
        self.data.erase(k)

    cpdef void clear_range(self, sparsekey start, sparsekey stop):
        self.data.erase(self.data.lower_bound(start), self.data.upper_bound(stop))

    cpdef sparsekey ndense(self):
        return self.data.size()

    cpdef sparsity(self):
        return 1 - self.ndense() / <float>self.size

    cpdef density(self):
        return self.ndense() / <float>self.size

    cpdef copy(self): 
        cdef SparseArray out = SparseArray(self.size, self.ref)
        cdef sparsekey k 
        cdef sparseval v 

        for k,v in self.data:
            out.data[k] = v

        return out

    def __richcmp__(self, other, int op):
       
        if isinstance(other, SparseArray):  
            if op == 2:
                return self.sparse_eq(other)
            elif op == 3:
                return self.sparse_eq(other).logical_not()
            else: 
                return self.sparse_cmp(other, op)

        elif  isinstance(other, Sequence):
            return self.dense_cmp(other, op)

        else:
            if op == 2:
                return self.eq_single(other)
            elif op == 3:
                return self.eq_single(other).logical_not()
            else:
                return self.cmp_single(other, op)
    
    cpdef SparseArray dense_cmp(self, other, op):
        cdef vector[sparseval] thisone = self.tolist()
        cdef SparseArray out = SparseArray(self.size, 0)

        cdef stlmap[sparsekey, sparseval].iterator it = self.data.begin()
        cdef int i = 0
        while it != self.data.end():
            while i < deref(it).first:
                if sparseval_cmp(other[i], self.ref, op):
                    out.set_item(i, True)
                inc(i)

            out.set_item(i, sparseval_cmp(deref(it).second, other[i], op))

            inc(it)
            inc(i)

        while i < len(other):
            if sparseval_cmp(other[i], self.ref, op):
                out.set_item(i, True)

            inc(i)

        return out

            
    cpdef SparseArray sparse_eq(self, SparseArray other):
        cdef SparseArray out = SparseArray(self.size, self.ref == other.ref)
        
        cdef sparsekey k
        cdef sparseval v

        # This can definitely be done faster
        for k,v in self.data:
            out.set_item(k, other.get_item(k) == v)

        for k,v in other.data:
            out.set_item(k, self.get_item(k) == v)

        return out

    cpdef SparseArray sparse_cmp(self, SparseArray other, int op):
        cdef SparseArray out = SparseArray(self.size, sparseval_cmp(self.ref, other.ref, op))

        cdef sparsekey k
        cdef sparseval v

        # This can definitely be done faster
        for k,v in self.data:
            out.set_item(k, sparseval_cmp(other.get_item(k), v, op))

        for k,v in other.data:
            out.set_item(k, sparseval_cmp(self.get_item(k), v, op))

        return out

    cpdef SparseArray eq_single(self, sparseval val):
        cdef SparseArray out = SparseArray(self.size, self.ref == val)

        cdef sparsekey k
        cdef sparseval v

        cdef stlmap[sparsekey, sparseval].iterator it = self.data.begin()
        while it != self.data.end():
            k = deref(it).first
            v = deref(it).second
            out.set_item(k, val == v)
            inc(it)

        return out

    cpdef SparseArray cmp_single(self, sparseval val, int op):
        cdef SparseArray out = SparseArray(self.size, sparseval_cmp(self.ref, val, op))

        cdef sparsekey k
        cdef sparseval v

        cdef stlmap[sparsekey, sparseval].iterator it = self.data.begin()
        while it != self.data.end():
            k = deref(it).first
            v = deref(it).second
            out.set_item(k, sparseval_cmp(v, val, op))
            inc(it)

        return out

    cpdef SparseArray logical_not(self):
        cdef SparseArray out = SparseArray(self.size, not self.ref)
        
        for k, v in self.data:
            out.set_item(k, not v) 
        
        return out


    cpdef vector[sparseval] tolist(self):
        cdef vector[sparseval] out = vector[sparseval](self.size, self.ref)
        cdef sparsekey k
        cdef sparseval v

        for k,v in self.data:
            out[k] = v

        return out
