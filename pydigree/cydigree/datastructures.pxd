from libc.stdint cimport uint32_t, int8_t, uint8_t

ctypedef uint32_t sparsekey

cdef inline bint compare(object a, object b, int op)

cdef class SparseArray:
    cdef readonly IntTree container
    cdef readonly int8_t refcode
    cdef readonly sparsekey size
    cdef inline sparsekey fix_index(self, sparsekey index)
    cdef object _get_single_item(self, sparsekey index)
    cdef SparseArray _get_slice(self, index)
    cdef _get_fancyidx(self, index)
    cdef _get_boolidx(self, index)
    cpdef void set_item(self, sparsekey index, int8_t value)
    cdef void _set_slice(self, sparsekey start, sparsekey stop, values)
    cdef void _set_slice_to_sparray(self, sparsekey start, sparsekey stop, SparseArray values)
    cdef void _set_boolidx(self, indices, values)
    cdef void _set_fancyidx(self, indices, value)
    cpdef _cmp_single(self, object value, int op)
    cpdef SparseArray _cmp_sequence(self, other, int op)
    cpdef SparseArray _cmp_sparray(self, SparseArray other, int op)
    cpdef bint any(self)
    cpdef bint all(self)
    cpdef SparseArray logical_not(self)
    cpdef double sparsity(self)

cdef struct IntTreeNode:
    IntTreeNode* left 
    IntTreeNode* right
    IntTreeNode* parent
    sparsekey key
    int8_t value 
    int8_t height

cdef class IntTree(object):
    cdef IntTreeNode* root
    cpdef bint empty(self)
    cpdef bint verify(self)
    cpdef void clear(self)
    cpdef NodeStack to_stack(self)
    cpdef sparsekey size(self)
    cpdef int8_t get(self, sparsekey key, int8_t default=*)
    cpdef int8_t find(self, sparsekey key)
    cdef NodeStack path_to_root(self, sparsekey key)
    cdef NodeStack path_to_node(self, sparsekey key)
    cpdef void insert(self, sparsekey key, int8_t value=*)
    cdef void rebalance_node(self, IntTreeNode* node)
    cpdef void delete(self, sparsekey key, bint silent=*)
    cdef void delleaf(self, IntTreeNode* node, NodeStack ancestors)
    cdef void del1childl(self, IntTreeNode* node, NodeStack ancestors)
    cdef void del1childr(self, IntTreeNode* node, NodeStack ancestors)
    cdef void del2child(self, IntTreeNode* node, NodeStack ancestors)
    cpdef void delrange(self, sparsekey start, sparsekey end)
    cpdef IntTree getrange(self, sparsekey start, sparsekey end)
    cpdef IntTree intersection(self, IntTree other)
    cpdef IntTree union(self, IntTree other)

cdef bint node_verify(IntTreeNode* node)
cdef int8_t node_balance(IntTreeNode* node)
cdef void update_node_height(IntTreeNode* node)

cdef void deltree(IntTreeNode* start)

cdef void rotate_right(IntTreeNode* root)
cdef void rotate_left(IntTreeNode* root)
cdef void rotate_double_left(IntTreeNode* root)
cdef void rotate_double_right(IntTreeNode* root)

cdef struct NodeStackItem:
    IntTreeNode* node
    NodeStackItem* following

cdef class NodeStack(object):
    cdef NodeStackItem* front
    cdef void push(self, IntTreeNode* node)
    cdef IntTreeNode* peek(self)
    cdef IntTreeNode* pop(self)
    cdef bint empty(self)

