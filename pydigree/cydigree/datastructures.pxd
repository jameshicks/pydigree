from libc.stdint cimport uint32_t, int8_t, uint8_t

ctypedef uint32_t sparse_key
ctypedef int8_t sparse_val 

cdef inline bint sparse_val_compare(sparse_val a, sparse_val b, int op)

cdef class SparseArray:
    cdef readonly IntTree container
    cdef readonly int8_t refcode
    cdef readonly sparse_key size
    cdef inline sparse_key fix_index(self, sparse_key index)
    cdef object _get_single_item(self, sparse_key index)
    cdef SparseArray _get_slice(self, index)
    cdef _get_fancyidx(self, index)
    cdef _get_boolidx(self, index)
    cpdef void set_item(self, sparse_key index, sparse_val value)
    cdef void _set_slice(self, sparse_key start, sparse_key stop, values)
    cdef void _set_slice_to_sparray(self, sparse_key start, sparse_key stop, SparseArray values)
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
    sparse_key key
    sparse_val value 
    int8_t height

cdef class IntTree(object):
    cdef IntTreeNode* root
    cpdef bint empty(self)
    cpdef bint verify(self)
    cpdef void clear(self)
    cpdef NodeStack to_stack(self)
    cpdef sparse_key size(self)
    cpdef sparse_val get(self, sparse_key key, sparse_val default=*)
    cpdef sparse_val find(self, sparse_key key)
    cpdef void insert(self, sparse_key key, sparse_val value=*)
    cdef void rebalance_node(self, IntTreeNode* node, IntTreeNode* parent)
    cpdef void delete(self, sparse_key key, bint silent=*)
    cdef void del2child(self, IntTreeNode* node)
    cpdef void delrange(self, sparse_key start, sparse_key end)
    cpdef IntTree getrange(self, sparse_key start, sparse_key end)
    cpdef IntTree intersection(self, IntTree other)
    cpdef IntTree union(self, IntTree other)

cdef bint node_verify(IntTreeNode* node)
cdef int8_t node_balance(IntTreeNode* node)
cdef void update_node_height(IntTreeNode* node)

cdef void deltree(IntTreeNode* start)

cdef void rotate_right(IntTreeNode* root, IntTreeNode* parent)
cdef void rotate_left(IntTreeNode* root, IntTreeNode* parent)
cdef void rotate_double_left(IntTreeNode* root, IntTreeNode* parent)
cdef void rotate_double_right(IntTreeNode* root, IntTreeNode* parent)

cdef struct NodeStackItem:
    IntTreeNode* node
    NodeStackItem* following

cdef class NodeStack(object):
    cdef NodeStackItem* front
    cdef void push(self, IntTreeNode* node)
    cdef IntTreeNode* peek(self)
    cdef IntTreeNode* pop(self)
    cdef bint empty(self)

