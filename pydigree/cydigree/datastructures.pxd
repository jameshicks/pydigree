from libc.stdint cimport uint32_t, int8_t, uint8_t

cdef inline bint compare(object a, object b, int op)

cdef class SparseArray:
    cdef readonly IntTree container
    cdef readonly int8_t refcode
    cdef readonly uint32_t size
    cdef inline uint32_t fix_index(self, uint32_t index)
    cdef object _get_single_item(self, uint32_t index)
    cdef SparseArray _get_slice(self, index)
    cdef _get_fancyidx(self, index)
    cdef _get_boolidx(self, index)
    cpdef void set_item(self, uint32_t index, int8_t value)
    cdef void _set_slice(self, uint32_t start, uint32_t stop, values)
    cdef void _set_slice_to_sparray(self, uint32_t start, uint32_t stop, SparseArray values)
    cdef void _set_boolidx(self, indices, values)
    cdef void _set_fancyidx(self, indices, value)
    cpdef _cmp_single(self, object value, int op)
    cpdef SparseArray _cmp_sequence(self, other, int op)
    cpdef SparseArray _cmp_sparray(self, SparseArray other, int op)
    cpdef bint any(self)
    cpdef bint all(self)
    cpdef SparseArray logical_not(self)
    cpdef double sparsity(self)

cdef class IntTreeNode(object):
    cdef readonly uint32_t key
    cdef public int8_t value
    cdef uint8_t height
    cdef public IntTreeNode left, right, parent
    cpdef bint is_leaf(self)
    cdef update_height(self)
    cpdef int8_t balance(self)
    cpdef bint is_balanced(self)

cdef class IntTree(object):
    cdef readonly IntTreeNode root
    cpdef bint empty(self)
    cpdef NodeStack to_stack(self)
    cpdef uint32_t size(self)
    cpdef int8_t get(self, uint32_t key, int8_t default=*)
    cpdef int8_t find(self, uint32_t key)
    cpdef IntTreeNode find_node(self, uint32_t key)
    cpdef NodeStack path_to_root(self, uint32_t key)
    cpdef NodeStack path_to_node(self, uint32_t key)
    cpdef void insert(self, uint32_t key, int8_t value=*)
    cdef rebalance_node(self, IntTreeNode node)
    cpdef void delete(self, uint32_t key, bint silent=*)
    cdef void delleaf(self, IntTreeNode node, NodeStack ancestors)
    cdef void del1childl(self, IntTreeNode node, NodeStack ancestors)
    cdef void del1childr(self, IntTreeNode node, NodeStack ancestors)
    cdef void del2child(self, IntTreeNode node, NodeStack ancestors)
    cpdef IntTreeNode min_node(self, start=*)
    cpdef uint32_t min(self)
    cpdef IntTreeNode max_node(self, start=*)
    cpdef delrange(self, uint32_t start, uint32_t end)
    cpdef getrange(self, uint32_t start, uint32_t end)
    cpdef uint32_t max(self)
    cpdef IntTree intersection(self, IntTree other)
    cpdef IntTree union(self, IntTree other)

cpdef IntTreeNode rotate_right(IntTreeNode root)
cpdef IntTreeNode rotate_left(IntTreeNode root)
cpdef IntTreeNode rotate_double_left(IntTreeNode root)
cpdef IntTreeNode rotate_double_right(IntTreeNode root)

cdef class NodeStack(object):
    cdef public NodeStackItem front
    cpdef void push(self, IntTreeNode val)
    cpdef IntTreeNode pop(self)
    cpdef IntTreeNode peek(self)
    cpdef bint empty(self)
    cpdef void reverse(self)
cdef class NodeStackItem(object):
    cdef NodeStackItem following
    cdef IntTreeNode obj
