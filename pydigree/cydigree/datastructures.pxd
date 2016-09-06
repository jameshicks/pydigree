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

cdef class IntTreeNode(object):
    cdef readonly sparsekey key
    cdef public int8_t value
    cdef uint8_t height
    cdef public IntTreeNode left, right, parent

cdef class IntTree(object):
    cdef readonly IntTreeNode root
    cpdef bint empty(self)
    cpdef NodeStack to_stack(self)
    cpdef sparsekey size(self)
    cpdef int8_t get(self, sparsekey key, int8_t default=*)
    cpdef int8_t find(self, sparsekey key)
    cpdef IntTreeNode find_node(self, sparsekey key)
    cpdef NodeStack path_to_root(self, sparsekey key)
    cpdef NodeStack path_to_node(self, sparsekey key)
    cpdef void insert(self, sparsekey key, int8_t value=*)
    cdef rebalance_node(self, IntTreeNode node)
    cpdef void delete(self, sparsekey key, bint silent=*)
    cdef void delleaf(self, IntTreeNode node, NodeStack ancestors)
    cdef void del1childl(self, IntTreeNode node, NodeStack ancestors)
    cdef void del1childr(self, IntTreeNode node, NodeStack ancestors)
    cdef void del2child(self, IntTreeNode node, NodeStack ancestors)
    cpdef IntTreeNode min_node(self, start=*)
    cpdef sparsekey min(self)
    cpdef IntTreeNode max_node(self, start=*)
    cpdef delrange(self, sparsekey start, sparsekey end)
    cpdef getrange(self, sparsekey start, sparsekey end)
    cpdef sparsekey max(self)
    cpdef IntTree intersection(self, IntTree other)
    cpdef IntTree union(self, IntTree other)

cpdef int8_t node_balance(IntTreeNode node)
cpdef void update_node_height(IntTreeNode node)

cpdef void rotate_right(IntTreeNode root)
cpdef void rotate_left(IntTreeNode root)
cpdef void rotate_double_left(IntTreeNode root)
cpdef void rotate_double_right(IntTreeNode root)

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
