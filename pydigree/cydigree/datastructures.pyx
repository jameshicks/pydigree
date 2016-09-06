from collections import Sequence
from libc.stdint cimport uint32_t, uint8_t, int8_t
from cpython.mem cimport PyMem_Malloc, PyMem_Free

import numpy as np

cdef inline bint compare(object a, object b, int op):
    # <   0
    # ==  2
    # >   4
    # <=  1
    # !=  3
    # >=  5

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

    def __init__(self, sparsekey size, refcode=None, initial=None):
        self.container = IntTree()
        self.refcode = refcode
        self.size = size

    def __len__(self):
        return self.size

    property values:
        def __get__(self):
            return [x.value for x in self.container.traverse()]


    property indices:
        def __get__(self):
            return [x.key for x in self.container.traverse()]

    cdef inline sparsekey fix_index(self, sparsekey index): 
        if index >= 0:
            return index 
        else:
            return self.size + index

    # Item Getting
    #

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self._get_slice(index)
        elif isinstance(index, Sequence) and type(index[0]) is bool:
            return self._get_boolidx(index)
        elif isinstance(index, Sequence) and type(index[0]) is int:
            return self._get_fancyidx(index)
        else:
            return self._get_single_item(index)

    cdef object _get_single_item(self, sparsekey index):
        index = self.fix_index(index)
        return self.container.get(index, self.refcode)

    cdef SparseArray _get_slice(self, index):
        cdef sparsekey start = index.start
        cdef sparsekey stop = index.stop 

        cdef SparseArray subarray = SparseArray(stop - start, self.refcode)
        subarray.container = self.container.getrange(start, stop)
        cdef NodeStack ns = subarray.container.to_stack()
        cdef IntTreeNode node = ns.pop()
        while node:
            node.key = node.key - start 
            node = ns.pop()
        return subarray

    cdef _get_fancyidx(self, index):
        cdef sparsekey i = 0
        cdef sparsekey nvals = len(index)
        cdef output = SparseArray(nvals, self.refcode)
        for i in range(nvals):
            output[i] = self[index[i]]

        return output

    cdef _get_boolidx(self, index):
        if len(index) != self.size:
            raise IndexError('Bool indices not same size as array')
        cdef sparsekey nvals = sum(index)
        cdef sparsekey i = 0
        cdef output = SparseArray(nvals, self.refcode)
        cdef bint included
        

        for i in range(nvals):
            included = index[i]
            if included:
                output[i] = self.container.get(i, self.refcode)
        return output

    # Item Setting
    #

    def __setitem__(self, index, object value):
        if isinstance(index, slice):
            self._set_slice(index.start, index.stop, value)
        elif isinstance(index, Sequence) and type(index[0]) is bool:
            self._set_boolidx(index, value)
        elif isinstance(index, Sequence) and type(index[0]) is int:
            self._set_fancyidx(index, value)
        else:
            self.set_item(index, value)

    cpdef void set_item(self, sparsekey index, int8_t value):
        index = self.fix_index(index)
        if value != self.refcode:
            self.container.insert(index, value)
        else:
            if not self.container.empty():
                self.container.delete(index)

    cdef void _set_slice(self, sparsekey start, sparsekey stop, values):
        cdef sparsekey i, nvals 

        if isinstance(values, SparseArray):
            self._set_slice_to_sparray(start, stop, values)
            return
        
        elif isinstance(values, Sequence) and not isinstance(values, str):
            nvals = len(values)
            if stop - start != nvals:
                raise IndexError('Value wrong shape for slice')

        self.container.delrange(start, stop)
        for i in range(start, stop):
            if isinstance(values, Sequence) and not isinstance(values, str):
                ival = values[i - start]
            else:
                ival = values
            if ival != self.refcode:
                self.container.insert(i, ival) 

    cdef void _set_slice_to_sparray(self, sparsekey start, sparsekey stop, SparseArray values):
        if stop - start != values.size:
            raise IndexError('Value wrong shape for slice')
        
        for node in values.container.traverse():
            self.container.insert(node.key + start, node.value)

    cdef void _set_boolidx(self, indices, values):
        cdef list trueidxs = [i for i,x in enumerate(indices) if x]
        self._set_fancyidx(trueidxs, values)

    cdef void _set_fancyidx(self, indices, value):
        cdef bint multivalue = isinstance(value, Sequence) and not isinstance(value, str)
        cdef sparsekey nidx = len(indices)
        cdef sparsekey nval 
        
        if multivalue:
            nval = len(value)
            if nidx != nval:
                raise IndexError('Indices and values different length')

        cdef sparsekey i
        for i in range(nidx):
            if multivalue:
                self[indices[i]] = value[i] 
            else:
                self[indices[i]] = value

    # Comparison methods
    #    
    cpdef _cmp_single(self, object value, int op):
        cdef SparseArray output = SparseArray(self.size, compare(self.refcode, value, op))
        cdef NodeStack s = self.container.to_stack()

        cdef IntTreeNode node = s.pop()
        while node:
            output[node.key] = compare(node.value, value, op)
            node = s.pop()

        return output

    cpdef SparseArray _cmp_sequence(self, other, int op):
        cdef SparseArray output = SparseArray(self.size, False)
        cdef sparsekey seqlen = len(other)
        cdef sparsekey i = 0
        cdef NodeStack densesites = self.container.to_stack()
        cdef IntTreeNode curdense = densesites.pop()
        
        for i in range(seqlen):
            if curdense is not None and i == curdense.key:
                output[i] = compare(curdense.value, other[i], op)
                curdense = densesites.pop()
            else:
                output[i] = compare(self.refcode, other[i], op)

        return output

    cpdef SparseArray _cmp_sparray(self, SparseArray other, int op):
        if self.size != other.size:
            raise IndexError('Cannot compare arrays of differing size')

        cdef NodeStack selfstack = self.container.to_stack()
        cdef NodeStack otherstack = other.container.to_stack()

        cdef IntTreeNode selfnode = selfstack.pop()
        cdef IntTreeNode othernode = otherstack.pop()

        cdef SparseArray output = SparseArray(self.size, self.refcode == other.refcode)
        while selfnode is not None and othernode is not None:
            
            if othernode is None or selfnode.key < othernode.key:
                output[selfnode.key] = compare(selfnode.value, other.refcode, op)
                selfnode = selfstack.pop()
            elif selfnode is None or selfnode.key > othernode.key:
                output[othernode.key] = compare(self.refcode, othernode.value, op) 
                othernode = otherstack.pop() 
            else:
                output[selfnode.key] = compare(selfnode.value, othernode.value, op)
                selfnode = selfstack.pop()
                othernode = otherstack.pop()

        return output

    def __richcmp__(self, value, op):
        if type(value) is SparseArray:
            return self._cmp_sparray(value, op)
        elif isinstance(value, Sequence) and not isinstance(value, str):
            return self._cmp_sequence(value, op)
        else:
            return self._cmp_single(value, op)
    
    # Logic functions
    #
    cpdef bint any(self):
        if self.refcode:
            return True
        
        for node in self.container.traverse():
            if node.value:
                return True
        
        return False

    cpdef bint all(self):
        if self.sparsity() > 0 and not self.refcode:
            return False
        
        for node in self.container.traverse():
            if not node.value:
                return False
        
        return True

    cpdef SparseArray logical_not(self):
        cdef SparseArray output = SparseArray(self.size, not self.refcode)
        cdef NodeStack s = self.container.to_stack()
        cdef IntTreeNode node = s.pop()

        while node:
            output[node.key] = not node.value
            node = s.pop()
        return output

    # Builders
    @staticmethod
    def from_dense(dense, refcode):
        cdef sparsekey i = 0
        cdef sparsekey n = len(dense)
        cdef object x
        
        cdef SparseArray output = SparseArray(n, refcode)

        for i in range(n):
            x = dense[i]
            if x != refcode:
                output[i] = x

        return output

    @staticmethod
    def from_numpy(dense, int8_t refcode):
        if not type(dense) is np.ndarray:
            raise ValueError("Dense data is not numpy array")
            
        dense8 = dense.astype(np.int8)
        cdef int8_t [:] denseview = dense8

        cdef int i = 0
        cdef int size = dense8.shape[0]

        cdef SparseArray output = SparseArray(size, refcode)

        for i in range(size):
            if denseview[i] != refcode:
                output[i] = denseview[i]

        return output

    # Misc
    #
    cpdef double sparsity(self):
        'Returns the proportion of sparse sites in the array'
        return 1 - <double>self.container.size() / self.size

    def tolist(self):
        cdef list output = [self.refcode] * self.size
        for node in self.container.traverse():
            output[node.key] = node.value

        return output

    def items(self):
        for node in self.container.traverse():
            yield node.key, node.value

############
############

cdef struct IntTreeNode2:
    IntTreeNode2* left 
    IntTreeNode2* right
    IntTreeNode2* parent
    sparsekey key
    int8_t value 
    int8_t height

cdef IntTreeNode2* new_node(sparsekey key, int8_t value):
    cdef IntTreeNode2* node = <IntTreeNode2*>PyMem_Malloc(sizeof(IntTreeNode2))
    if not node:
        raise MemoryError("Couldnt alloc memory for node")

    node.key = key
    node.value = value
    node.left = NULL
    node.right = NULL
    node.height = 0

    return node

cdef void del_node(IntTreeNode2* node):
    PyMem_Free(node)

############
############

cdef class IntTreeNode(object):
    'An IntTree node'

    def __init__(self, sparsekey key, int8_t value=0):
        self.key = key
        self.value = value
        self.height = 0
        self.left = None
        self.right = None
        self.parent = None

    def __repr__(self):
        return 'IntTreeNode({})'.format(self.key)

cdef class IntTree(object):
    def __init__(self):
        self.root = None

    @staticmethod
    def from_keys(keys):
        tree = IntTree()
        for key in keys:
            tree.insert(key)

        return tree

    @staticmethod
    def from_pairs(pairs):
        tree = IntTree()
        for key, value in pairs:
            tree.insert(key, value)

        return tree

    def __contains__(self, sparsekey key):
        node = self.root
        while node:
            if key > node.key:
                node = node.right
            elif key < node.key:
                node = node.left
            elif key == node.key:
                return True

        return False

    def __nonzero__(self):
        return not self.empty()

    def __len__(self):
        return self.size()

    cpdef bint empty(self):
        return self.root is None

    def traverse(self, reverse=False):
        if not self.root:
            return

        if reverse:
            yield from self._traverse_reverse()
        else:
            yield from self._traverse()

    def _traverse(self):
        cdef NodeStack s = NodeStack()
        cdef IntTreeNode node = self.root
        
        while (not s.empty()) or (node is not None):
            if node is not None:
                s.push(node)
                node = node.left
            else:
                node = s.pop()
                yield node 
                node = node.right

    def _traverse_reverse(self):
        cdef NodeStack s = NodeStack()
        cdef IntTreeNode node = self.root
        while (not s.empty()) or (node is not None):
            if node is not None:
                s.push(node)
                node = node.right
            else:
                node = s.pop()
                yield node 
                node = node.left

    cpdef NodeStack to_stack(self):
        cdef NodeStack s = NodeStack()
        cdef NodeStack out = NodeStack()
        cdef IntTreeNode node = self.root

        while (not s.empty()) or (node is not None):
            if node is not None:
                s.push(node)
                node = node.right
            else:
                node = s.pop()
                out.push(node)
                node = node.left
        return out

    cpdef sparsekey size(self):
        if self.root is None: 
            return 0
        cdef NodeStack s = NodeStack()
        cdef sparsekey tot = 0
        cdef IntTreeNode node = self.root

        while (not s.empty()) or (node is not None):
            if node is not None:
                s.push(node)
                node = node.right
            else:
                node = s.pop()
                tot += 1 
                node = node.left
        return tot

    def keys(self):
        yield from (node.key for node in self.traverse())

    cpdef int8_t get(self, sparsekey key, int8_t default=0):
        if not self.root:
            return default

        cdef IntTreeNode node = self.root
        while node is not None:
            if key > node.key:
                node = node.right
            elif key < node.key:
                node = node.left
            else:
                return node.value

        return default

    cpdef int8_t find(self, sparsekey key):
        if not self.root:
            raise KeyError('Node not found')

        cdef IntTreeNode node = self.root
        while node is not None:
            if key > node.key:
                node = node.right
            elif key < node.key:
                node = node.left
            else:
                return node.value

        raise KeyError('Node not found')


    cpdef IntTreeNode find_node(self, sparsekey key):
        cdef IntTreeNode node = self.root
        while node is not None:

            if node.key == key:
                return node

            elif key < node.key:
                node = node.left
            elif key > node.key:
                node = node.right
        raise KeyError('Key not found: {}'.format(key))

    cpdef NodeStack path_to_root(self, sparsekey key):
        s = NodeStack()
        cur_node = self.root
        while cur_node:
            s.push(cur_node)
            if key > cur_node.key:
                cur_node = cur_node.right
            elif key < cur_node.key:
                cur_node = cur_node.left
            elif key == cur_node.key:
                return s
        raise KeyError('Node not found {}'.format(key))

    cpdef NodeStack path_to_node(self, sparsekey key):
        s = self.path_to_root(key)
        s.reverse()
        return s

    cpdef void insert(self, sparsekey key, int8_t value=0):
        cdef IntTreeNode new_node = IntTreeNode(key, value)

        if self.root is None:
            self.root = new_node
            return

        cdef IntTreeNode cur_node = self.root
        
        while True:
            if new_node.key > cur_node.key:
                if cur_node.right is not None:
                    cur_node = cur_node.right
                else:
                    cur_node.right = new_node
                    new_node.parent = cur_node
                    break

            elif new_node.key < cur_node.key:
                if cur_node.left is not None:
                    cur_node = cur_node.left
                else:
                    cur_node.left = new_node
                    new_node.parent = cur_node
                    break
            else:
                cur_node.value = value
                # We don't need to rebalance if the key was already in the tree
                return 

        cdef IntTreeNode parent = new_node

        while parent:
            self.rebalance_node(parent)
            parent = parent.parent

    cdef rebalance_node(self, IntTreeNode node):
        update_node_height(node)

        cdef int8_t balance = node_balance(node)
        
        if -1 <= balance <= 1:
            return

        elif balance > 1:
            if node.right is not None and node_balance(node.right) < 0:
                new_root = node.right.left
                rotate_double_left(node)

            else:
                new_root = node.right
                rotate_left(node)

        elif balance < -1:
            if node.left is not None and node_balance(node.left) > 0:
                new_root = node.left.right
                rotate_double_right(node)

            else:
                new_root = node.left
                rotate_right(node)

        else:
            pass
        if node is self.root:
            self.root = new_root


    cpdef void delete(self, sparsekey key, bint silent=True):
        if self.root is None:
            raise KeyError('Tree is empty')

        try:
            ancestors = self.path_to_root(key)
        except KeyError:
            if silent:
                return
            else:
                raise
        node = ancestors.pop()


        if node.right is None and node.left is None:
            self.delleaf(node, ancestors)
        elif node.right is None:
            self.del1childl(node, ancestors)
        elif node.left is None:
            self.del1childr(node, ancestors)
        else:
            self.del2child(node, ancestors)

    cdef void delleaf(self, IntTreeNode node, NodeStack ancestors):
        cdef IntTreeNode ancestor = ancestors.peek()
        if node is self.root:
            self.root = None
            return

        elif node.key > ancestor.key:
            ancestor.right = None
        else: 
            ancestor.left = None

        for ancestor in ancestors:
            self.rebalance_node(ancestor)

    cdef void del1childl(self, IntTreeNode node, NodeStack ancestors):
        if self.root is node:
            self.root = node.left
            update_node_height(node.left)

            return

        cdef IntTreeNode ancestor = ancestors.peek()
        cdef IntTreeNode child = node.left
        child.parent = ancestor
        if child.key > ancestor.key:
            ancestor.right = child
        else:
            ancestor.left = child

        for ancestor in ancestors:
            self.rebalance_node(ancestor)

    cdef void del1childr(self, IntTreeNode node, NodeStack ancestors):
        if self.root is node:
            self.root = node.right
            update_node_height(node.right)
        
        cdef IntTreeNode ancestor = ancestors.peek()
        cdef IntTreeNode child = node.right    
        child.parent = ancestor
        if child.key > ancestor.key:
            ancestor.right = child
        else:
            ancestor.left = child

        for ancestor in ancestors:
            self.rebalance_node(ancestor)

    cdef void del2child(self, IntTreeNode node, NodeStack ancestors):
        # Find a replacement for the node to be deleted
        cdef IntTreeNode replacement = node.left
        while replacement.right:
            replacement = replacement.right

        path_to_replacement = self.path_to_root(replacement.key)
        path_to_replacement.pop() # Remove replacement from stack

        self.delete(replacement.key)

        # The parent of the node to be deleted
        cdef IntTreeNode direct_ancestor = ancestors.peek()

        replacement.left, replacement.right = node.left, node.right 
        replacement.parent = direct_ancestor
        update_node_height(replacement)

        if node is self.root:
            self.root = replacement
        
        elif replacement.key > direct_ancestor.key:
            direct_ancestor.right = replacement
        else:
            direct_ancestor.left = replacement

        for ancestor in path_to_replacement:
            self.rebalance_node(ancestor)

    cpdef IntTreeNode min_node(self, start=None):
        if not self.root:
            raise KeyError('Tree empty!')
        if not start:
            node = self.root
        else:
            node = start

        while node.left is not None:
            node = node.left

        return node

    cpdef sparsekey min(self):
        return self.min_node().key

    cpdef IntTreeNode max_node(self, start=None):
        if not self.root:
            raise KeyError('Tree empty!')
        if not start:
            node = self.root
        else:
            node = start

        while node.right is not None:
            node = node.right

        return node

    cpdef delrange(self, sparsekey start, sparsekey end):
        'Deletes keys where start <= key < stop'
        cdef NodeStack delstack = NodeStack()
        cdef NodeStack s = NodeStack()
        cdef IntTreeNode node = self.root

        while (not s.empty()) or (node is not None):
            if node is not None:
                s.push(node)
                node = node.left
            else:
                node = s.pop()
                if start <= node.key < end:
                    delstack.push(node)
                elif node.key >= end:
                    break 
                node = node.right

        for delnode in delstack:
            self.delete(delnode.key)

    cpdef getrange(self, sparsekey start, sparsekey end):
        cdef IntTree ntree = IntTree()
        cdef NodeStack s = NodeStack()
        cdef IntTreeNode node = self.root

        while (not s.empty()) or (node is not None):
            if node is not None:
                s.push(node)
                node = node.left
            else:
                node = s.pop()
                if start <= node.key < end:
                    ntree.insert(node.key, node.value)
                elif node.key >= end:
                    break 
                node = node.right

        return ntree
                     

    cpdef sparsekey max(self):
        return self.max_node().key

    cpdef IntTree intersection(self, IntTree other):
        cdef NodeStack t1 = self.to_stack()
        cdef NodeStack t2 = other.to_stack()

        a, b = t1.pop(), t2.pop()

        cdef IntTree ntree = IntTree()
        while a is not None and b is not None:
            if a.key == b.key:
                ntree.insert(a.key)
                a, b = t1.pop(), t2.pop()
            elif a.key < b.key:
                a = t1.pop()
            else:
                b = t2.pop()

        return ntree

    cpdef IntTree union(self, IntTree other):
        cdef NodeStack t1 = self.to_stack()
        cdef NodeStack t2 = other.to_stack()

        a, b = t1.pop(), t2.pop()

        ntree = IntTree()

        while a is not None or b is not None:
            if a is None:
                ntree.insert(b.key)
                b = t2.pop()

            elif b is None:
                ntree.insert(a.key)
                a = t1.pop()

            elif a.key == b.key:
                ntree.insert(a.key)
                a, b = t1.pop(), t2.pop()

            elif a.key < b.key:
                ntree.insert(a.key)
                a = t1.pop()

            else:
                ntree.insert(b.key)
                b = t2.pop()

        return ntree

# Node manipulation functions
cpdef int8_t node_balance(IntTreeNode node):
    cdef uint8_t lefth = (node.left.height) if node.left is not None else 0
    cdef uint8_t righth = (node.right.height) if node.right is not None else 0
    cdef int8_t balance = righth - lefth
    return balance

cpdef void update_node_height(IntTreeNode node):
        lheight = node.left.height if node.left is not None else 0
        rheight = node.right.height if node.right is not None else 0
        node.height = max(lheight, rheight) + 1

cpdef void rotate_right(IntTreeNode root):
    pivot = root.left
    root.left = pivot.right
    if root.left is not None:
        root.left.parent = root

    pivot.right = root
    pivot.parent = root.parent
    root.parent = pivot

    if pivot.parent is None:
        pass
    elif pivot.parent.left is root:
        pivot.parent.left = pivot
    elif pivot.parent.right is root:
        pivot.parent.right = pivot

    update_node_height(root)
    update_node_height(pivot)


cpdef void rotate_left(IntTreeNode root):

    pivot = root.right
    root.right = pivot.left

    if root.right is not None:
        root.right.parent = root

    pivot.left = root
    pivot.parent = root.parent
    root.parent = pivot

    if pivot.parent is None:
        pass
    elif pivot.parent.left is root:
        pivot.parent.left = pivot
    elif pivot.parent.right is root:
        pivot.parent.right = pivot
    else:
        raise KeyError
    update_node_height(root)
    update_node_height(pivot)


cpdef void rotate_double_left(IntTreeNode root):
    rotate_right(root.right)
    rotate_left(root)


cpdef void rotate_double_right(IntTreeNode root):
    rotate_left(root.left)
    rotate_right(root)



#########



cdef class NodeStack(object):
    'A first-in last-out datastructure for IntTree Nodes'
    def __init__(self, starts=None):
        self.front = None
        if starts:
            for val in starts:
                self.push(val)

    def __bool__(self):
        return self.front is not None

    cpdef void push(self, IntTreeNode val):
        cdef NodeStackItem item = NodeStackItem(val, following=self.front)
        self.front = item

    cpdef IntTreeNode pop(self):
        popped = self.front
        if popped is None:
            return None

        self.front = popped.following
        return popped.obj

    cpdef IntTreeNode peek(self):
        if self.front:
            return self.front.obj
        else:
            return None

    cpdef bint empty(self):
        return self.front is None

    cpdef void reverse(self):
        l = list(self)
        for x in l:
            self.push(x)

    def __iter__(self):
        item = self.pop()
        while item is not None:
            yield item
            item = self.pop()

    def _to_list(self):
        l = []
        f = self.front
        while f:
            l.append(f.obj) 
            f = f.following

        return l


cdef class NodeStackItem(object):
    def __init__(self, IntTreeNode obj, following=None):
        self.obj = obj
        self.following = following

def print_sizes():
    print('IntTreeNode: {}'.format(sizeof(IntTreeNode)))
    print('IntTree: {}'.format(sizeof(IntTree)))
    print('SparseArray: {}'.format(sizeof(SparseArray)))
    print('NodeStackItem: {}'.format(sizeof(NodeStackItem)))
    print('NodeStack: {}'.format(sizeof(NodeStack)))