from collections import Sequence
from libc.stdint cimport uint32_t, uint8_t, int8_t
from libc.stdio cimport printf
from cpython.mem cimport PyMem_Malloc, PyMem_Free

import numpy as np

DEF MAX_HEIGHT=30

cdef inline bint sparse_val_compare(sparse_val a, sparse_val b, int op):
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

    def __init__(self, sparse_key size, refcode=None, initial=None):
        self.container = IntTree()
        self.refcode = refcode
        self.size = size

    def __len__(self):
        return self.size

    def copy(self):
        cdef SparseArray output = SparseArray.from_items(self.items(), 
                                                         self.size, 
                                                         self.refcode)
        return output

    property values:
        def __get__(self):
            return list(self.container.values())


    property indices:
        def __get__(self):
            return list(self.container.keys())

    def keys(self):
        return list(self.container.keys())

    def values(self):
        return list(self.container.values())

    cdef inline sparse_key fix_index(self, sparse_key index): 
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

    cdef object _get_single_item(self, sparse_key index):
        index = self.fix_index(index)
        return self.container.get(index, self.refcode)

    cdef SparseArray _get_slice(self, index):
        cdef sparse_key start = index.start
        cdef sparse_key stop = index.stop 

        cdef SparseArray subarray = SparseArray(stop - start, self.refcode)
        subarray.container = self.container.getrange(start, stop)
        cdef NodeStack ns = subarray.container.to_stack()
        cdef IntTreeNode* node = ns.pop()
        while node:
            node.key = node.key - start 
            node = ns.pop()
        return subarray

    cdef _get_fancyidx(self, index):
        cdef sparse_key i = 0
        cdef sparse_key nvals = len(index)
        cdef output = SparseArray(nvals, self.refcode)
        for i in range(nvals):
            output[i] = self[index[i]]

        return output

    cdef _get_boolidx(self, index):
        if len(index) != self.size:
            raise IndexError('Bool indices not same size as array')
        cdef sparse_key nvals = sum(index)
        cdef sparse_key i = 0
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

    cpdef void set_item(self, sparse_key index, sparse_val value):
        index = self.fix_index(index)
        if value != self.refcode:
            self.container.insert(index, value)
        else:
            if not self.container.empty():
                self.container.delete(index)

    cdef void _set_slice(self, sparse_key start, sparse_key stop, values):
        cdef sparse_key i, nvals 

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
                ival = <sparse_val>values[i - start]
            else:
                ival = <sparse_val>values
            
            if ival != self.refcode:
                self.container.insert(i, ival) 

    cdef void _set_slice_to_sparray(self, sparse_key start, sparse_key stop, SparseArray values):
        if stop - start != values.size:
            raise IndexError('Value wrong shape for slice')
        
        toclear = [idx for idx in self.keys() if start <= idx <= stop]
        for idx in toclear:
            self.container.delete(idx)

        # TODO: Make efficient instead of using python traversal
        for k, v in values.items():
            self.container.insert(k + start, v)

    cdef void _set_boolidx(self, indices, values):
        cdef list trueidxs = [i for i,x in enumerate(indices) if x]
        self._set_fancyidx(trueidxs, values)

    cdef void _set_fancyidx(self, indices, value):
        cdef bint multivalue = isinstance(value, Sequence) and not isinstance(value, str)
        cdef sparse_key nidx = len(indices)
        cdef sparse_key nval 
        
        if multivalue:
            nval = len(value)
            if nidx != nval:
                raise IndexError('Indices and values different length')

        cdef sparse_key i
        for i in range(nidx):
            if multivalue:
                self[<sparse_key>indices[i]] = <sparse_val>value[i] 
            else:
                self[<sparse_key>indices[i]] = <sparse_val>value

    # Comparison methods
    #    
    cpdef _cmp_single(self, object value, int op):
        cdef SparseArray output = SparseArray(self.size, sparse_val_compare(self.refcode, value, op))
        cdef NodeStack s = self.container.to_stack()

        cdef IntTreeNode* node = s.pop()
        while node:
            output.set_item(node.key, sparse_val_compare(node.value, value, op))
            node = s.pop()

        return output

    cpdef SparseArray _cmp_sequence(self, other, int op):
        cdef SparseArray output = SparseArray(self.size, False)
        cdef sparse_key seqlen = len(other)
        cdef sparse_key i = 0
        cdef NodeStack densesites = self.container.to_stack()
        cdef IntTreeNode* curdense = densesites.pop()
        
        for i in range(seqlen):
            if curdense != NULL and i == curdense.key:
                output.set_item(i,  sparse_val_compare(curdense.value, other[i], op))
                curdense = densesites.pop()
            else:
                output.set_item(i,  sparse_val_compare(self.refcode, other[i], op))

        return output

    cpdef SparseArray _cmp_sparray(self, SparseArray other, int op):
        if self.size != other.size:
            raise IndexError('Cannot sparse_val_compare arrays of differing size')

        cdef NodeStack selfstack = self.container.to_stack()
        cdef NodeStack otherstack = other.container.to_stack()

        cdef IntTreeNode* selfnode = selfstack.pop()
        cdef IntTreeNode* othernode = otherstack.pop()

        cdef SparseArray output = SparseArray(self.size, self.refcode == other.refcode)
        while selfnode != NULL and othernode != NULL:
            
            if othernode == NULL or selfnode.key < othernode.key:
                output.set_item(selfnode.key, sparse_val_compare(selfnode.value, other.refcode, op))
                selfnode = selfstack.pop()
            elif selfnode == NULL or selfnode.key > othernode.key:
                output.set_item(othernode.key, sparse_val_compare(self.refcode, othernode.value, op)) 
                othernode = otherstack.pop() 
            else:
                output.set_item(selfnode.key, sparse_val_compare(selfnode.value, othernode.value, op))
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
        
        cdef sparse_val value
        for value in self.container.values():
            if value:
                return True
        
        return False

    cpdef bint all(self):
        if self.sparsity() > 0 and not self.refcode:
            return False
        
        cdef sparse_val value
        for value in self.container.values():
            if not value:
                return False
        
        return True

    cpdef SparseArray logical_not(self):
        cdef SparseArray output = SparseArray(self.size, not self.refcode)
        cdef NodeStack s = self.container.to_stack()
        cdef IntTreeNode* node = s.pop()

        while node:
            output.set_item(node.key, not node.value)
            node = s.pop()
        return output

    # Builders
    @staticmethod
    def from_dense(dense, refcode):
        cdef sparse_key i = 0
        cdef sparse_key n = len(dense)
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
                output.set_item(i, denseview[i])

        return output

    @staticmethod
    def from_items(items, size, refcode):
        cdef SparseArray output = SparseArray(size, refcode)
        for k,v in items:
            output.set_item(k, v)

        return output

    # Misc
    #
    cpdef double sparsity(self):
        'Returns the proportion of sparse sites in the array'
        return 1 - <double>self.container.size() / self.size

    def tolist(self):
        cdef list output = [self.refcode] * self.size
        # TODO: Make efficient instead of using python traversal
        for k, v in zip(self.container.keys(), self.container.values()):
            output[k] = v

        return output

    def items(self):
        # TODO: Make efficient instead of using python traversal
        for k, v in zip(self.container.keys(), self.container.values()):
            yield k,v

############
############

cdef IntTreeNode* new_node(sparse_key key, sparse_val value):
    cdef IntTreeNode* node = <IntTreeNode*>PyMem_Malloc(sizeof(IntTreeNode))
    if not node:
        raise MemoryError("Couldnt alloc memory for node")

    node.key = key
    node.value = value
    node.left = NULL
    node.right = NULL
    node.height = 0

    return node

cdef void del_node(IntTreeNode* node):
    PyMem_Free(node)

############
############

cdef class IntTree(object):
    def __init__(self):
        self.root = NULL

    def __dealloc__(self):
        self.clear()

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

    def __contains__(self, sparse_key key):
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

    cpdef bint verify(self):
        if self.empty():
            return True

        return node_verify(self.root)

    cpdef bint empty(self):
        return self.root == NULL

    cpdef void clear(self):
        'Removes all nodes from tree'
        if not self.root:
            return 
        deltree(self.root)
        self.root = NULL

    def traverse(self):
        if self.empty():
            return

        cdef NodeStack s = NodeStack()
        cdef IntTreeNode* node = self.root
        
        while (not s.empty()) or (node != NULL):
            if node != NULL:
                s.push(node)
                node = node.left
            else:
                node = s.pop()
                yield node.key 
                node = node.right

    def keys(self):
        if self.empty():
              return

        cdef NodeStack s = NodeStack()
        cdef IntTreeNode* node = self.root
        
        while (not s.empty()) or (node != NULL):
            if node != NULL:
                s.push(node)
                node = node.left
            else:
                node = s.pop()
                yield node.key 
                node = node.right

    def values(self):
        if self.empty():
            return

        cdef NodeStack s = NodeStack()
        cdef IntTreeNode* node = self.root
        
        while (not s.empty()) or (node != NULL):
            if node != NULL:
                s.push(node)
                node = node.left
            else:
                node = s.pop()
                yield node.value 
                node = node.right

    cpdef NodeStack to_stack(self):
        cdef NodeStack s = NodeStack()
        cdef NodeStack out = NodeStack()
        cdef IntTreeNode* node = self.root

        if not node:
            return s

        while (not s.empty()) or (node != NULL):
            if node != NULL:
                s.push(node)
                node = node.right
            else:
                node = s.pop()
                out.push(node)
                node = node.left
        return out

    cpdef sparse_key size(self):
        if self.empty(): 
            return 0

        cdef NodeStack s = NodeStack()
        cdef sparse_key tot = 0
        cdef IntTreeNode* node = self.root

        while (not s.empty()) or (node != NULL):
            if node != NULL:
                s.push(node)
                node = node.right
            else:
                node = s.pop()
                tot += 1 
                node = node.left
        return tot

    cpdef sparse_val get(self, sparse_key key, sparse_val default=0):
        if self.empty():
            return default

        cdef IntTreeNode* node = self.root
        while node != NULL:
            if key > node.key:
                node = node.right
            elif key < node.key:
                node = node.left
            else:
                return node.value

        return default

    cpdef sparse_val find(self, sparse_key key):
        if self.empty():
            raise KeyError('Node not found')

        cdef IntTreeNode* node = self.root
        while node:
            if key > node.key:
                node = node.right
            elif key < node.key:
                node = node.left
            else:
                return node.value

        raise KeyError('Node not found')

    cpdef void insert(self, sparse_key key, sparse_val value=0):
        cdef IntTreeNode* inserted = new_node(key, value)
        if self.empty():
            self.root = inserted
            return

        cdef IntTreeNode* cur_node = self.root
        cdef IntTreeNode* stack[MAX_HEIGHT]
        cdef uint8_t depth = 0

        while True:
            if inserted.key > cur_node.key:
                if cur_node.right != NULL:
                    stack[depth] = cur_node
                    depth += 1 

                    cur_node = cur_node.right
                else:
                    cur_node.right = inserted
                    break

            elif inserted.key < cur_node.key:
                if cur_node.left != NULL:
                    stack[depth] = cur_node
                    depth += 1 

                    cur_node = cur_node.left
                else:
                    cur_node.left = inserted
                    break
            else:
                cur_node.value = value
                # We don't need to rebalance if the key was already in the tree
                del_node(inserted)
                return 

        stack[depth] = inserted
        cdef IntTreeNode* ancestor = inserted
        cdef IntTreeNode* ancestor_ancestor = NULL
        while depth:
            depth -= 1
            ancestor = stack[depth]
            ancestor_ancestor = stack[depth-1] if depth > 0 else NULL
            self.rebalance_node(ancestor, ancestor_ancestor)

    cdef void rebalance_node(self, IntTreeNode* node, IntTreeNode* parent):
        update_node_height(node)

        cdef int8_t balance = node_balance(node)
        
        if -1 <= balance <= 1:
            return

        elif balance > 1:
            if node.right != NULL and node_balance(node.right) < 0:
                new_root = node.right.left
                rotate_double_left(node, parent)

            else:
                new_root = node.right
                rotate_left(node, parent)

        elif balance < -1:
            if node.left != NULL and node_balance(node.left) > 0:
                new_root = node.left.right
                rotate_double_right(node, parent)

            else:
                new_root = node.left
                rotate_right(node, parent)

        else:
            pass
        if node is self.root:
            self.root = new_root


    cpdef void delete(self, sparse_key key, bint silent=True):
        cdef IntTreeNode* stack[MAX_HEIGHT]
        cdef uint8_t depth = 0

        cdef IntTreeNode* node = self.root
        while node:

            if key > node.key:
                stack[depth] = node
                depth += 1
                node = node.right
            
            elif key < node.key:
                stack[depth] = node
                depth += 1
                node = node.left
            else:
                stack[depth] = node
                break
        else:
            if not silent:
                raise KeyError('Node not found') 
        
        cdef IntTreeNode* parent = stack[depth - 1] if depth > 0 else NULL

        if node.right == NULL and node.left == NULL: # LEAF 
            if not parent:
                self.root = NULL
            elif node.key > parent.key:
                parent.right = NULL
            else:
                parent.left = NULL

            del_node(node)

        elif node.right == NULL: # 1 Child on left
            if not parent: # node is the root
                self.root = node.left
            elif node.key > parent.key:
                parent.right = node.left
            else:
                parent.left = node.left

            del_node(node)

        elif node.left == NULL: # 1 Child on right
            if not parent:
                self.root = node.right
            elif node.key > parent.key: 
                parent.right = node.right
            else:
                parent.left = node.right

 
            del_node(node)

        else:
            self.del2child(node)
            return
        
        cdef IntTreeNode* ancestor
        cdef IntTreeNode* ancestor_ancestor
        while depth:
            depth -= 1
            ancestor = stack[depth]
            ancestor_ancestor = stack[depth-1] if depth > 0 else NULL
            self.rebalance_node(ancestor, ancestor_ancestor) 

    cdef void del2child(self, IntTreeNode* node):

        cdef IntTreeNode* replacement = node.left
        while replacement.right:
            replacement = replacement.right

        cdef sparse_key key = replacement.key
        cdef sparse_val value = replacement.value

        self.delete(key)

        node.key = key
        node.value = value


    cpdef void delrange(self, sparse_key start, sparse_key end):
        'Deletes keys where start <= key < stop'
        cdef NodeStack delstack = NodeStack()
        cdef NodeStack s = NodeStack()
        cdef IntTreeNode* node = self.root

        while (not s.empty()) or (node != NULL):
            if node != NULL:
                s.push(node)
                node = node.left
            else:
                node = s.pop()
                if start <= node.key < end:
                    delstack.push(node)
                elif node.key >= end:
                    break 
                node = node.right

        cdef IntTreeNode* delnode = delstack.pop() 
        while delnode:
            self.delete(delnode.key)
            delnode = delstack.pop()


    cpdef IntTree getrange(self, sparse_key start, sparse_key end):
        cdef IntTree ntree = IntTree()
        cdef NodeStack s = NodeStack()
        cdef IntTreeNode* node = self.root

        while (not s.empty()) or (node != NULL):
            if node:
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

    cpdef IntTree intersection(self, IntTree other):
        cdef NodeStack t1 = self.to_stack()
        cdef NodeStack t2 = other.to_stack()

        cdef IntTreeNode* a = t1.pop()
        cdef IntTreeNode* b = t2.pop()

        cdef IntTree ntree = IntTree()
        
        while a != NULL and b != NULL:
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

        cdef IntTreeNode* a = t1.pop()
        cdef IntTreeNode* b = t2.pop()

        cdef IntTree ntree = IntTree()

        while a != NULL or b != NULL:
            if a == NULL:
                ntree.insert(b.key)
                b = t2.pop()

            elif b == NULL:
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
cdef bint node_verify(IntTreeNode* node):
    if not -1 <= node_balance(node) <= 1:
        return False
    if node.left != NULL and not (node.key > node.left.key):
        return False
    if node.right != NULL and not (node.key < node.right.key):
        return False

    cdef bint l = node_verify(node.left) if node.left else True 
    cdef bint r = node_verify(node.right) if node.right else True
    return l and r

cdef int8_t node_balance(IntTreeNode* node):
    cdef uint8_t lefth = (node.left.height) if node.left != NULL else 0
    cdef uint8_t righth = (node.right.height) if node.right != NULL else 0
    cdef int8_t balance = righth - lefth
    return balance

cdef void update_node_height(IntTreeNode* node):
        lheight = node.left.height if node.left != NULL else 0
        rheight = node.right.height if node.right != NULL else 0
        node.height = max(lheight, rheight) + 1


cdef void deltree(IntTreeNode* node):
    """
    Recursively deletes all nodes in tree a node's subtree. Faster than 
    deleting each node individually because it does not have to rebalance the 
    tree after each deletion.
    """ 
    if node == NULL:
        return

    # Traverse in post-order removing nodes
    # TODO: rewrite in iterative style
    deltree(node.left)
    deltree(node.right)

    node.left = NULL
    node.right = NULL
    del_node(node)

cdef void rotate_right(IntTreeNode* root, IntTreeNode* parent):
    cdef IntTreeNode* pivot = root.left
    root.left = pivot.right
    pivot.right = root

    update_node_height(root)
    update_node_height(pivot)

    if not parent:
        return
    elif parent.key > root.key:
        parent.left = pivot
    else:
        parent.right = pivot

    update_node_height(parent)

cdef void rotate_left(IntTreeNode* root, IntTreeNode* parent):
    cdef IntTreeNode* pivot = root.right
    root.right = pivot.left
    pivot.left = root 

    update_node_height(root)
    update_node_height(pivot)

    if not parent:
        return
    elif parent.key > root.key:
        parent.left = pivot
    else:
        parent.right = pivot

    update_node_height(parent)


cdef void rotate_double_left(IntTreeNode* root, IntTreeNode* parent):
    rotate_right(root.right, root)
    rotate_left(root, parent)


cdef void rotate_double_right(IntTreeNode* root, IntTreeNode* parent):
    rotate_left(root.left, root)
    rotate_right(root, parent)



#########


cdef class NodeStack:
    def __cinit__(self):
        self.front = NULL

    cdef void push(self, IntTreeNode* node):
        cdef NodeStackItem* item = <NodeStackItem*>PyMem_Malloc(sizeof(NodeStackItem))
        
        item.node = node
        item.following = self.front
        
        self.front = item

    cdef IntTreeNode* peek(self):
        if self.front:
            return self.front.node
        else:
            return NULL

    cdef IntTreeNode* pop(self):
        if not self.front:
            return NULL
        
        item = self.front 
        self.front = item.following
        
        cdef IntTreeNode* node = item.node
        PyMem_Free(item)
        return node 

    cdef bint empty(self):
        return self.front == NULL


def print_sizes():
    print('IntTreeNode: {}'.format(sizeof(IntTreeNode)))
    print('IntTree: {}'.format(sizeof(IntTree)))
    print('SparseArray: {}'.format(sizeof(SparseArray)))
    print('NodeStackItem: {}'.format(sizeof(NodeStackItem)))
    print('NodeStack: {}'.format(sizeof(NodeStack)))