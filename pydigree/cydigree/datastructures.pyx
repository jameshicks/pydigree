from collections import Sequence

from libc.stdint cimport uint32_t, uint8_t, int8_t

cdef class SparseArray:
    # I imagine at some point I'll have to rewrite this whole thing
    # as a red-black tree or something, but until then we'll just bsearch
    # a sorted list

    cdef readonly IntTree container
    cdef readonly object refcode
    cdef readonly Py_ssize_t size

    def __init__(self, size, refcode=None, initial=None):
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

    cdef inline Py_ssize_t fix_index(self, Py_ssize_t index): 
        if index >= 0:
            return index 
        else:
            return self.size + index

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self._get_slice(index)
        elif isinstance(index, Sequence) and type(index[0]) is bool:
            return self._get_boolidx(index)
        elif isinstance(index, Sequence) and type(index[0]) is int:
            return self._get_fancyidx(index)
        else:
            return self._get_single_item(index)

    cdef object _get_single_item(self, index):
        index = self.fix_index(index)
        return self.container.get(index, self.refcode)

    cdef SparseArray _get_slice(self, index):
        cdef uint32_t start = index.start
        cdef uint32_t stop = index.stop 

        cdef SparseArray subarray = SparseArray(stop - start, self.refcode)
        subarray.container = self.container.getrange(start, stop)

        return subarray

    cdef _get_fancyidx(self, index):
        cdef Py_ssize_t i = 0
        cdef Py_ssize_t nvals = len(index)
        cdef output = SparseArray(nvals, self.refcode)
        for i in range(nvals):
            output[i] = self[index[i]]

        return output

    cdef _get_boolidx(self, index):
        if len(index) != self.size:
            raise IndexError('Bool indices not same size as array')
        cdef uint32_t nvals = sum(index)
        cdef uint32_t i = 0
        cdef output = SparseArray(nvals, self.refcode)
        cdef bint included
        

        for i in range(nvals):
            included = index[i]
            if included:
                output[i] = self.container.get(i, self.refcode)
        return output

    def __setitem__(self, index, object value):
        if isinstance(index, slice):
            self._set_slice(index.start, index.stop, value)
        else:
            self._set_single(index, value)

    cdef void _set_single(self, uint32_t index, object value):
        index = self.fix_index(index)
        if value != self.refcode:
            self.container.insert(index, value)
        else:
            self.container.delete(index)

    cdef void _set_slice(self, uint32_t start, uint32_t stop, values):
        cdef Py_ssize_t i, nvals 

        if isinstance(values, SparseArray):
            self._set_slice_to_sparray(start, stop, values)
        elif isinstance(values, Sequence):
            nvals = len(values)
            if stop - start != nvals:
                raise IndexError('Value wrong shape for slice')

        self.container.delrange(start, stop)
        for i in range(start, stop):
            if isinstance(values, Sequence):
                ival = values[i - start]
            else:
                ival = values
            if ival != self.refcode:
                self.container.insert(i, ival) 

    cdef void _set_slice_to_sparray(self, uint32_t start, uint32_t stop, SparseArray values):
        if stop - start != values.size:
            raise IndexError('Value wrong shape for slice')
        
        for node in values.container.traverse():
            self.container.insert(node.key, node.value)

    def tolist(self):
        cdef list output = [self.refcode] * self.size
        for node in self.container.traverse():
            output[node.key] = node.value

        return output 

############
############


cdef class IntTreeNode(object):
    'An IntTree node'
    cdef readonly uint32_t key
    cdef public object value
    cdef uint8_t height
    cdef public IntTreeNode left, right, parent

    def __init__(self, uint32_t key, value=None):
        self.key = key
        self.value = value
        self.height = 0
        self.left = None
        self.right = None
        self.parent = None

    def __repr__(self):
        return 'IntTreeNode({})'.format(self.key)

    cpdef bint is_leaf(self):
        return self.left is None and self.right is None

    cdef update_height(self):
        lheight = self.left.height if self.left is not None else 0
        rheight = self.right.height if self.right is not None else 0
        self.height = max(lheight, rheight) + 1

    @property
    def children(self):
        return self.left, self.right

    @children.setter
    def children(self, value):
        self.left, self.right = value

    cpdef int8_t balance(self):
        cdef uint8_t lefth = (self.left.height) if self.left is not None else 0
        cdef uint8_t righth = (self.right.height) if self.right is not None else 0
        return righth - lefth

    cpdef bint is_balanced(self):
        return -1 <= self.balance() <= 1

cdef class IntTree(object):
    cdef readonly IntTreeNode root
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

    def __contains__(self, uint32_t key):
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
        return self.root is not None

    def __len__(self):
        return self.size()

    def traverse(self, reverse=False):
        if not self.root:
            return

        if reverse:
            yield from self._traverse_reverse()
        else:
            yield from self._traverse()

    def _traverse(self):
        s = NodeStack()
        node = self.root
        while (not s.empty()) or (node is not None):
            if node is not None:
                s.push(node)
                node = node.left
            else:
                node = s.pop()
                yield node 
                node = node.right

    def _traverse_reverse(self):
        s = NodeStack()
        node = self.root
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

    cpdef uint32_t size(self):
        if self.root is None: 
            return 0
        cdef NodeStack s = NodeStack()
        cdef uint32_t tot = 0
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

    cpdef object get(self, uint32_t key, default=None):
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

    cpdef object find(self, uint32_t key):
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


    cpdef IntTreeNode find_node(self, uint32_t key):
        cdef IntTreeNode node = self.root
        while node is not None:

            if node.key == key:
                return node

            elif key < node.key:
                node = node.left
            elif key > node.key:
                node = node.right
        raise KeyError('Key not found: {}'.format(key))

    cpdef NodeStack path_to_root(self, uint32_t key):
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
        raise KeyError('Node not found'.format(key))

    cpdef NodeStack path_to_node(self, uint32_t key):
        s = self.path_to_root(key)
        s.reverse()
        return s

    cpdef void insert(self, uint32_t key, object value=None):
        cdef IntTreeNode new_node = IntTreeNode(key, value)

        if self.root is None:
            self.root = new_node
            return

        cur_node = self.root
        
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

        parent = new_node

        while parent:
            self.rebalance_node(parent)
            parent = parent.parent

    cdef rebalance_node(self, IntTreeNode node):
        node.update_height()

        if node.is_balanced():
            return

        if node.balance() > 1:
            if node.right is not None and node.right.balance() < 0:
                new_root = node.right.left
                rotate_double_left(node)

            else:
                new_root = node.right
                rotate_left(node)

        elif node.balance() < -1:
            if node.left is not None and node.left.balance() > 0:
                new_root = node.left.right
                rotate_double_right(node)

            else:
                new_root = node.left
                rotate_right(node)

        else:
            pass
        if node is self.root:
            self.root = new_root


    cpdef void delete(self, uint32_t key):
        if self.root is None:
            raise KeyError('Tree is empty')

        ancestors = self.path_to_root(key)
        node = ancestors.pop()


        if node.is_leaf():
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
            node.left.update_height()
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
            node.right.update_height()
        
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

        replacement.children = node.children 
        replacement.parent = direct_ancestor
        replacement.update_height()

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

    cpdef uint32_t min(self):
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

    cpdef delrange(self, uint32_t start, uint32_t end):
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

    cpdef getrange(self, uint32_t start, uint32_t end):
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
                     

    cpdef uint32_t max(self):
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


cpdef IntTreeNode rotate_right(IntTreeNode root):
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

    root.update_height()
    pivot.update_height()


cpdef IntTreeNode rotate_left(IntTreeNode root):

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
    root.update_height()
    pivot.update_height()


cpdef IntTreeNode rotate_double_left(IntTreeNode root):
    rotate_right(root.right)
    rotate_left(root)


cpdef IntTreeNode rotate_double_right(IntTreeNode root):
    rotate_left(root.left)
    rotate_right(root)



#########



cdef class NodeStack(object):
    'A first-in last-out datastructure for IntTree Nodes'
    cdef public NodeStackItem front
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
    cdef NodeStackItem following
    cdef IntTreeNode obj

    def __init__(self, IntTreeNode obj, following=None):
        self.obj = obj
        self.following = following