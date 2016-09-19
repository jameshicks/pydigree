from libc.stdio cimport printf
from libc.string cimport memcpy
from libc.stdint cimport uint32_t, uint8_t, int8_t
from cpython.mem cimport PyMem_Malloc, PyMem_Free

DEF BINSIZE=16
DEF MAX_DEPTH = 32

ctypedef  uint32_t variantkey

cdef struct VariantTreeNode:
    VariantTreeNode* left
    VariantTreeNode* right
    uint32_t key
    uint8_t height
    int8_t values[BINSIZE]

cdef VariantTreeNode* new_node(uint32_t major_key, int8_t refcode):
    cdef VariantTreeNode* node = <VariantTreeNode*>PyMem_Malloc(sizeof(VariantTreeNode))
    if not node:
        printf("Cant allocate memory for node\n")
    node.key = major_key
    node.left = NULL
    node.right = NULL
    node.height = 0

    cdef int i
    for i in range(BINSIZE):
        node.values[i] = refcode
    return node

cdef void delete_node(VariantTreeNode* node):
    PyMem_Free(node)

cdef class VariantTree:
    cdef readonly int8_t refcode
    cdef VariantTreeNode* root

    def __init__(self, int8_t refcode=0):
        self.refcode = refcode
        self.root = NULL

    def __dealloc__(self):
        deltree(self.root)

    cpdef bint empty(self):
        return self.root == NULL

    cpdef bint verify(self):
        if self.empty():
            return True

        return node_verify(self.root)

    cpdef int _binsize(self):
        return BINSIZE

    cpdef keys(self):

        cdef int minor_key

        cdef NodeStack s = NodeStack()
        cdef VariantTreeNode* node = self.root
        
        outp = []
        while (not s.empty()) or (node != NULL):
            if node != NULL:
                s.push(node)
                node = node.left
            else:
                node = s.pop()
                for minor_key in range(BINSIZE):
                    if node.values[minor_key] != self.refcode:
                        outp.append(node.key * BINSIZE  + minor_key)
                node = node.right

        return outp

    cpdef values(self):

        cdef int minor_key

        cdef NodeStack s = NodeStack()
        cdef VariantTreeNode* node = self.root
        
        outp = []
        while (not s.empty()) or (node != NULL):
            if node != NULL:
                s.push(node)
                node = node.left
            else:
                node = s.pop()
                for minor_key in range(BINSIZE):
                    if node.values[minor_key] != self.refcode:
                        outp.append(node.values[minor_key])
                node = node.right

        return outp

    cpdef node_count(self):
        return treesize(self.root)

    cpdef int8_t get_item(self, variantkey key):
        cdef uint32_t major_key = key // BINSIZE
        cdef int8_t minor_key = key % BINSIZE       

        cdef VariantTreeNode* node = self.root
        while node:

            if node.key > major_key:
                node = node.left
            elif node.key < major_key:
                node = node.right
            else:
                return node.values[minor_key]
        return self.refcode

    cpdef void set_item(self, variantkey key, int8_t value):
        if value == self.refcode:
            self.clear_item(key)
            return 
        
        cdef VariantTreeNode* tbi
        cdef uint32_t major_key = key // BINSIZE
        cdef int8_t minor_key = key % BINSIZE
        
        if not self.root:
            tbi = new_node(major_key, self.refcode)
            tbi.values[minor_key] = value
            self.root = tbi
            return 

        cdef VariantTreeNode* stack[MAX_DEPTH]  

        cdef int depth = 0
        cdef VariantTreeNode* cur_node = self.root
        while cur_node:
            stack[depth] = cur_node
            depth += 1

            if cur_node.key < major_key:
                cur_node = cur_node.right
            elif cur_node.key > major_key:
                cur_node = cur_node.left
            else:
                cur_node.values[minor_key] = value
                return

        cur_node = stack[depth-1]
        tbi = new_node(major_key, self.refcode)
        tbi.values[minor_key] = value
        
        if major_key > cur_node.key:
            cur_node.right = tbi
        else:
            cur_node.left = tbi


        cdef VariantTreeNode* node
        cdef VariantTreeNode* nodes_parent
        while depth:
            depth -= 1
            node = stack[depth]
            nodes_parent = stack[depth-1] if depth > 0 else NULL
            self.rebalance_node(node, nodes_parent)

    cpdef void clear_item(self, variantkey key):
        cdef uint32_t major_key = key // BINSIZE
        cdef int8_t minor_key = key % BINSIZE
        
        cdef VariantTreeNode* cur_node = self.root
        while cur_node:
            if cur_node.key > major_key:
                cur_node = cur_node.left
            elif cur_node.key < major_key:
                cur_node = cur_node.right
            else:
                break
        else:
            return

        cur_node.values[minor_key] = self.refcode

        # Check if node is empty
        cdef int i
        for i in range(BINSIZE):
            # if node is empty, remove it
            if cur_node.values[i] != self.refcode:
                return
        self._delnode(major_key)

    cdef void rebalance_node(self, VariantTreeNode* node, VariantTreeNode* parent):
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

    
    cdef void _delnode(self, variantkey major_key):
        cdef VariantTreeNode* stack[MAX_DEPTH]
        cdef int depth = 0

        cdef VariantTreeNode* node = self.root

        while node:
            stack[depth] = node
            depth += 1
            if node.key > major_key:
                node = node.left
            elif node.key < major_key:
                node = node.right
            else:
                break
        else:
            return
        
        depth -= 1 
        
        cdef VariantTreeNode* par = stack[depth-1] if depth else NULL 
        
        if node.left and node.right: # 2 children
            self._del2child(node)
            return
        
        elif node.left: # Only 1 child, left
            if not par:
                self.root = node.left
            elif node.key > par.key:
                par.left = node.left
            else:
                par.right = node.left

        elif node.right: # Only 1 child, right
            if not par:
                self.root = node.right
            elif node.key > par.key:
                par.left = node.right
            else:
                par.right = node.right
        
        else: # Leaf node
            if not par:
                self.root = NULL
            elif node.key > par.key:
                par.right = NULL
            else:
                par.left = NULL 

        delete_node(node)

        while depth:
            depth -= 1
            node = stack[depth]
            par = stack[depth-1] if depth else NULL
            self.rebalance_node(node, par) 

    cdef void _del2child(self, VariantTreeNode* node):
        cdef VariantTreeNode* sucessor = node.right
        while sucessor:
            sucessor = sucessor.left

        cdef uint32_t tmpkey = sucessor.key
        memcpy(&node.values, &sucessor.values, BINSIZE*sizeof(int8_t))

        self._delnode(sucessor.key)
        node.key = tmpkey

    cpdef VariantTree getrange(self, uint32_t start, uint32_t stop):
        cdef uint32_t start_major_key = start // BINSIZE
        cdef uint32_t stop_major_key = stop // BINSIZE

        cdef uint32_t start_minor_key = start % BINSIZE
        cdef uint32_t stop_minor_key = stop % BINSIZE

        cdef NodeStack s = NodeStack()
        cdef VariantTreeNode* node = self.root

        cdef uint32_t localstart, localstop, minor_key, fullidx

        cdef VariantTree outp = VariantTree()

        while (not s.empty()) or (node != NULL):
            if node != NULL:
                s.push(node)
                node = node.left
            else:
                node = s.pop()
                if start_major_key <= node.key <= stop_major_key:
                    localstart = start_minor_key if node.key == start_major_key else 0
                    localstop = stop_minor_key if node.key == stop_minor_key else BINSIZE

                    for minor_key in range(localstart, localstop):
                        fullidx = node.key * BINSIZE  + minor_key
                        outp.set_item(fullidx, node.values[minor_key])
                
                elif node.key > stop_major_key:
                    break

                node = node.right

        return outp

    cpdef void clearrange(self, uint32_t start, uint32_t stop):
        pass

    #cdef setrange(self):
    #    pass

    #cdef intersection(self):
    #    pass

    #cdef union(self):
    #    pass

# Node manipulation functions
cdef bint node_verify(VariantTreeNode* node):
    if not -1 <= node_balance(node) <= 1:
        return False
    if node.left != NULL and not (node.key > node.left.key):
        return False
    if node.right != NULL and not (node.key < node.right.key):
        return False

    cdef bint l = node_verify(node.left) if node.left else True 
    cdef bint r = node_verify(node.right) if node.right else True
    return l and r

cdef int8_t node_balance(VariantTreeNode* node):
    cdef uint8_t lefth = (node.left.height) if node.left != NULL else 0
    cdef uint8_t righth = (node.right.height) if node.right != NULL else 0
    cdef int8_t balance = righth - lefth
    return balance

cdef void update_node_height(VariantTreeNode* node):
        lheight = node.left.height if node.left != NULL else 0
        rheight = node.right.height if node.right != NULL else 0
        node.height = max(lheight, rheight) + 1

cdef uint32_t treesize(VariantTreeNode* node):
    if node == NULL:
        return 0
    return 1 + treesize(node.left) + treesize(node.right)


cdef void deltree(VariantTreeNode* node):
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
    delete_node(node)

cdef void rotate_right(VariantTreeNode* root, VariantTreeNode* parent):
    cdef VariantTreeNode* pivot = root.left
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

cdef void rotate_left(VariantTreeNode* root, VariantTreeNode* parent):
    cdef VariantTreeNode* pivot = root.right
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


cdef void rotate_double_left(VariantTreeNode* root, VariantTreeNode* parent):
    rotate_right(root.right, root)
    rotate_left(root, parent)


cdef void rotate_double_right(VariantTreeNode* root, VariantTreeNode* parent):
    rotate_left(root.left, root)
    rotate_right(root, parent)

###
cdef struct NodeStackItem:
    VariantTreeNode* node
    NodeStackItem* following

cdef class NodeStack:
    cdef NodeStackItem* front
    def __cinit__(self):
        self.front = NULL

    cdef void push(self, VariantTreeNode* node):
        cdef NodeStackItem* item = <NodeStackItem*>PyMem_Malloc(sizeof(NodeStackItem))
        
        item.node = node
        item.following = self.front
        
        self.front = item

    cdef VariantTreeNode* peek(self):
        if self.front:
            return self.front.node
        else:
            return NULL

    cdef VariantTreeNode* pop(self):
        if not self.front:
            return NULL
        
        item = self.front 
        self.front = item.following
        
        cdef VariantTreeNode* node = item.node
        PyMem_Free(item)
        return node 

    cdef bint empty(self):
        return self.front == NULL
