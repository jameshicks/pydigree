from libc.stdio cimport printf
from libc.string cimport memcpy
from libc.stdint cimport uint32_t, uint8_t, int8_t
from cpython.mem cimport PyMem_Malloc, PyMem_Free

DEF BINSIZE=16
DEF MAX_DEPTH = 32

ctypedef  uint32_t variantkey

cdef struct VariantTreeNode:
    VariantTreeNode* link[2]
    uint32_t key
    uint8_t height
    int8_t values[BINSIZE]

cdef VariantTreeNode* new_node(uint32_t major_key, int8_t refcode):
    cdef VariantTreeNode* node = <VariantTreeNode*>PyMem_Malloc(sizeof(VariantTreeNode))
    if not node:
        printf("Cant allocate memory for node\n")
    node.key = major_key
    node.link[0] = NULL
    node.link[1] = NULL
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

    def __contains__(self, variantkey key):
        'Returns true if a variant is in the tree'
        cdef uint32_t major_key = key // BINSIZE
        cdef uint32_t minor_key = key % BINSIZE

        if self.empty():
            return False

        cdef VariantTreeNode* node = self.root
        while node:
            if node.key == major_key:
                return node.values[minor_key] != self.refcode
            else:
                node = node.link[node.key < major_key]

        return False
        
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
                node = node.link[0]
            else:
                node = s.pop()
                for minor_key in range(BINSIZE):
                    if node.values[minor_key] != self.refcode:
                        outp.append(node.key * BINSIZE  + minor_key)
                node = node.link[1]

        return outp

    cpdef values(self):

        cdef int minor_key

        cdef NodeStack s = NodeStack()
        cdef VariantTreeNode* node = self.root
        
        outp = []
        while (not s.empty()) or (node != NULL):
            if node != NULL:
                s.push(node)
                node = node.link[0]
            else:
                node = s.pop()
                for minor_key in range(BINSIZE):
                    if node.values[minor_key] != self.refcode:
                        outp.append(node.values[minor_key])
                node = node.link[1]

        return outp

    cpdef node_count(self):
        return treesize(self.root)

    cpdef int8_t get_item(self, variantkey key):
        cdef uint32_t major_key = key // BINSIZE
        cdef int8_t minor_key = key % BINSIZE       

        cdef VariantTreeNode* node = self.root
        while node:
            if node.key == major_key:
                return node.values[minor_key]
            else:
                node = node.link[node.key < major_key]
            
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

            if cur_node.key == major_key:
                cur_node.values[minor_key] = value
                return
            else:
                cur_node = cur_node.link[cur_node.key < major_key]

        cur_node = stack[depth-1]
        tbi = new_node(major_key, self.refcode)
        tbi.values[minor_key] = value
        
        if major_key > cur_node.key:
            cur_node.link[1] = tbi
        else:
            cur_node.link[0] = tbi


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

            if cur_node.key == major_key:
                break
            else:
                cur_node = cur_node.link[cur_node.key < major_key]

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
            if node.link[1] != NULL and node_balance(node.link[1]) < 0:
                new_root = node.link[1].link[0]
                rotate_double_left(node, parent)

            else:
                new_root = node.link[1]
                rotate_left(node, parent)

        elif balance < -1:
            if node.link[0] != NULL and node_balance(node.link[0]) > 0:
                new_root = node.link[0].link[1]
                rotate_double_right(node, parent)

            else:
                new_root = node.link[0]
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

            if node.key == major_key:
                break
            else:
                node = node.link[node.key < major_key]

        else:
            return
        
        depth -= 1 
        
        cdef VariantTreeNode* par = stack[depth-1] if depth else NULL
        cdef VariantTreeNode* target 
        
        if node.link[0] and node.link[1]: # 2 children
            self._del2child(node, par)
            return
        
        target = NULL if node_is_leaf(node) else node.link[node.link[0] == NULL]
   
        if not par: # Node is root
            self.root = target
        else:
            par.link[par.key < node.key] = target 

        delete_node(node)

        while depth:
            depth -= 1
            node = stack[depth]
            par = stack[depth-1] if depth else NULL
            self.rebalance_node(node, par) 

    cdef void _del2child(self, VariantTreeNode* node, VariantTreeNode* par):
        # Find successor (smallest in the right subtree)
        cdef VariantTreeNode* successor = node.link[1]
        while successor:
            if not successor.link[0]:
                break
            successor = successor.link[0]

        # Get the path to the successor from the root
        cdef VariantTreeNode* path[MAX_DEPTH]
        cdef int depth = 0
        cdef VariantTreeNode* cur_node = self.root
        while cur_node:
            path[depth] = cur_node
            depth += 1

            if cur_node == successor:
                depth -= 1
                break
            else:
                cur_node = cur_node.link[cur_node.key < successor.key]

        # Remove the successor node
        cdef VariantTreeNode* successor_par = path[depth-1]

        cdef VariantTreeNode* target 
        if node_is_leaf(successor):
            target = NULL
        else:
            target = successor.link[successor.link[0] == NULL] 

        successor_par.link[successor_par.key < successor.key] = target

        # Place the successor where the node was
        successor.link[0] = node.link[0]
        successor.link[1] = node.link[1]
        if par == NULL:
            self.root = successor
        else:
            par.link[par.key < node.key] = successor

        cdef int i = 0
        for i in range(depth):
            if path[i] == node:
                path[i] = successor
                break

        # Rebalance nodes all the way back up tree
        cdef VariantTreeNode* cur_node_par
        while depth:
            depth -= 1
            cur_node = path[depth]
            cur_node_par = path[depth-1] if depth else NULL
            self.rebalance_node(cur_node, cur_node_par)

        delete_node(node)

    cpdef VariantTree get_range(self, uint32_t start, uint32_t stop):
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
                node = node.link[0]
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

                node = node.link[1]

        return outp

    cpdef void clear_range(self, uint32_t start, uint32_t stop):
        cdef uint32_t start_major_key = start // BINSIZE
        cdef uint32_t stop_major_key = stop // BINSIZE

        cdef uint32_t start_minor_key = start % BINSIZE
        cdef uint32_t stop_minor_key = stop % BINSIZE

        cdef NodeStack s = NodeStack()
        cdef VariantTreeNode* node = self.root

        cdef uint32_t localstart, localstop, minor_key, fullidx

        cdef NodeStack delstack = NodeStack()


        while (not s.empty()) or (node != NULL):
            if node != NULL:
                s.push(node)
                node = node.link[0]
            else:
                node = s.pop()

                if node.key < start_major_key:
                    pass
                elif node.key > stop_major_key:
                    break
                elif start_major_key < node.key < stop_major_key:
                    delstack.push(node)
                else:
                    localstart = start_minor_key if node.key == start_major_key else 0
                    localstop = stop_minor_key if node.key == stop_major_key else BINSIZE

                    for minor_key in range(localstart, localstop):
                        node.values[minor_key] = self.refcode

                    # Delte the node if it's empty
                    for minor_key in range(0, BINSIZE):
                        if node.values[minor_key] != self.refcode:
                            break
                    else:
                        delstack.push(node)

                node = node.link[1]

        cdef VariantTreeNode* delnode = delstack.pop() 
        while delnode:
            self._delnode(delnode.key)
            delnode = delstack.pop()

    #cdef setrange(self):
    #    pass

    #cdef intersection(self):
    #    pass

    #cdef union(self):
    #    pass

# Node manipulation functions
cdef bint node_is_leaf(VariantTreeNode* node):
    return node.link[0] == NULL and node.link[1] == NULL

cdef bint node_verify(VariantTreeNode* node):
    if not -1 <= node_balance(node) <= 1:
        return False
    if node.link[0] != NULL and not (node.key > node.link[0].key):
        return False
    if node.link[1] != NULL and not (node.key < node.link[1].key):
        return False

    cdef bint l = node_verify(node.link[0]) if node.link[0] else True 
    cdef bint r = node_verify(node.link[1]) if node.link[1] else True
    return l and r

cdef int8_t node_balance(VariantTreeNode* node):
    cdef uint8_t lefth = (node.link[0].height) if node.link[0] != NULL else 0
    cdef uint8_t righth = (node.link[1].height) if node.link[1] != NULL else 0
    cdef int8_t balance = righth - lefth
    return balance

cdef void update_node_height(VariantTreeNode* node):
        lheight = node.link[0].height if node.link[0] != NULL else 0
        rheight = node.link[1].height if node.link[1] != NULL else 0
        node.height = max(lheight, rheight) + 1

cdef uint32_t treesize(VariantTreeNode* node):
    if node == NULL:
        return 0
    return 1 + treesize(node.link[0]) + treesize(node.link[1])


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
    deltree(node.link[0])
    deltree(node.link[1])

    node.link[0] = NULL
    node.link[1] = NULL
    delete_node(node)

cdef void rotate_right(VariantTreeNode* root, VariantTreeNode* parent):
    cdef VariantTreeNode* pivot = root.link[0]
    root.link[0] = pivot.link[1]
    pivot.link[1] = root

    update_node_height(root)
    update_node_height(pivot)

    if not parent:
        return
    elif parent.key > root.key:
        parent.link[0] = pivot
    else:
        parent.link[1] = pivot

    update_node_height(parent)

cdef void rotate_left(VariantTreeNode* root, VariantTreeNode* parent):
    cdef VariantTreeNode* pivot = root.link[1]
    root.link[1] = pivot.link[0]
    pivot.link[0] = root 

    update_node_height(root)
    update_node_height(pivot)

    if not parent:
        return
    elif parent.key > root.key:
        parent.link[0] = pivot
    else:
        parent.link[1] = pivot

    update_node_height(parent)


cdef void rotate_double_left(VariantTreeNode* root, VariantTreeNode* parent):
    rotate_right(root.link[1], root)
    rotate_left(root, parent)


cdef void rotate_double_right(VariantTreeNode* root, VariantTreeNode* parent):
    rotate_left(root.link[0], root)
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
