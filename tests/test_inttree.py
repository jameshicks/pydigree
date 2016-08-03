from nose.tools import assert_raises
from pydigree.cydigree.datastructures import IntTreeNode, IntTree, NodeStack
from pydigree.cydigree.datastructures import rotate_right, rotate_left
from pydigree.cydigree.datastructures import rotate_double_left, rotate_double_right

def test_NodeStack():
    s = NodeStack([IntTreeNode(x) for x in (1,2,3,4)])
    assert bool(s)
    assert not bool(NodeStack())
    assert list(x.key for x in s) == [4,3,2,1]

    s = NodeStack([IntTreeNode(x) for x in (1,2,3,4)])
    s.reverse()
    assert list(x.key for x in s) == [1,2,3,4]

def test_right_rotation():
    a = IntTreeNode(1)
    b = IntTreeNode(2)
    c = IntTreeNode(3)

    c.left = b
    b.left = a
    b.parent = c
    a.parent = b

    # BEFORE:
    #      c
    #     /
    #    b
    #   /
    #  a
    rotate_right(c)

    # AFTER:
    #     b
    #    / \
    #   a   c
    assert a.parent is b
    assert c.parent is b
    assert b.parent is None
    assert b.left is a
    assert b.right is c
    assert a.is_leaf()
    assert c.is_leaf() 

    ######
    a = IntTreeNode(1)
    b = IntTreeNode(2)
    c = IntTreeNode(3)
    d = IntTreeNode(4)

    d.left = c
    c.left = b
    b.left = a
    c.parent = d
    b.parent = c
    a.parent = b

    # BEFORE:
    #        d
    #       /
    #      c
    #     /
    #    b
    #   /
    #  a

    rotate_right(d)

    # AFTER:
    #      c
    #     / \
    #    b   d
    #   /
    #  a 
    assert c.children == (b,d)
    assert b.left is a
    assert b.right is None
    assert a.is_leaf()
    assert d.parent is c and b.parent is c
    assert c.parent is None
    ####
    # Before
    #
    #     5
    #    / \
    #   3   7
    #  / \
    # 2   4 
    two = IntTreeNode(2)
    three = IntTreeNode(3)
    four = IntTreeNode(4)
    five = IntTreeNode(5)
    seven = IntTreeNode(7)

    two.parent, three.left = three, two
    four.parent, three.right = three, four
    three.parent, five.left = five, three
    seven.parent, five.right = five, seven

    rotate_right(five)
    # After
    #    3
    #   / \
    #  2   5
    #     / \
    #    4   7
    assert five.parent is three
    assert two.parent is three
    assert three.parent is None
    assert four.parent is five
    assert seven.parent is five
    assert three.children == (two, five)
    assert five.children == (four, seven)
    assert two.is_leaf() and four.is_leaf() and seven.is_leaf()

def test_left_rotation():
    a = IntTreeNode(1)
    b = IntTreeNode(2)
    c = IntTreeNode(3)
    d = IntTreeNode(4)

    a.right = b
    b.right = c
    c.parent = b
    b.parent = a


    # BEFORE:
    # a
    #  \
    #   b 
    #    \
    #     c
    rotate_left(a)

    # AFTER:
    #     b
    #    / \
    #   a   c
    assert a.parent is b
    assert c.parent is b
    assert b.parent is None
    assert b.left is a
    assert b.right is c
    assert a.is_leaf()
    assert c.is_leaf() 

    ####
    # Before
    #    3
    #   / \
    #  2   5
    #     / \
    #    4   7 
    two = IntTreeNode(2)
    three = IntTreeNode(3)
    four = IntTreeNode(4)
    five = IntTreeNode(5)
    seven = IntTreeNode(7)

    three.children = two, five
    five.children = four, seven
    four.parent, seven.parent = five, five
    two.parent, five.parent = three, three

    rotate_left(three)
    # After
    #
    #     5
    #    / \
    #   3   7
    #  / \
    # 2   4
    assert five.children == (three, seven)
    assert three.parent is five and seven.parent is five 
    assert three.children == (two, four)
    assert two.parent is three and four.parent is three
    assert five.parent is None
    assert two.is_leaf() and four.is_leaf() and seven.is_leaf()

def test_rotate_doubleleft():
    #BEFORE
    #  a
    #   \
    #    c
    #   /
    #  b
    a = IntTreeNode(1)
    b = IntTreeNode(2)
    c = IntTreeNode(3)
    z = IntTreeNode(0)

    a.parent = z
    z.right = a
    a.right = c
    c.left = b
    b.parent = c
    c.parent = a

    rotate_double_left(a)
    # AFTER
    #   b
    #  / \
    # a   c
    #

    assert a.parent is c.parent is b
    assert b.children == (a,c)
    assert b.parent is z
    assert z.right is b


def test_rotate_doubleright():
    # BEFORE
    #   c
    #  /
    # a
    #  \
    #   b

    a = IntTreeNode(1)
    b = IntTreeNode(2)
    c = IntTreeNode(3)

    a.parent = c
    b.parent = a
    c.left = a
    a.right = b

    rotate_double_right(c)
    # AFTER
    #   b
    #  / \
    # a   c
    #

    assert a.parent is c.parent is b
    assert b.children == (a,c)
    assert b.parent is None

def test_insertion():
    tree = IntTree()
    tree.insert(10)
    assert tree.root.key == 10

    tree.insert(5)
    assert tree.root.left.key == 5

    tree.insert(8)

    assert tree.root.key == 8 
    assert tree.root.left.key == 5
    assert tree.root.right.key == 10

    tree.insert(7)
    assert tree.root.left.right.key == 7

    tree.insert(6)
    assert tree.root.left.key == 6
    assert tree.root.left.left.key == 5
    assert tree.root.left.right.key == 7

def test_path():
    tree = IntTree.from_keys([69, 60, 22, 91, 19, 71, 96, 27, 84, 43])
    assert tree.root.key == 60

    expected = [84,91,71,60]
    observed = [x.key for x in tree.path_to_root(84)]
    assert observed == expected

    expected = expected[::-1]
    observed = [x.key for x in tree.path_to_node(84)]
    assert observed == expected

def test_traverse():
    keys = [10, 5, 3, 18, 2]
    tree = IntTree.from_keys([10, 5, 3, 18, 2])
    assert [x.key for x in tree.traverse()] == list(sorted(keys))
    assert [x.key for x in tree.traverse(reverse=True)] == list(sorted(keys, reverse=True))
    assert [x.key for x in tree.to_stack()] == list(sorted(keys)) 


def test_minmax():
    tree = IntTree()

    rvals = [10, 5, 3, 18, 2, 100, 4123, 4393014, 49310]
    for i, rval in enumerate(rvals):
        tree.insert(rval)

    assert tree.min() == 2
    assert tree.max() == 4393014

    emtree = IntTree()

def test_del():
    rvals = [48, 23, 74, 3, 44, 64, 98, 41, 56, 91]


    # Delete leaf
    tree = IntTree.from_keys([10,5,25,3,8])
    assert tree.size() == 5
    assert tree.root.key == 10
    assert tree.root.left.key, tree.root.right.key == (5,25)
    assert tree.root.right.is_leaf()

    tree.delete(8)
    assert tree.root.key == 10
    assert tree.root.left.key, tree.root.right.key == (5,25)
    assert tree.root.right.is_leaf()
    assert tree.root.left.left.key == 3
    assert tree.root.left.right is None
    assert 8 not in {x.key for x in tree.traverse()}

    # Delete node with one child
    tree = IntTree.from_keys(rvals)
    tree.delete(64)
    assert tree.root.right.left.key == 56

    # Symmetric one child case (we need to add another node)
    tree = IntTree.from_keys(rvals)
    tree.insert(45)
    tree.delete(45)
    assert tree.root.left.right.right is None

    # Delete two nodes 
    tree = IntTree.from_keys([10, 5, 25, 3, 8])
    assert tree.root.key == 10
    assert tree.root.left.key == 5
    assert tree.root.right.key == 25
    assert tree.root.left.left.key == 3

    tree.delete(5)
    assert tree.root.key == 10
    assert tree.root.left.key == 3
    #
    # TODO: test this correctly (nose.assert_raises doesnt work here)
    # try:
    #     tree.delete(900)
    #     assert False
    # except KeyError:
    #     pass

    # emtree = IntTree()
    # try:
    #     emtree.delete(1)
    #     assert False
    # except KeyError:
    #     pass

def test_intersect():
    t1 = IntTree.from_keys([1, 3, 5, 7, 9])
    t2 = IntTree.from_keys([3,6,7])
    # assert 0
    t_intersect = t1.intersection(t2)
    # assert 0
    assert [x.key for x in t_intersect.traverse()] == [3,7]

def test_union():
    t1 = IntTree.from_keys([1, 3, 5, 7, 9])
    t2 = IntTree.from_keys([3,6,7])

    t_union = t1.union(t2)
    assert [x.key for x in t_union.traverse()] == [1,3,5,6,7,9] 

def test_selfbalancing():
    tree = IntTree()
    rvals = [7680, 1027, 2564, 4103, 6671, 4118, 5143, 6680, 5144, 5146, 6682, 
             1560, 31, 2079, 546, 4643, 3625, 6188, 4656, 9776, 7732, 1591, 
             8760, 65, 6722, 68, 71, 8265, 8778, 1103, 81, 2642, 9299, 4692, 
             7251, 9809, 2136, 5720, 4700, 7260, 613, 6762, 4203, 5739, 3694, 
             5238, 4727, 6262, 2680, 5240, 4219, 2173, 1661, 8838, 4744, 1161, 
             9866, 4234, 1676, 6799, 9875, 3740, 7837, 5629, 9892, 1700, 9386, 
             1197, 2734, 7855, 6651, 8384, 2242, 1222, 200, 4816, 7382, 3803, 
             4317, 6366, 7389, 1251, 4839, 8938, 7403, 748, 751, 4848, 1268, 
             5878, 247, 4343, 2298, 9978, 6908, 5886, 256, 3845, 1798, 263, 
             8968, 6921, 6410, 3849, 5393, 2325, 6421, 4374, 3352, 6935, 282, 
             6939, 1815, 8987, 5409, 3875, 1315, 5418, 1323, 1838, 5426, 7988, 
             9015, 7481, 826, 4924, 6973, 1859, 6468, 3911, 4424, 1871, 1879, 
             3415, 7513, 3932, 5469, 8546, 1890, 9061, 7018, 7531, 2926, 4462, 
             9587, 7544, 893, 383, 5505, 8065, 7058, 7571, 2967, 1944, 6553, 
             9625, 2968, 5016, 1431, 6559, 5535, 6048, 6562, 1955, 9122, 4518, 
             2982, 3500, 3508, 6072, 7103, 6082, 963, 6084, 3523, 455, 7627, 
             8654, 979, 7637, 6101, 8150, 8664, 7642, 3035, 2525, 480, 2020, 
             5092, 4071, 1512, 3567, 2032, 6132, 4086, 7161, 4091, 1533]

    for i, rval in enumerate(rvals):
        tree.insert(rval)


    for node in tree.traverse():
        assert node.is_balanced() 

    assert tree.find_node(31).key == 31
    assert tree.find_node(6973).key == 6973
    assert tree.find_node(9015).key == 9015
    assert tree.size() == len(rvals)
    # print(tree.size())

    for k in [5409, 3875, 1315, 5418, 1323, 1838, 7103, 6082, 963, 6084]:
        tree.delete(k)
        for x in tree.traverse():
            assert x.is_balanced()

def test_special():
    # __contains__
    tree = IntTree.from_keys([10, 5, 8, 3, 20])
    assert 20 in tree
    assert 10 in tree
    assert 15 not in tree
    assert 1000 not in tree

    emtree = IntTree()
    assert 1 not in emtree

    # # __nonzero__
    # if not tree:
    #     assert False

    # if emtree:
    #     assert False

    # __len__
    assert len(emtree) == 0
    assert len(tree) == 5