from nose.tools import assert_raises
from pydigree.cydigree.datastructures import IntTree

def test_retrieval():
    tree = IntTree.from_pairs(zip([5,10,15,20], [1,2,3,4]))
    assert tree.get(10) == 2
    assert tree.get(5) == 1
    assert tree.get(100, -1) == -1

    assert tree.find(10) == 2
    assert tree.find(5) == 1

    keys = range(50)
    a = IntTree.from_pairs(zip(keys, keys))

    b = a.getrange(25, 35)
    assert list(b.keys()) == [25, 26, 27, 28, 29, 30, 31, 32, 33, 34]


# def test_traverse():
#     keys = [10, 5, 3, 18, 2]
#     tree = IntTree.from_keys([10, 5, 3, 18, 2])
#     assert [x.key for x in tree.traverse()] == list(sorted(keys))
#     assert [x.key for x in tree.traverse(reverse=True)] == list(sorted(keys, reverse=True))
#     assert [x.key for x in tree.to_stack()] == list(sorted(keys)) 


# def test_minmax():
#     tree = IntTree()

#     rvals = [10, 5, 3, 18, 2, 100, 4123, 4393014, 49310]
#     for i, rval in enumerate(rvals):
#         tree.insert(rval)

#     assert tree.min() == 2
#     assert tree.max() == 4393014

#     emtree = IntTree()

# def test_del():
#     rvals = [48, 23, 74, 3, 44, 64, 98, 41, 56, 91]


#     # Delete leaf
#     tree = IntTree.from_keys([10,5,25,3,8])
#     assert tree.size() == 5
#     assert tree.root.key == 10
#     assert tree.root.left.key, tree.root.right.key == (5,25)

#     # Is tree.root.right a leaf?
#     assert tree.root.right.right is None and tree.root.right.left is None

#     tree.delete(8)
#     assert tree.root.key == 10
#     assert tree.root.left.key, tree.root.right.key == (5,25)
    
#     # Is tree.root.right a leaf?
#     assert tree.root.right.right is None and tree.root.right.left is None
#     assert tree.root.left.left.key == 3
#     assert tree.root.left.right is None
#     assert 8 not in {x.key for x in tree.traverse()}

#     # Delete node with one child
#     tree = IntTree.from_keys(rvals)
#     tree.delete(64)
#     assert tree.root.right.left.key == 56

#     # Symmetric one child case (we need to add another node)
#     tree = IntTree.from_keys(rvals)
#     tree.insert(45)
#     tree.delete(45)
#     assert tree.root.left.right.right is None

#     # Delete two nodes 
#     tree = IntTree.from_keys([10, 5, 25, 3, 8])
#     assert tree.root.key == 10
#     assert tree.root.left.key == 5
#     assert tree.root.right.key == 25
#     assert tree.root.left.left.key == 3

#     tree.delete(5)
#     assert tree.root.key == 10
#     assert tree.root.left.key == 3
#     #
#     # TODO: test this correctly (nose.assert_raises doesnt work here)
#     # try:
#     #     tree.delete(900)
#     #     assert False
#     # except KeyError:
#     #     pass

#     # emtree = IntTree()
#     # try:
#     #     emtree.delete(1)
#     #     assert False
#     # except KeyError:
#     #     pass

#     tree = IntTree.from_keys([100, 50, 200])
#     tree.delete(100)
#     assert list(tree.keys()) == [50,200]

#     tree = IntTree.from_keys([0, 5, 10, 15, 25, 30, 35, 40, 45, 50])
#     tree.delrange(10,40)
#     assert list(tree.keys()) == [0,5,40,45,50]

#     tree = IntTree.from_keys([0, 5, 10, 15, 25, 30, 35, 40, 45, 50])
#     assert tree.size() == 10
#     tree.clear()
#     assert tree.size() == 0

# def test_intersect():
#     t1 = IntTree.from_keys([1, 3, 5, 7, 9])
#     t2 = IntTree.from_keys([3,6,7])
#     t_intersect = t1.intersection(t2)
#     assert [x.key for x in t_intersect.traverse()] == [3,7]

# def test_union():
#     t1 = IntTree.from_keys([1, 3, 5, 7, 9])
#     t2 = IntTree.from_keys([3,6,7])

#     t_union = t1.union(t2)
#     assert [x.key for x in t_union.traverse()] == [1,3,5,6,7,9] 

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

    assert tree.verify() 
    print('Verified after insertion')
    for k in [5409, 3875, 1315, 5418, 1323, 1838, 7103, 6082, 963, 6084]:
        print('TEST DELETING {}'.format(k))
        tree.delete(k)
        assert tree.verify()
        print('VERIFIED')
    print('SELFBALANCING')
# def test_special():
#     # __contains__
#     tree = IntTree.from_keys([10, 5, 8, 3, 20])
#     assert 20 in tree
#     assert 10 in tree
#     assert 15 not in tree
#     assert 1000 not in tree

#     emtree = IntTree()
#     assert 1 not in emtree

#     # # __nonzero__
#     # if not tree:
#     #     assert False

#     # if emtree:
#     #     assert False

#     # __len__
#     assert len(emtree) == 0
#     assert len(tree) == 5