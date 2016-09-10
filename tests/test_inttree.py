from nose.tools import assert_raises
import random
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
    assert list(b.traverse()) == [25, 26, 27, 28, 29, 30, 31, 32, 33, 34]

def test_ranges():
    tree = IntTree.from_keys([1,3,5,7,9,11,13,15])
    assert list(tree.getrange(5,11).traverse()) == [5,7,9]

    tree.delrange(5,11)
    assert list(tree.traverse()) == [1,3,11,13,15]

def test_intersect():
    t1 = IntTree.from_keys([1, 3, 5, 7, 9])
    t2 = IntTree.from_keys([3, 6, 7])
    t_intersect = t1.intersection(t2)
    assert list(t_intersect.traverse()) == [3,7]


def test_union():
    t1 = IntTree.from_keys([1, 3, 5, 7, 9])
    t2 = IntTree.from_keys([3,6,7])

    t_union = t1.union(t2)
    assert [x for x in t_union.traverse()] == [1,3,5,6,7,9] 

def test_selfbalancing():
    random.seed(100)
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

    delvals = rvals[:]
    random.shuffle(delvals)

    for k in delvals:
        tree.delete(k)
        assert tree.verify()

def test_special():
    # __contains__
    tree = IntTree.from_keys([10, 5, 8, 3, 20])
    assert 20 in tree
    assert 10 in tree
    assert 15 not in tree
    assert 1000 not in tree

    emtree = IntTree()
    assert 1 not in emtree

    # __nonzero__
    assert bool(tree)

    assert not bool(emtree)

    # __len__
    assert len(emtree) == 0
    assert len(tree) == 5