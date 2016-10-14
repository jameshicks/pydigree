from pydigree.cydigree.varianttree import VariantTree

def test_basics():
    vt = VariantTree()
    assert vt.refcode == 0
    assert vt.node_count() == 0 
    assert vt.get_item(10) == 0
    assert 10 not in vt
    assert not list(vt.keys())
    assert not list(vt.values())

    vt.set_item(10, 1)
    assert vt.node_count() == 1
    assert vt.get_item(10) == 1
    assert list(vt.keys()) == [10]
    assert list(vt.values()) == [1]
    assert 10 in vt

    vt.set_item(100, 1)
    assert vt.node_count() == 2
    assert vt.get_item(100) == 1
    assert list(vt.keys()) == [10, 100]
    assert list(vt.values()) == [1,1]
    assert 100 in vt

    vt.set_item(11, 2)
    assert vt.node_count() == 2
    assert vt.get_item(11) == 2
    assert list(vt.keys()) == [10,11,100]
    assert list(vt.values()) == [1,2,1]
    assert 11 in vt

    vt.clear_item(11)
    assert vt.node_count() == 2
    assert vt.get_item(11) == 0
    assert list(vt.keys()) == [10,100]
    assert list(vt.values()) == [1,1]
    assert 11 not in vt

    vt.clear_item(100)
    assert vt.node_count() == 1
    assert 100 not in vt

    vt.clear_item(10)
    assert vt.node_count() == 0
    assert 10 not in vt

    import random
    random.seed(1)
    rvals = list(set([random.randint(1,10000) for x in range(1000)]))
    vt = VariantTree()
    assert vt.empty()

    for rval in rvals:
        assert rval not in vt
        vt.set_item(rval, 1)
        assert rval in vt

    assert vt.node_count() == len({x // vt._binsize() for x in rvals})
    random.shuffle(rvals)

    for rval in rvals:
        assert rval in vt
        vt.clear_item(rval)
        assert rval not in vt

    assert vt.node_count() == 0
    assert vt.empty()


def test_selfbalancing():
    size = 200000
    import random
    random.seed(1)
    randvals = [random.randint(1, 10000000) for x in range(size)]

    vt = VariantTree()
    for i,x in enumerate(randvals):
        vt.set_item(x,1)
        if i % 1000 == 0:
            assert vt.verify()

    random.shuffle(randvals)

    for i,x in enumerate(randvals):
        vt.clear_item(x)
        if i % 1000 == 0:
            assert vt.verify()

def test_ranges():
    testvals = [(1,1), (100,1), (25, 5), (50, 5), (75, 5)]

    vt = VariantTree()
    for k,v in testvals: 
        vt.set_item(k,v)

    g = vt.get_range(10,90)
    assert g.keys() == [25, 50, 75]
    assert g.values() == [5,5,5]

    vt.clear_range(20,70)
    assert vt.keys() == [1,75,100]
    assert vt.values() == [1,5,1]
    assert vt.node_count() == 3