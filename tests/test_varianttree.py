from pydigree.cydigree.varianttree import VariantTree

def test_basics():
	vt = VariantTree()
	assert vt.refcode == 0
	assert vt.node_count() == 0 
	assert vt.get_item(10) == 0
	assert not list(vt.keys())
	assert not list(vt.values())

	vt.set_item(10, 1)
	assert vt.node_count() == 1
	assert vt.get_item(10) == 1
	assert list(vt.keys()) == [10]
	assert list(vt.values()) == [1]

	vt.set_item(100, 1)
	assert vt.node_count() == 2
	assert vt.get_item(100) == 1
	assert list(vt.keys()) == [10, 100]
	assert list(vt.values()) == [1,1]

	vt.set_item(11, 2)
	assert vt.node_count() == 2
	assert vt.get_item(11) == 2
	assert list(vt.keys()) == [10,11,100]
	assert list(vt.values()) == [1,2,1]
	
	vt.clear_item(11)
	assert vt.node_count() == 2
	assert vt.get_item(11) == 0
	assert list(vt.keys()) == [10,100]
	assert list(vt.values()) == [1,1]

	vt.clear_item(100)
	assert vt.node_count() == 1

	vt.clear_item(10)
	assert vt.node_count() == 0

def test_selfbalancing():
	size = 200000
	import random
	random.seed(1)
	randvals = [random.randint(1, 10000000) for x in range(size)]

	vt = VariantTree()
	for i,x in enumerate(randvals):
		vt.set_item(x,1)
		if i % 100 == 0:
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

	g = vt.getrange(10,90)
	assert g.keys() == [25, 50, 75]
	assert g.values() == [5,5,5]
