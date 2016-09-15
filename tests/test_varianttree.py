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