
def maketestdata():
	from pydigree.datastructures import SortedPairContainer
	testdata = [(1,'a'), (5, 'b'), (6,'c'), (9,'d')]
	t =  SortedPairContainer(testdata)
	return t

def test_sparsecontainer_slice():
	t = maketestdata()
	assert t[4:6] == [(5,'b'), (6,'c')]

def test_sparsecontainer_get():
	from nose.tools import assert_raises
	t = maketestdata()
	assert t[5] == 'b'

	def should_raise_keyerror(idx):
		return t[idx]

	assert_raises(KeyError, should_raise_keyerror, 4)
	
def test_sparsecontainer_contains():
	t = maketestdata()
	assert 5 in t
	assert 2 not in t

def test_sparsecontainer_getindex(): 
	t = maketestdata()
	assert t.getindex(5) == 1

def test_sparsecontainer_delitem():
	t = maketestdata()

	assert 5 in t
	del t[5]
	assert 5 not in t