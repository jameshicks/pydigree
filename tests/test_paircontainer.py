
def test_sparsecontainer_slice():
	from pydigree.datastructures import SortedPairContainer
	testdata = [(1,'a'), (5, 'b'), (6,'c'), (9,'d')]
	t =  SortedPairContainer(testdata)

	assert t[4:6] == [(5,'b'), (6,'c')]

def test_sparsecontainer_get():
	from pydigree.datastructures import SortedPairContainer
	from nose.tools import assert_raises
	testdata = [(1,'a'), (5, 'b'), (6,'c'), (9,'d')]
	t = SortedPairContainer(testdata)
	assert t[5] == 'b'

	def should_raise_keyerror(idx):
		return t[idx]

	assert_raises(KeyError, should_raise_keyerror, 4)