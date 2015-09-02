from pydigree.genotypes import ChromosomeTemplate

def test_chromosometemplate():
	# Test the marker finder
	c = ChromosomeTemplate()
	for i in xrange(1,100):
		c.add_genotype(map_position=i, bp=(i*1e6))

	assert c.closest_marker(0) == 0
	assert c.closest_marker(5000001) == 4
	assert c.closest_marker(5999999) == 5
	assert c.closest_marker(1e10) == c.nmark() - 1 