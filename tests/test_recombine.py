from pydigree.recombination import recombine
from pydigree.genotypes import Alleles
import numpy as np

def test_recombine():
	a = Alleles(np.zeros(10))
	b = Alleles(np.ones(10))
	m = np.arange(1,100,10)

	n = recombine(a, b, m)
	assert len(n) == len(m)
	assert type(n) == type(a) == type(b)
