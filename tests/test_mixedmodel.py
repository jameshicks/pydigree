import numpy as np
import pydigree as pyd

from pydigree.stats.mixedmodel.mixedmodel import make_incidence_matrix

def test_make_incidence_matrix():
	phenlab = 'testvar'
	inds = [pyd.Individual(None, i) for i in xrange(6)]
	phens = [1,2,3,1,2,3]
	for ind, phen in zip(inds, phens):
		ind.phenotypes[phenlab] = phen

	observed = make_incidence_matrix(inds, phenlab)
	expected = np.array([1,0,0,0,1,0,0,0,1] * 2).reshape(-1,3)
	assert (observed==expected).all()
