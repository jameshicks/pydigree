from math import log
from scipy import stats 

def LikelihoodRatioTest(null_model, alt_model):
	"""
	Compares two nested models by likelihood ratio test

	:returns: Result of test
	:rtype: LikelhoodRatioTestResult
	"""

	chisq = -2.0 * null_model.loglikelihood() + 2.0 * alt_model.loglikelihood()
	df = null_model.df - alt_model.df # Null model has more DFs
	res = LikelihoodRatioTestResult(chisq, df, stats.chi2, len(null_model.observations()))
	return res

class LikelihoodRatioTestResult(object):
	def __init__(self, statistic, df, distribution, n):
		self.statistic = statistic
		self.df = df
		self.distribution = distribution
		self.n = n
	
	@property
	def pvalue(self):
		return self.distribution.sf(self.statistic, self.df)

	@property
	def lod(self):
		return self.statistic / (2.0 * log(10.0))