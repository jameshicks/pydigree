from math import log
from scipy import stats 

def LikelihoodRatioTest(null_model, alt_model):
	chisq = -2.0 * null_model.loglikelihood() + 2.0 * alt_model.loglikelihood()
	df = null_model.df - alt_model.df # Null model has more DFs
	res = LikelihoodRatioTestResult(chisq, df, stats.chi2)
	return res

class LikelihoodRatioTestResult(object):
	def __init__(self, statistic, df, distribution):
		self.statistic = statistic
		self.df = df
		self.distribution = distribution
	
	@property
	def pvalue(self):
		return self.distribution.sf(self.statistic, self.df)

	@property
	def lod(self):
		return self.statistic / (2.0 * log(10.0))