from math import log
from scipy import stats 

def LikelihoodRatioTest(null_model, alt_model):
	chisq = -2.0 * null_model.loglikelihood() + 2.0 * alt_model.loglikelihood()
	df = alt_model.df - null_model.df 
	res = LikelihoodRatioTestResult(chisq)

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