"Methods for statistical testing"

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
    """
    The result of a likelihood ratio test.

    :ivar statistic: test statistic
    :ivar df: degrees of freedom
    :ivar distribution: distribution of test statistic
    :ivar n: sample size

    :type statistic: float
    :type df: int
    :type distribution: scipy probability distribution
    :type n: int
    """

    def __init__(self, statistic, df, distribution, n):
        self.statistic = statistic
        self.df = df
        self.distribution = distribution
        self.n = n
    
    @property
    def pvalue(self):
        """
        P-value for the test

        :rtype: float
        """
        return self.distribution.sf(self.statistic, self.df)

    @property
    def lod(self):
        """
        LOD score (log10 LR) of the result

        :rtype float:
        """
        return self.statistic / (2.0 * log(10.0))
