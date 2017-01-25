from math import log
from scipy.stats import chi2

from pydigree.stats.mixedmodel import MixedModel, RandomEffect


class VarianceComponentsLinkage(object):

    def __init__(self, pedigrees, outcome=None, fixed_effects=None,
                 ibd_matrix=None, null_model=None, joint=False, verbose=False,
                 maximization='Average Information'):
        self.pedigrees = pedigrees
        self.null_model = null_model
        self.alternative_model = None
        self.outcome = outcome

        if fixed_effects is None:
            fixed_effects = []
        self.fixed_effects = fixed_effects

        self.analysis_individuals = [x for x in pedigrees.individuals if
                                     self.outcome in x.phenotypes]
        if self.fixed_effects:
            for effect in self.fixed_effects:
                self.analysis_individuals = [x for x in self.analysis_individuals
                                             if effect in x.phenotypes]
        self.ibd_matrix = ibd_matrix
        self.joint = joint
        self.verbose = verbose
        self.maximization = maximization

    def fit_null_model(self):
        null_model = MixedModel(self.pedigrees,
                                outcome=self.outcome,
                                fixed_effects=self.fixed_effects)
        null_model.add_genetic_effect()
        null_model.fit_model()
        null_model.maximize(method=self.maximization, verbose=self.verbose)
        return null_model

    def fit_alternative_model(self):
        ibd_model = MixedModel(self.pedigrees,
                               outcome=self.outcome,
                               fixed_effects=self.fixed_effects)
        ibd_model.add_genetic_effect()

        ranef = RandomEffect(self.analysis_individuals,
                             'IBD',
                             incidence_matrix='eye',
                             covariance_matrix=self.ibd_matrix)
        ibd_model.add_random_effect(ranef)

        ibd_model.fit_model()
        ibd_model.maximize(verbose=self.verbose, method=self.maximization)
        return ibd_model

    def fit(self):
        if not self.joint:
            available_peds = {ind.full_label[0]
                              for ind in self.pedigrees.individuals}
            results = []
            for pedigree_label in available_peds:
                ped = self.pedigrees[pedigree_label]
                mod = VarianceComponentsLinkage(self.pedigrees,
                                                outcome=self.outcome,
                                                fixed_effects=self.fixed_effects,
                                                joint=True,
                                                maximization=self.maximization)
                result = mod.fit()
                results.append(result)
            lodsum = sum(x.lod for x in results)
            obj = VarianceComponentsLinkageResult(lod=lodsum)
            return obj

        if not self.null_model:
            self.null_model = self.fit_null_model()

        if not self.alternative_model:
            self.alternative_model = self.fit_alternative_model()

        result = VarianceComponentsLinkageResult(null_llik=self.null_model.loglikelihood(),
                                                 alt_llik=self.alternative_model.loglikelihood())
        return result


class VarianceComponentsLinkageResult(object):

    def __init__(self, null_llik=None, alt_llik=None, lod=None):
        if lod is not None:
            self.lod = lod

        else:
            self.llik_null = null_llik
            self.llik_alt = alt_llik

            self.chisq = -2.0 * self.llik_null + 2.0 * self.llik_alt
            self.pvalue = 1 - chi2.cdf(self.chisq, 1)
            self.lod = self.chisq / (2.0 * log(10.0))
