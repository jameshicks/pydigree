"A linear mixed effects model"

from itertools import product
import copy

import numpy as np
from numpy.linalg import LinAlgError

from scipy.sparse import csc_matrix, issparse
from scipy.sparse import eye as sparseeye
from scipy.linalg import inv as scipy_inv
from scipy.linalg import pinv

from pydigree.stats.mixedmodel.likelihood import makeP
from pydigree.stats.mixedmodel.likelihood import full_loglikelihood
from pydigree.stats.mixedmodel.likelihood import REML, ML

from pydigree.stats.mixedmodel.maximization import newtonlike_maximization
from pydigree.stats.mixedmodel.maximization import expectation_maximization
# from pydigree.stats.mixedmodel.maximization import minque
from pydigree.stats.mixedmodel.maximization import grid_search
from pydigree.stats.mixedmodel.maximization import MLEResult


def is_genetic_effect(effect):
    """
    Is this effect a genetic effect?

    :rtype: bool
    """
    return effect in set(['additive', 'dominance', 'mitochondrial'])


def inv(M):
    "Invert a matrix. If sparse, convert to dense first"
    if issparse(M):
        M = M.todense()
    return scipy_inv(M)


def make_incidence_matrix(individuals, effect_name):
    if effect_name.lower() == 'residual':
        incidence_matrix = sparseeye(len(individuals))

    elif is_genetic_effect(effect_name):
        incidence_matrix = sparseeye(len(individuals))

    else:
        levels = sorted({ind.phenotypes[effect_name] for ind in individuals})

        # Missing values are not a valid level
        levels = [x for x in levels if x is not None]

        nlevels = len(levels)

        # Calculate which individual has which level
        gen = (ind.phenotypes[effect_name] == level for ind, level in
               product(individuals, levels))
        Z = np.fromiter(gen, dtype=np.uint8)

        # Shout out to scipy for both not documenting reshape on any of their
        # sparse matrix objects and also not making them take the same number
        # of arguments
        Z = Z.reshape(-1, nlevels)

        # Check for missing values and complain about them!
        # Kind of hard to read but heres how it goes:
        # Check if any of the rows are all zero.
        if (Z == 0).all(axis=1).any():
            raise LinAlgError('Missing values in random effect')

        incidence_matrix = csc_matrix(Z)

    return incidence_matrix


class RandomEffect(object):
    "A random effect in a mixed model"
    __slots__ = ['label', 
                 'variance_component',
                 'incidence_matrix', 
                 'covariance_matrix', 
                 'levels', 
                 'V_i']

    def __init__(self, individuals, label, variance=None,
                 incidence_matrix=None, covariance_matrix=None, levels=None):
        """
        Create the random effect.

        :param individuals: Individuals included
        :param label: name of the effect
        :param variance: variance associated with the effect
        :param incidence_matrix: incidence matrix for random effect
        :param covariance_matrix: covariance matrix for random effect
        :param levels: levels of random effect
        :type individuals: iterable
        :type label: string
        :type variance: float
        :type incidence_matrix: matrix
        :type covariance_matrix: matrix 
        """
        nobs = len(individuals)

        self.label = label
        self.variance_component = variance
        
        if isinstance(incidence_matrix, str) and incidence_matrix == 'eye':
            self.incidence_matrix = sparseeye(nobs, nobs)
        elif incidence_matrix is None:
            self.incidence_matrix = make_incidence_matrix(individuals,
                                                          self.label)
        else:
            self.incidence_matrix = incidence_matrix

        if covariance_matrix is None:
            # Number of levels of random effects is the number of
            # columns in the incidence matrix
            nlevel = self.incidence_matrix.shape[1]
            self.covariance_matrix = sparseeye(nlevel, nlevel)
        else:
            # Covariance matrices are square
            if covariance_matrix.shape[0] != covariance_matrix.shape[1]:
                raise LinAlgError('Covariance matrix not square')
            if covariance_matrix.shape[0] != self.incidence_matrix.shape[1]:
                raise LinAlgError('Incidence and covariance matrix '
                                  'not conformable')
            self.covariance_matrix = covariance_matrix

        if not levels:
            self.levels = ['L{}'.format(i) for i in 
                            range(self.incidence_matrix.shape[1])]
        else:
            if len(levels) != incidence_matrix.shape[1]:
                raise ValueError('Number of levels not correct')
            self.levels = levels

        self.V_i = self.Z * self.G * self.Z.T

    def __repr__(self):
        return 'Random Effect: {}'.format(self.label)

    @property
    def nlevels(self):
        """
        The number of levels of the random effect

        :rtype: int
        """
        return len(self.levels)

    # Convenience properties for linear algebra
    @property
    def sigma(self):
        """ 
        Convenience property for returning the variance of the component 
        
        :rtype: float
        """
        return self.variance_component

    @property
    def Z(self):
        "Convenience property for returning the incidence matrix"
        return self.incidence_matrix

    @property
    def G(self):
        "Convenience property for returning the covariance_matrix"
        return self.covariance_matrix

class MixedModel(object):

    """
    Fits linear models in the form of y = X * b + sum(Z_i * u_i) + e, where:
      y is the vector of outcomes
      X is a design matrix of fixed effects
      b is the vector of coefficients correstponding to those fixed effects
      Z_i is an incidence matrix corresponding to random effect i
      u_i is a vector of values corresponding to random effect i
      e is a vector of errors
    """

    def __init__(self, pedigrees, outcome=None, fixed_effects=None,
                 random_effects=None, covariance_matrices=None, only=None):
        self.mle = None
        self.pedigrees = pedigrees

        self.outcome = outcome
        self.fixed_effects = fixed_effects if fixed_effects else []
        self.obs = []

        if only is not None:
            self.only = frozenset(only)

        if not random_effects:
            self.random_effects = []
        else:
            if not all(isinstance(x, RandomEffect) for x in random_effects):
                raise ValueError(
                    'Random effects must be of class RandomEffect')
            self.random_effects = random_effects

        # Extra variance component for residual variance. Works like
        # any other random effect.
        residual = RandomEffect(self.observations(), 'Residual')
        self.random_effects.append(residual)

        self.V = None
        self.X = None
        self.y = None
        self.Zlist = None
        self.beta = None

    def copy(self):
        """ 
        Returns a deep copy of the model 

        :returns: A copy of the model object
        :rtype: MixedModel
        """
        # We want to avoid copying pedigree and individual data, so
        # we'll set the pedigrees attribute to None for a sec, and then
        # change it back
        peds = self.pedigrees
        self.pedigrees = None

        newmm = copy.deepcopy(self)

        newmm.pedigrees = peds
        self.pedigrees = peds

        return newmm

    def fit_model(self):
        """
        Builds X, Z, Y, and R for the model

        
        :returns: void
        """
        self.y = self._makey()
        self.X = self._makeX()
        self.Zlist = self._makeZs()

        need_vcs = not all(x is not None for x in self.variance_components)
        if need_vcs:
            # Not a great way to start but you need to start with something
            need_vcs = True
            vcs = [0] * len(self.random_effects)
            vcs[-1] = np.var(self.y)
            self.set_variance_components(vcs)

        self.V = self._makeV()
        self.beta = self._makebeta()

        if need_vcs:
            vcs[-1] = np.var(self.y - self.X * self.beta)

    def _fit_results(self):
        self.V = self._makeV()
        self.beta = self._makebeta()

    def clear_model(self):
        """ Clears all parameters from the model """
        self.random_effects = []
        self.fixed_effects = []
        self.Zlist = []
        self.X = None
        self.y = None
        self.beta = None
        self.V = None
        self.obs = None

    def observations(self):
        """
        Returns a list of the fully observed individuals in the model
        Fully observed individuals have observations for each fixed effect and
        and observation for the outcome variable.

        :returns: the fully observed individuals
        :rtype: list of Individuals
        """

        def has_all_fixefs(ind, effects):
            if not effects:
                return True
            for effect in effects:
                if effect not in ind.phenotypes:
                    return False
                elif ind.phenotypes[effect] is None:
                    return False
            return True

        def has_outcome(ind):
            try:
                return ind.phenotypes[self.outcome] is not None
            except KeyError:
                return False

        obs = [x for x in self.pedigrees.individuals
               if (has_all_fixefs(x, self.fixed_effects) and 
                   has_outcome(x) and x)]
        
        if not self.obs:
            self.obs = obs
        
        return obs

    def nobs(self):
        """
        :returns: the number of fully observed individuals in the model
        :rtype: integer
        """
        return len(self.observations())

    @property
    def variance_components(self):
        """
        The current variances associated with each random effect

        :rtype: list of floats
        """
        return [x.sigma for x in self.random_effects]

    @property
    def covariance_matrices(self):
        """
        The covariance matrices associated with each random effect

        :rtype: list of matrices
        """
        return [x.covariance_matrix for x in self.random_effects]

    @property
    def R(self):
        "Covariance matrix of the residual variance"
        return self.random_effects[-1].covariance_matrix

    @property
    def P(self):
        "Projection matrix"
        return makeP(self.X, inv(self.V))

    def residual_variance(self):
        """ 
        Returns the variance in y not accounted for by random effects 
        
        :rtype: float
        """
        return self.random_effects[-1].sigma


    def _makey(self):
        """ Prepares the vector of outcome variables for model estimation """
        obs = self.observations()
        return np.matrix([x.phenotypes[self.outcome] for x in obs]).transpose()

    def _makeX(self):
        """
        Builds the design matrix for the fixed effects in the mixed model.
        Includes a column of ones for the intercept.

        Returns: matrix
        """
        obs = self.observations()
        xmat = [[1] * len(obs)]
        for phen in self.fixed_effects:
            xmat.append([ob.phenotypes[phen] for ob in obs])
        X = np.matrix(list(zip(*xmat)))
        return X

    def _makeZs(self):
        """
        Makes the incidence matrix for random effects

        :rtype: A list of numpy matrices
        """
        Zlist = [ranef.Z for ranef in self.random_effects]

        return [csc_matrix(Z) for Z in Zlist]

    def _makeV(self, vcs=None):
        if vcs is None and (not all(x is not None for x in self.variance_components)):
            raise ValueError('Variance components not set')
        
        if vcs is None:
            variance_components = self.variance_components
        else:
            variance_components = vcs

        V = sum(sigma * Z * A * Z.T for sigma, Z, A in
                zip(variance_components,
                     self.Zlist,
                     self.covariance_matrices))

        return V

    def _makebeta(self):
        """
        Calculates BLUEs for the fixed effects portion of the model

        Reference:
        McCulloch & Seale. Generalized, Linear, and Mixed Models. (2001)
        Equation 6.24
        """
        vinv = inv(self.V.todense())
        return pinv(self.X.T * vinv * self.X) * self.X.T * vinv * self.y

    def set_outcome(self, outcome):
        """ Sets the outcome for the mixed model """
        self.outcome = outcome
        self.y = self._makey()

    def add_fixed_effects(self, effect):
        """ Adds a fixed effect to the model """
        self.fixed_effects.append(effect)
        self.X = self._makeX()

    def add_random_effect(self, effect):
        """ Adds a random effect to the model """
        if not isinstance(effect, RandomEffect):
            raise ValueError('Random effect must be of type RandomEffect')
        self.random_effects.insert(-1, effect)
        self.Zlist = self._makeZs()

    def add_genetic_effect(self, kind='additive'):
        """
        Adds a genetic effect to the model as a random effect

        :param kind: type of effect to add
        :type kind: 'additive' or 'dominance'
        """
        inds = [x.full_label for x in self.observations()]
        peds = self.pedigrees

        if kind.lower() == 'additive':
            covmat = peds.additive_relationship_matrix(inds)
        elif kind.lower() == 'dominance':
            covmat = peds.dominance_relationship_matrix(inds)
        else:
            raise NotImplementedError(
                'Nonadditive/dominance genetic effects not implemented')
        effect = RandomEffect(
            self.observations(), kind, covariance_matrix=covmat)
        self.add_random_effect(effect)

    def set_variance_components(self, variance_components):
        """
        Manually set variance components for each random effect in the model.
        Useful if you know a priori, say a heritability, and just want to
        predict breeding values for the trait.

        :param variance_components: variances associated with each random effect
        :type variance_components: iterable of numerics
        """
        if not all(x is not None for x in variance_components):
            raise ValueError('Not all variance components are specified')
        for sigma, ranef in zip(variance_components, self.random_effects):
            ranef.variance_component = sigma

    def maximize(self, method="Average Information", restricted=False,
                 starts=None, verbose=False):
        """
        Finds the optimal values for variance components in the model using
        provided optimization methods.

        :param restricted: Uses REML estimation
        :param starts: starting values for the variance components
        :param method: maximization method
        :param verbose: output maximization progress
        :type restricted: bool
        :type method: string
        :type starts: iterable of numerics
        :type verbose: bool:  
        """

        if (isinstance(self.mle, MLEResult) and
                self.maximized.method == method):
            return
        self.fit_model()

        if starts is None:
            starts = self._starting_variance_components()

        likefunc = REML if restricted else ML
        llik = likefunc(self, info=method)
        llik.set_parameters(starts)

        # if method.lower().startswith('minque'):
        #     mle = minque(self, value=0, verbose=verbose, starts=starts)

        if method.lower() in {'em', 'emreml', 'expectation-maximization'}:
            mle = expectation_maximization(self, llik, verbose=verbose)

        elif method.lower() == 'grid':
            mle = grid_search(self, llik, nevals=20, oob=False)

        else:
            mle = newtonlike_maximization(self, llik, verbose=verbose)

        self.mle = mle
        self.set_variance_components(mle.parameters)
        self.fit_model()

        # Get the full loglikelihood at the REML maximimum so we
        # can use it later
        self.mle.full_loglikelihood = full_loglikelihood(self.y, self.V,
                                                         self.X, self.beta)

    @property
    def maximized(self):
        """
        Has the model been maximized? 

        :rtype: bool
        """
        return isinstance(self.mle, MLEResult)

    def loglikelihood(self, restricted=False, vcs=None, vmat=None):
        """
        Returns the loglikelihood of the model with the current model 
        parameters

        :returns: loglikelihood
        :rtype: float
        """
        if self.mle is not None and vcs is None and vmat is None:
            if restricted:
                return self.mle.restricted_loglikelihood
            else:
                return self.mle.full_loglikelihood

        if vcs is not None:
            V = self._makeV(vcs=vcs)
        elif vmat is not None:
            V = vmat
        elif self.V is None:
            self.V = self._makeV()
        else:
            V = self.V
        if not restricted:
            return full_loglikelihood(self.y, V, self.X, self.beta)
        else:
            return REML(self).loglikelihood()

    @property
    def df(self):
        '''
        The number of observations minus the number of fixed effects, minus
        the number of non-residual random effects

        :rtype: integer
        '''
        return self.nobs() - self.X.shape[1] - len(self.random_effects) + 1

    @property
    def bic(self):
        """ 
        Calculates the Bayesian Information Criterion (BIC) for the model

        :rtype: float
        """

        if not self.maximized:
            raise ValueError('Model not maximized!')
        # Add 1 because the intercept has to be estimated
        nparam = len(self.fixed_effects) + len(self.random_effects) + 1
        n = self.nobs()
        loglike = self.loglikelihood()

        return -2 * loglike + nparam * np.log(n)

    def blup(self, idx):
        """
        Get the BLUPs for a random effect

        :param idx: index of effect
        :type idx: int

        :rtype: np.array
        """
        rf = self.random_effects[idx]
        res = (self.y - self.X * self.beta) 
        blups = rf.G * rf.Z.T * inv(self.V.todense()) * res
        return np.array(blups.T)[0]

    def summary(self):
        """ 
        Prints a summary of the current model 

        :rtype: void
        """
        if not self.maximized:
            raise ValueError('Model not maximized!')

        self._fit_results()

        print()
        print('Linear mixed model fit by {}'.format(self.mle.method))
        print()
        print('Fixed effects:')
        fixefnames = ['(Intercept)'] + self.fixed_effects
        betas = self.beta.T.tolist()[0]
        print('\t'.join(['Name', 'Estimate']))
        for name, beta in zip(fixefnames, betas):
            print('\t'.join(q for q in [name, '{:5.3f}'.format(beta)]))
        print()
        print('Variance components:')
        print('\t'.join(['Component', 'Variance', '% Variance']))
        totalvar = sum(self.variance_components)
        for effect, vc in zip(self.random_effects, self.variance_components):
            print('\t'.join(v for v in [effect.label,
                                        '{:5.3f}'.format(vc),
                                        '{:5.3f}'.format(100 * vc / totalvar)]))
        print()
        print('Observations: {}'.format(self.nobs()))
        print('Loglikelihood: {:10.2f}'.format(self.loglikelihood()))
        print('BIC: {:10.3f}'.format(self.bic))
        print()


    def _variance_after_fixefs(self):
        return np.var(self.y - self.X * self.beta)

    def _starting_variance_components(self, kind='equal'):
        """
        Starting variance components in optimization.
        Valid values:

        'ols': Starting values are all 0 except residual, which is 
            var(y - X*Beta)
        'EM': the starting values are the variance components after
            100 iterations of expectation-maximization REML (started from all
            equal values).
        'equal': Chooses all variance components (including residual)
            to be equal.


        :param kind: the method to find starting values
        :type kind: string

        :returns: variance components
        :rtype: numpy array of floats
        """
        # 'minque0': Starting values are those from MINQUE with all weights
        #     set equal to 0 except for the residual variance, which is set
        #     to 1. This is the default method used by SAS's PROC MIXED.
        # 'minque1': Starting values are those from MINQUE with all weights
        #     set equal to 1
        # if kind.lower() == 'minque0':
        #     return minque(self, value=0, return_after=1, return_vcs=True)

        # if kind.lower() == 'minque1':
        #     return minque(self, value=1, return_after=1, return_vcs=True)

        # if kind.lower() == 'minquemean':
        #     zero = minque(self, value=0, return_after=1, return_vcs=True)
        #     one = minque(self, value=1, return_after=1, return_vcs=True)
        #     return (zero + one) / 2.0

        if kind.lower() == 'ols':
            vcs_start = np.zeros(len(self.random_effects))
            vcs_start[-1] = self._variance_after_fixefs()
            return vcs_start

        if kind.lower() == 'equal':
            v = self._variance_after_fixefs()
            n = len(self.random_effects)
            vcs_start = [v/float(n)] * n
            return vcs_start

        if kind.lower() == 'em':
            starts = self._starting_variance_components('equal')
            vcs_start = expectation_maximization(self,
                                                 REML(self),
                                                 starts=starts,
                                                 return_after=100)
            return vcs_start.parameters
        else:
            raise ValueError('Unknown method: {}'.format(kind))
