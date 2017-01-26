"Classes for specifying the mating structure in a population"

from pydigree.individual import Individual
import numpy as np


class MatingStructure(object):
    """
    A class representing how a population mates, for when fully random mating 
    may not be desired.
    """

    def __init__(self):
        "Create the structure"
        self.cliques = []

    def next_generation(self, pop, gensize):
        """
        Create the next generation by mating people using the current structure

        :pop: Population to be advanced
        :gensize: Number of children to generate
        :type pop: Population
        :type gensize: int 
        :rtype: list of Individuals
        :returns: progeny
        """
        if not self.cliques:
            self.cliques = self.form_cliques(pop)

        nclique = len(self.cliques)

        # Randomly choose one of the cliques to be the parents for this child
        parents = np.random.randint(0, nclique, gensize)

        # Choose the the sexes each child
        sexes = np.random.randint(0, 2, gensize)

        progeny = [self.cliques[parents[i]].mate(pop=pop, 
                                                 sex=sexes[i],
                                                 label=i)
                   for i in range(gensize)]

        return progeny


class MatingClique(object):

    """
    A class that holds a group of individuals that can form offspring

    :ivar pop: current population
    :ivar males: the males in the clique
    :ivar females: the females in the clique
    :type pop: Population
    """

    def __init__(self, pop, males=None, females=None):
        """
        Create a clique
        
        :param pop: The population to draw the clique individuals are 
            drawn from
        :param males: Males in the clique
        :param females: Females in the clique
        """
        self.pop = pop
        self.males = males if males is not None else []
        self.females = females if females is not None else []

    def children_possible(self):
        """
        Are enough males and females to have a child

        :rtype: bool 
        """
        return len(self.males) > 0 and len(self.females) > 0

    def get_male(self):
        "Chooses a random male from the clique"
        if len(self.males) == 1:
            return self.males[0]
        else:
            return np.random.choice(self.males)

    def get_female(self):
        "Chooses a random female from the clique"
        if len(self.females) == 1:
            return self.females[0]
        else:
            np.random.choice(self.females)

    def mate(self, pop=None, sex=None, label=None):
        """
        Generate offspring from the clique. If more than one father or mother 
        is available, they are randomly selected.

        :param pop: Population for the new individual
        :param sex: Sex of the offspring if specified, otherwise offspring 
            randomly sex is chosen
        :param label: Label for the offspring individual
        :type pop: Population
        :type sex: 0,1

        :returns: offspring individual
        :rtype: Individual
        """
        if not self.children_possible():
            raise ValueError("Children not possible from this clique")

        fa = self.get_male()
        ma = self.get_female()

        if sex is None:
            sex = np.random.randint(0, 2)

        child = Individual(pop, label, fa, ma, sex)
        return child


class RandomMating(MatingStructure):
    """
    A MatingStructure representing purely random mating
    """

    def next_generation(self, pop, gensize):
        """
        Create individuals for the next generation by random mating

        :param pop: Parent population
        :param gensize: Size of next generation
        :type pop: Population
        :type gensize: int
        """

        males = pop.males()
        females = pop.females()

        fathers = np.random.randint(0, len(males), gensize)
        mothers = np.random.randint(0, len(females), gensize)
        sexes = np.random.randint(0, 2, gensize)

        progeny = [Individual(pop, i,
                              males[fathers[i]], females[mothers[i]],
                              sexes[i])
                   for i in range(gensize)]

        return progeny


class MonogamousMating(MatingStructure):
    """
    A mating structure where pairs of individals from the opposite sex 
    pair up and form childen together, and not with anyone else
    """
    def form_cliques(self, pop):
        """
        Creates the monogamous pairs

        :param pop: The population we're working with
        :returns: The cliques
        :rtype: List of MatingCliques
        """
        fathers = pop.males()
        mothers = pop.females()

        np.random.shuffle(fathers)
        np.random.shuffle(mothers)

        return [MatingClique(pop, [pa], [ma]) for pa, ma 
                in zip(fathers, mothers)]
