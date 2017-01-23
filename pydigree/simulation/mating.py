from pydigree.individual import Individual
import numpy as np


class MatingStructure(object):

    def __init__(self):
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

        progeny = [self.cliques[parents[i]].mate(pop=pop, sex=sexes[i], label=i)
                   for i in range(gensize)]

        return progeny


class MatingClique(object):

    """
    A class that holds a group of individuals that can form offspring
    """
    males = []
    females = []

    def __init__(self, pop, males, females):
        self.pop = pop
        self.males = males
        self.females = females

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
        :param sex: Sex of the offspring if specified, otherwise randomly chosen
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


class RandomMating(object):

    def next_generation(self, pop, gensize):
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

    def form_cliques(self, pop):
        """
        Creates the monogamous pairs
        """
        fathers = pop.males()
        mothers = pop.females()

        np.random.shuffle(fathers)
        np.random.shuffle(mothers)

        return [MatingClique(pop, [pa], [ma]) for pa, ma in zip(fathers, mothers)]
