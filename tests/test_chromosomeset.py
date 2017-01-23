import pydigree

def test_counters():
    nchrom = 50 
    pop = pydigree.Population()
    for i in range(nchrom):
        c = pydigree.ChromosomeTemplate()
        c.add_genotype(0.1, 0)
        pop.add_chromosome(c)

    assert pop.chromosomes.nchrom() == nchrom
    assert pop.chromosomes.nloci() == nchrom

def test_randomloci():
    nchrom = 100
    nloc = 50

    pop = pydigree.Population()
    for i in range(nchrom):
        c = pydigree.ChromosomeTemplate()
        c.add_genotype(0.1, 0)
        pop.add_chromosome(c)

    for i in range(100):
        locs = list(pop.chromosomes.select_random_loci(nloc))
        locs.sort()

        # all the chromosomes have 1 marker
        assert all(x[1] == 0 for x in locs)

        # all the chromosomes are valid indices
        assert all(0 <= x[0] < nchrom for x in locs)

        # No duplicates!
        for i in range(1, nloc):
            assert locs[i] != locs[i-1]
