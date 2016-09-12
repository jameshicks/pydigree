import numpy as np
import os
from pydigree.io.vcf import read_vcf
from pydigree.cydigree.vcfparse import vcf_allele_parser

testdir = os.path.dirname(os.path.abspath(__file__))
TESTDATA_DIR = os.path.join(testdir, 'test_data', 'vcf')

def test_vcf():
    testvcf = os.path.join(TESTDATA_DIR, 'test.vcf')
    
    pop = read_vcf(testvcf)
    assert len(pop.individuals) == 3
    assert len(pop.chromosomes) == 2

    # Test individual with good genotypes
    assert (pop['NA00001'].genotypes[1][0].todense() == np.array([0, 0, 1, 0, 0, 0, 1])).all()
    assert (pop['NA00001'].genotypes[1][1].todense() == np.array([0, 0, 2, 0, 1, 1, 1])).all()
    # assert not pop['NA00001'].genotypes[1][0].missing.all()

    # Test individual with bad genotypes.
    # assert pop['NA00003'].genotypes[1][0].missing.all()


    # Test for FILTER == PASS
    # pop = read_vcf(testvcf)
    # assert pop.chromosomes[1].nmark() == 6

    pop = read_vcf(testvcf, freq_info='AF')
    assert list(pop.chromosomes[0].frequencies) == [0.5]
    diff = (pop.chromosomes[1].frequencies - np.array([0.5, 0.17, 0.333, 0, 0, 0, 0]) )
    assert (diff < 0.001).all()

def test_vcf_allele_parser():
    a = "0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,."
    expected = list(reversed([(2,1), (4,1), (5, 1)])) # stacks are FILO
    observed = vcf_allele_parser(a, 0)
    assert observed.tolist() == expected

