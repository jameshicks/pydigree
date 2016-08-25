import numpy as np
import os
from pydigree.io.vcf import read_vcf, vcf_allele_parser

testdir = os.path.dirname(os.path.abspath(__file__))
TESTDATA_DIR = os.path.join(testdir, 'test_data', 'vcf')

def test_vcf():
    testvcf = os.path.join(TESTDATA_DIR, 'test.vcf')
    
    pop = read_vcf(testvcf)
    assert len(pop.individuals) == 3
    assert len(pop.chromosomes) == 2

    # Test individual with good genotypes
    assert (pop['NA00001'].genotypes[1][0].todense() == ['0', '0', '1', '0', '0', '0', '1']).all()
    assert (pop['NA00001'].genotypes[1][1].todense() == ['0', '0', '2', '0', '1', '1', '1']).all()
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

def test_vcf_alleleparser():
    assert vcf_allele_parser('./.') == ('.', '.')
    assert vcf_allele_parser('1/1') == ('1', '1')
    assert vcf_allele_parser('2/1') == ('2', '1')
    assert vcf_allele_parser('1|2') == ('1', '2')
    assert vcf_allele_parser('10/1') == ('10', '1')
    assert vcf_allele_parser('1|10') == ('1', '10')
    assert vcf_allele_parser('10/10') == ('10','10')