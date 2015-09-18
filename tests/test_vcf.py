import os
from pydigree.io.vcf import read_vcf, vcf_allele_parser

testdir = os.path.dirname(os.path.abspath(__file__))
TESTDATA_DIR = os.path.join(testdir, 'test_data', 'vcf')

def test_vcf():
    testvcf = os.path.join(TESTDATA_DIR, 'test.vcf')
    
    pop = read_vcf(testvcf, minqual=0)
    assert len(pop.individuals) == 3
    assert len(pop.chromosomes) == 2

    # Test individual with good genotypes
    assert (pop['NA00001'].genotypes[1][0].todense() == ['0', '0', '1', '0', '0', '0', '1']).all()
    assert (pop['NA00001'].genotypes[1][1].todense() == ['0', '0', '2', '0', '1', '1', '1']).all()
    assert not pop['NA00001'].genotypes[1][0].missing.all()

    # Test individual with bad genotypes.
    assert pop['NA00003'].genotypes[1][0].missing.all()

    # Test variant level quality filtering 
    pop = read_vcf(testvcf, minqual=20)
    assert pop.chromosomes[1].nmark() == 6

    # Test for FILTER == PASS
    pop = read_vcf(testvcf, minqual=0, require_pass=True)
    assert pop.chromosomes[1].nmark() == 6

def test_vcf_alleleparser():
    assert vcf_allele_parser('./.') == ('.', '.')
    assert vcf_allele_parser('1/1') == ('1', '1')
    assert vcf_allele_parser('2/1') == ('2', '1')
    assert vcf_allele_parser('1|2') == ('1', '2')
    assert vcf_allele_parser('10/1') == ('10', '1')
    assert vcf_allele_parser('1|10') == ('1', '10')
    assert vcf_allele_parser('10/10') == ('10','10')