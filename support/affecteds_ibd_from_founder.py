#!/usr/bin/env python
# Marks affected individuals as IBD from a common specified ancestor
# Usage: <script> ped fam ancestor chr pos allele

import sys
import pydigree

peds = pydigree.read_ped(sys.argv[1])
fam, ancestor, chr, pos, allele = sys.argv[2:]
print ' '.join(['GENOTYPE', fam, ancestor, chr, pos, allele, 'P', 'set'])

for ind in peds.individuals:
    if not ind.phenotypes['affected']:
        continue
    print ' '.join(['IBD', fam, ind.label, ancestor, chr, pos, 'P'])
