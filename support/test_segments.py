import sys

from itertools import izip

import numpy as np

import pydigree as pyd
from pydigree.io import smartopen
from pydigree.sgs.sgs import intervals_to_array
from pydigree.ibs import ibs

replicate = sys.argv[1]
ms = int(sys.argv[2])
prefix='null'

peds = pyd.io.plink.read_plink('{}-{}.ped'.format(prefix, replicate), '{}.map'.format(prefix))
ped = peds['1']
s = pyd.sgs.sgs_population(ped, seed_size=ms)


with smartopen('{}-{}.ibd.gz'.format(prefix, replicate)) as f:
    trueibd = {}
    for line in f:
        fam, id1, id2, ibd_states = line.strip().split(None, 3)
        trueibd[frozenset({id1,id2})] = np.array([int(x) for x in ibd_states.split()])

a = intervals_to_array(s[frozenset({ped['7'],ped['8']})][0], ped.chromosomes[0].nmark())
b = trueibd[frozenset({'7','8'})]


genos1 = izip(*ped['7'].genotypes[0])
genos2 = izip(*ped['8'].genotypes[0])
identical = [ibs(x,y) for x,y in izip(genos1, genos2)]

from pydigree.common import table, runs


for start, stop in runs(list(a), lambda x: x>0, ms):
    print 'Predicted segment: {}-{}'.format(start, stop)

print

for start, stop in runs(list(b), lambda x: x>0, 2):
    print 'True IBD Segment: {}-{}'.format(start, stop)

correct_calls = a == b
print 'Accuracy: {}'.format(correct_calls.sum() / float(correct_calls.shape[0]))
