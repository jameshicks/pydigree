from itertools import izip

import numpy as np

import pydigree as pyd
from pydigree.io import smartopen
from pydigree.sgs.sgs import intervals_to_array

ms=1000

peds = pyd.io.plink.read_plink('null-1.ped','null.map')
ped = peds['1']
s = pyd.sgs.sgs_population(ped, min_seg=ms)


with smartopen('null-1.ibd.gz') as f:
    trueibd = {}
    for line in f:
        fam, id1, id2, ibd_states = line.strip().split(None, 3)
        trueibd[frozenset({id1,id2})] = np.array([int(x) for x in ibd_states.split()])

a = intervals_to_array(s[frozenset({ped['3'],ped['4']})][0], ped.chromosomes[0].nmark())
b = trueibd[frozenset({'3','4'})]

from pydigree.common import table, runs
print table(a)
print table(b)
print

for start, stop in runs(a, lambda x: x>0, ms):
    print 'Predicted segment: {}-{}'.format(start, stop)

print

for start, stop in runs(b, lambda x: x>0, ms):
    print 'True IBD Segment: {}-{}'.format(start, stop)

#for i,data  in enumerate(izip(a,b)):
#    pred, actual = data
#    print i, pred, actual
