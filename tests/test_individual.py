

import os
import glob
from pydigree.io import read_ped

PEDDIR = os.path.abspath(os.path.join(os.path.abspath(__file__), '..', '..', 'sample_pedigrees'))
peds = {}
for filename in glob.glob(PEDDIR + '/*ped'):
    ped = list(read_ped(filename).pedigrees)[0]
    pedname = os.path.basename(filename)[:-4]
    peds[pedname] = ped

def test_parents():
	ped = peds['fullsib']
	ped['1'].parents() == (None, None)
	ped['3'].parents() == (ped['1'], ped['2'])

def test_isfounder():
	ped = peds['first_cousins']
	assert ped['1'].is_founder()
	assert not ped['3'].is_founder()
	assert ped['4'].is_founder()

def test_depth():
	ped = peds['first_cousins']
	# Parental generation
	assert ped['1'].depth == 0
	# Marryin founder
	assert ped['4'].depth == 0
	# F1
	assert ped['3'].depth == 1
	# F2
	assert ped['7'].depth == 2
