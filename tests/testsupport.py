import os
import glob
from pydigree.io import read_ped


def getpeds():
	''' Read the pedigrees from sample_pedigrees and return a dict '''
	PEDDIR = os.path.abspath(os.path.join(os.path.abspath(__file__), '..', '..', 'sample_pedigrees'))
	peds = {}
	for filename in glob.glob(PEDDIR + '/*ped'):
		ped = list(read_ped(filename).pedigrees)[0]
		pedname = os.path.basename(filename)[:-4]
		peds[pedname] = ped
	return peds
