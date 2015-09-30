#!/usr/bin/env python

import pydigree
import itertools

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--file', required=True,
                    help='Pedigree file for kinship/inbreeding calculation')
args = parser.parse_args()

ped = pydigree.io.read_ped(args.file)

for pedigree in ped:
    lab = pedigree.label
    ids = sorted([x.label for x in pedigree])
    # itertools combinations with replacement isn't available in 2.6
    # So we'll do this as a two step process
    # 1) Pairwise kinship coefs
    for x, y in itertools.combinations(ids, 2):
        k = pedigree.kinship(x, y)
        print lab, x, y, k
    # 2) Inbreeding coefs.
    for x in ids:
        print lab, x, x, pedigree.inbreeding(x)
