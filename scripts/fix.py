#!/usr/bin/env python

import sys

with open(sys.argv[1]) as f:
    for line in f:
        fam,id,fa,mo,sex,aff = line.strip().split()
        if aff == '1': aff='2'
        elif aff == '2': aff='1'
        else: aff ='0'
        print ' '.join([fam,id,fa,mo,sex,aff])
