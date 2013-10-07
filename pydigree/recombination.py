#!/usr/bin/env python

from array import array
from _pydigree import recombine_haldane


def recombine(chr1, chr2, map):
    newchrom = recombine_haldane(chr1, chr2, map)
    if isinstance(chr1, array) and isinstance(chr2, array):
        if chr1.typecode != chr2.typecode:
            raise ValueError('Chromosomes have two different typecodes!')
        newchrom = array(chr1.typecode, newchrom)
    return newchrom
