from __future__ import division

import os
import glob
from pydigree.io import read_ped
from testsupport import getpeds


def test_commonancestors():
    from pydigree.paths import common_ancestors
    peds = getpeds()
    ped = peds['fullsib']
    # Two founders have no common ancestors
    assert common_ancestors(ped['1'], ped['2']) == set()

    # Parent and child no common ancestors
    assert common_ancestors(ped['1'], ped['3']) == set()

    # Full sibs have parents in common
    assert common_ancestors(ped['3'], ped['4']) == {ped['1'], ped['2']}
    
    # Test first cousins
    ped = peds['first_cousins']
    assert common_ancestors(ped['7'], ped['8']) == {ped['1'], ped['2']}

def test_pathkinship():
    from pydigree.paths import kinship
    peds = getpeds()
    ped = peds['fullsib']
    
    # None should return 0
    assert kinship(ped['1'], None) == 0
    # Self
    assert kinship(ped['1'], ped['1']) == 1/2
    # Unrelated founders 
    assert kinship(ped['1'], ped['2']) == 0
    # Full-sib
    assert kinship(ped['3'], ped['4']) == 1/4
    # Parent-child
    assert kinship(ped['1'], ped['3']) == 1/4
    
    ped = peds['first_cousins']
    # First cousins
    assert kinship(ped['7'], ped['8']) == 1/16
    # Avuncular (related)
    assert kinship(ped['7'], ped['5']) == 1/8
    # Avuncular (marry-in)
    assert kinship(ped['7'], ped['6']) == 0
    # Grandparent
    assert kinship(ped['7'], ped['1']) == 1/8

    ped = peds['half_sibs']
    assert kinship(ped['4'], ped['5']) == 1/8 

def test_pathfraternity():
    from pydigree.paths import fraternity
    peds = getpeds()
    ped = peds['fullsib']

    # None should return 0                                                                                                                                                                                                                                                                 
    assert fraternity(ped['1'], None) == 0
    # Unrelated founders                                                                                                                                                                                                                                                                   
    assert fraternity(ped['1'], ped['2']) == 0
    # Full-sib                                                                                                                                                                                                                                                                             
    assert fraternity(ped['3'], ped['4']) == 1/4
    # Parent-child                                                                                                                                                                                                                                                                         
    assert fraternity(ped['1'], ped['3']) == 0

    ped = peds['first_cousins']
    # First cousins                                                                                                                                                                                                                                                                        
    assert fraternity(ped['7'], ped['8']) == 0
    # Avuncular (related)                                                                                                                                                                                                                                                                  
    assert fraternity(ped['7'], ped['5']) == 0
    # Avuncular (marry-in)                                                                                                                                                                                                                                                                 
    assert fraternity(ped['7'], ped['6']) == 0
    # Grandparent                                                                                                                                                                                                                                                                          
    assert fraternity(ped['7'], ped['1']) == 0

    ped = peds['half_sibs']
    assert fraternity(ped['4'], ped['5']) == 0