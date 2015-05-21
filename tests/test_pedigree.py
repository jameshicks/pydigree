from __future__ import division

from pydigree.pedigree import Pedigree
from nose.tools import raises
from testsupport import getpeds

def test_recursivekinship():
    peds = getpeds()
    ped = peds['fullsib']
    
    # Self
    assert ped.kinship('1', '1') == 1/2
    # Unrelated founders 
    assert ped.kinship('1', '2') == 0
    # Full-sib
    assert ped.kinship('3', '4') == 1/4
    # Parent-child
    assert ped.kinship('1', '3') == 1/4
    
    ped = peds['first_cousins']
    # First cousins
    assert ped.kinship('7', '8') == 1/16
    # Avuncular (related)
    assert ped.kinship('7', '5') == 1/8
    # Avuncular (marry-in)
    assert ped.kinship('7', '6') == 0
    # Grandparent
    assert ped.kinship('7', '1') == 1/8

    ped = peds['half_sibs']
    assert ped.kinship('4', '5') == 1/8 

def test_recursivefraternity():
    peds = getpeds()
    ped = peds['fullsib']

    # Unrelated founders                                                                                                                                                                                                                                                                   
    assert ped.fraternity('1', '2') == 0
    # Full-sib                                                                                                                                                                                                                                                                             
    assert ped.fraternity('3', '4') == 1/4
    # Parent-child                                                                                                                                                                                                                                                                         
    assert ped.fraternity('1', '3') == 0

    ped = peds['first_cousins']
    # First cousins                                                                                                                                                                                                                                                                        
    assert ped.fraternity('7', '8') == 0
    # Avuncular (related)                                                                                                                                                                                                                                                                  
    assert ped.fraternity('7', '5') == 0
    # Avuncular (marry-in)                                                                                                                                                                                                                                                                 
    assert ped.fraternity('7', '6') == 0
    # Grandparent                                                                                                                                                                                                                                                                          
    assert ped.fraternity('7', '1') == 0

    ped = peds['half_sibs']
    assert ped.fraternity('4', '5') == 0

def test_inbreeding():
    peds = getpeds()
    ped = peds['first_cousin_child']
    assert ped.inbreeding('9') == 1/16
    # Everyone else is outbred
    for x in ped.individuals:
        if x.id == '9':
            continue
        assert ped.inbreeding(x.id) == 0


@raises(NotImplementedError)
def test_ld():
    ped = Pedigree()
    ped.ld()