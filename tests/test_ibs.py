from pydigree import ibs


def test_ibs():
    assert ibs( (2,2), (2,2) ) == 2
    assert ibs( (2,1), (2,2) ) == 1
    assert ibs( (1,1), (2,2) ) == 0
    assert ibs( (0,0), (1,1) ) == None
