from pydigree import product, cumsum


def test_product():
    assert product([1,2,3,4,5]) == 120
    assert product([0,1,2,3,4]) ==  0


def test_cumsums():
    assert cumsum([]) == []
    assert cumsum([1,2,3]) == [1,3,6]

