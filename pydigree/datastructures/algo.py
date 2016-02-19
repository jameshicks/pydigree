def binsearch(x, a, lo=None, hi=None, key=None):
    """ 
    A binary search that allows for arbitrary transformations with a function
    object 'key'. 

    Arguments:
    x: a list/array structure
    a: the value (after transformation) searched for
    key: a function to transform the list elements

    Returns: The index of a
    Raises: KeyError if key(a) not in x
    """
    if not x:
        raise KeyError

    if key is None:
        key = lambda x: x

    if lo is not None:
        idx_min = lo
    else:
        idx_min = 0 
    
    if hi is not None:
        idx_max = hi
    else:
        idx_max = len(x) - 1 

    
    while idx_min <= idx_max:
        idx_test = int((idx_min + idx_max) / 2)
        value_at_index = key(x[idx_test])

        if value_at_index > a: 
            idx_max = idx_test - 1
        elif value_at_index < a:
            idx_min = idx_test + 1
        elif value_at_index == a: 
            return idx_test
    raise KeyError

def binsearch_left(x, a, lo=None, hi=None, key=None):
    "Find the rightmost insertion point for a in x where sort is maintained"
    if not x:
        return 0

    if key is None:
        key = lambda x: x

    if lo is not None:
        idx_min = lo
    else:
        idx_min = 0 
    
    if hi is not None:
        idx_max = hi
    else:
        idx_max = len(x) - 1 

    
    while idx_min <= idx_max:
        idx_test = int((idx_min + idx_max) / 2)
        value_at_index = key(x[idx_test])

        if value_at_index >= a: 
            idx_max = idx_test - 1
        else:
            idx_min = idx_test + 1
    return idx_min


def binsearch_right(x, a, lo=None, hi=None, key=None):
    "Find the rightmost insertion point for a in x where sort is maintained"
    if not x: 
        return 0

    if key is None:
        key = lambda x: x

    if lo is not None:
        idx_min = lo
    else:
        idx_min = 0 
    
    if hi is not None:
        idx_max = hi
    else:
        idx_max = len(x) - 1 

    
    while idx_min <= idx_max:
        idx_test = int((idx_min + idx_max) / 2)
        value_at_index = key(x[idx_test])

        if value_at_index > a: 
            idx_max = idx_test - 1
        else:
            idx_min = idx_test + 1

    return idx_min