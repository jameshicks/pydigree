import sys
from itertools import takewhile

desired_tokens = {'cdef', 'cpdef', 'ctypedef', '@property'}
spaces_per_level = 4 
with open(sys.argv[1]) as f:
    for line in f:
        if not line.strip():
            continue

        l = line.split()
        
        if l[0] not in desired_tokens:
            continue

        spaces = len(list(takewhile(lambda x: x == ' ', line)))
        if spaces / spaces_per_level < 2:
            print(line.rstrip(' \t\n:'))