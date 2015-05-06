def read_kinship(filename):
    kindict = {}
    with open(filename) as f:
        for line in f:
            fam, ida, idb, phi = line.strip().split()
            kindict[frozenset({(fam, ida), (fam, idb)})] = float(phi)
    return kindict
