from libc.stdlib cimport atoi
from libc.stdio cimport sscanf

def read_map_to_nested_list(mapfile):
    def parseline(line):
        chrom, snp, cm, pos = line.strip().split()
        cm = float(cm)
        pos = int(pos)
        return [chrom, snp, cm, pos]
    with open(mapfile) as f:
        return [parseline(x) for x in f]

cdef int haplocheck(a):
    return a.endswith('.0') or a.endswith('.1')

def read_germline(matchfile, mapfile, int haploid=0,
    double markerdensitylimit=1e-4, double minsegmentlength=500000):

    cdef unsigned long markersinshare
    cdef double density
    cdef int start, stop

    # Get map info 
    cmap = read_map_to_nested_list(mapfile)
    positions = (x[3] for x in cmap)
    posd =  dict([(y,x) for x,y in enumerate(positions)])

    with open(matchfile) as sharef:
        shared = {}
        
        # Test the first line to see if we're in a haploid file
        line = sharef.readline()
        l = line.strip().split()

        if haploid and not haplocheck(l[1]):
            raise ValueError('Not a `germline --haploid` formatted file')

        sharef.seek(0)
        for line in sharef:
            l = line.strip().split()
            print l
            # Germline haploid output gives which haplotype is shared.
            # If we're just counting up (like we would for spairs) we
            # can just lop off the haplotype identifiers, and the algorithm
            # will just add another 1 when it comes over an overlapping region
            if haploid:
                l[1] = l[1][:-2]
                l[3] = l[3][:-2]

            ind1 = l[0], l[1]
            ind2 = l[2], l[3]
            
            pair = frozenset( (ind1,ind2) )

            start = atoi(l[5])
            stop = atoi(l[6])

            istart = posd[start]
            istop = posd[stop]
            markersinshare = istop - istart + 1 # Both bounds are shared
            density = markersinshare / <double>(stop - start)

            if (stop - start) < minsegmentlength:
                continue

            if density < markerdensitylimit:
                continue

            segment = istart, istop
            try:
                shared[pair].append(segment)
            except KeyError:
                shared[pair] = [segment]

    return shared

def read_kinship(filename):
    kindict = {}
    with open(filename) as f:
        for line in f:
            fam, ida, idb, phi = line.strip().split()
            kindict[frozenset({ (fam, ida), (fam, idb)})] = float(phi)
    return kindict
