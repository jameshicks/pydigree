def read_map_to_nested_list(mapfile):
    def parseline(line):
        chrom, snp, cm, pos = line.strip().split()
        cm = float(cm)
        pos = int(pos)
        return [chrom, snp, cm, pos]
    with open(mapfile) as f:
        return [parseline(x) for x in f]

def haplocheck(a):
    return a.endswith('.0') or a.endswith('.1')

def read_germline(matchfile, mapfile, haploid=0,
        markerdensitylimit=1e-4, minsegmentlength=500000):

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
        for i,line in enumerate(sharef):
            l = line.strip().split()

            # Germline haploid output gives which haplotype is shared.
            # If we're just counting up (like we would for spairs) we
            # can just lop off the haplotype identifiers, and the algorithm
            # will just add another 1 when it comes over an overlapping region
            if haploid:
                l[1] = l[1][:-2]
                l[3] = l[3][:-2]

            ind1 = tuple(l[0:2])
            ind2 = tuple(l[2:4])
            
            pair = frozenset({ind1,ind2})

            start = int(l[5])
            stop = int(l[6])

            istart = posd[start]
            istop = posd[stop]
            markersinshare = istop - istart
            density = markersinshare / (stop - start)

            if (stop - start) < minsegmentlength:
                continue

            if density < markerdensitylimit:
                continue

            if pair not in shared:
                shared[pair] = []
            
            shared[pair].append((istart, istop))

    return shared
