def write_sgs(data, filename):
    with open(filename, 'w') as o:
        for pair, intervals in data.items():
            for i, chromintervals in enumerate(intervals):
                if not chromintervals:
                    continue
                ind1, ind2 = pair
                l = [ind1.pedigree.label, ind1.label,
                     ind2.pedigree.label, ind2.label, i] + chromintervals
                o.write('\t'.join(l) + '\n')
