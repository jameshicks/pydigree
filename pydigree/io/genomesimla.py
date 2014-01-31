import pydigree

def read_gs_chromosome_template(templatef):
    with open(templatef) as f:
        f.readline()  # The label and
        f.readline()  # the number of markers, both of which we dont need.
        c = pydigree.Chromosome()
        last_cm = 0
        for line in f:
            if line == '\n':
                continue
            label, majf, minf, cm, bp = line.strip().split()
            cm = float(cm)
            last_cm += cm
            c.add_genotype(float(minf), last_cm, label=label, bp=bp)
    return c

