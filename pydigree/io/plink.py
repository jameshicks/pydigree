from itertools import izip, chain, imap

from pydigree import Chromosome

def read_map(mapfile):
    last_chr = None
    chroms = []
    chromosome = Chromosome()
    with open(mapfile) as f:
        for line in f:
            line = line.strip().split()
            chr, label, cm, pos = line
            if int(pos) < 0:
                continue
            if chr != last_chr:
                chroms.append(chromosome)
                chromosome = Chromosome()
            chromosome.add_genotype(None, cm, label=label, bp=pos)
    return chroms

def write_ped(pedigrees, pedfile, mapfile=None, genotypes=True, delim=' '):
    with open(pedfile, 'w') as f:
        for pedigree in pedigrees:
            for ind in pedigree:
                aff = {1: '2', 0: '1', None: '-9'}[ind.phenotypes['affected']]
                outline = [pedigree.label, ind.id,
                           ind.father.id if ind.father is not None else '0',
                           ind.mother.id if ind.mother is not None else '0',
                           1 if ind.sex == 0 else 2,
                           aff]
                if genotypes:
                    g = []
                    for chroma, chromb in ind.genotypes:
                        g.extend(chain.from_iterable(izip(chroma, chromb)))
                    outline.extend(g)
                outline = imap(str, outline)
                f.write(delim.join(outline) + '\n')
    if not mapfile:
        return
    with open(mapfile, 'w') as f:
        for ci, chromosome in enumerate(pedigrees.chromosomes):
            for mi, marker in enumerate(chromosome._iinfo):
                label, cm, mb, frequency = marker
                if not mb:
                    mb = int(cm * 10e6)
                if not label:
                    label = 'SNP%s-%s' % (chromosome, mi)
                f.write('\t'.join(str(x) for x
                                  in [chromosome, label, cm, mb]) + '\n')
