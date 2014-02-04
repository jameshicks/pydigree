from itertools import izip, chain, imap

import pydigree


def read_map(mapfile):
    last_chr, last_pos = None, 0
    chroms = []
    chromosome = pydigree.Chromosome()
    with open(mapfile) as f:
        for line in f:
            line = line.strip().split()
            chr, label, cm, pos = line
            if int(pos) < 0:
                continue
            if chr != last_chr:
                chroms.append(chromosome)
                chromosome = pydigree.Chromosome()
            if pos < last_pos:
                raise ValueError('Map file not sorted')
            chromosome.add_genotype(None, cm, label=label, bp=pos)
            last_chr, last_pos = chr, pos
    return chroms


def write_ped(pedigrees, pedfile, mapfile=None, genotypes=True, delim=' ',
              predicate=None):
    """
    write_ped writes data in a plink-format PED file, and optionally a
    plink-format map file.

    Arguments
    ------

    pedigrees: An object of class PedigreeCollection containing what you
        want to output
    pedfile: a string giving the name out the file to output to.
    mapfile: the name of a mapfile to output, if you want to output one.
        an object that evaluates as False or None will skip the mapfile
    genotypes: Should genotypes be output True/False
    delim: Field seperator
    predicate: Which inputs to include in the output file. If not specified
        all are output. If the string is 'affected', only affected
        individuals are output. If the string is 'phenotyped', all individuals
        with phenotype information are output. Any other value of predicate
        must be a function to perform on the individual that evaluates to
        True/False for whether the individual should be output.

    Returns: Nothing
    """
    if not predicate:
        predicate = lambda x: True
    elif predicate == 'affected':
        predicate = lambda x: x.phenotypes['affected'] == 1
    elif predicate == 'phenotyped':
        predicate = lambda x: x.phenotypes['affected'] in set([0, 1])
    elif not callable(predicate):
        raise ValueError('Not a valid predicate!')

    afflab = {1: '2', 0: '1', None: '-9'}

    with open(pedfile, 'w') as f:
        for pedigree in pedigrees:
            for ind in pedigree:
                if not predicate(ind):
                    continue
                # Get the phenotype code
                aff = afflab[ind.phenotypes['affected']]
                # Prepare the 6-column identifier
                outline = [pedigree.label, ind.id,
                           ind.father.id if ind.father is not None else '0',
                           ind.mother.id if ind.mother is not None else '0',
                           1 if ind.sex == 0 else 2,
                           aff]
                # Get the genotypes in the format we need them
                if genotypes:
                    g = []
                    for chroma, chromb in ind.genotypes:
                        g.extend(chain.from_iterable(izip(chroma, chromb)))
                    outline.extend(g)
                # Make strings
                outline = imap(str, outline)
                # Write it out
                f.write(delim.join(outline))
                f.write('\n')

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
