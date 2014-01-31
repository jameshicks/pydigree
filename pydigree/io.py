from itertools import izip, chain, imap

from pydigree.common import *
from pydigree.population import Population
from pydigree.individual import Individual
from pydigree.pedigree import Pedigree
from pydigree.pedigreecollection import PedigreeCollection
from pydigree.chromosome import Chromosome


def read_ped(filename, population=None, delimiter=None, affected_labels=None):
    """
    Reads a plink format pedigree file, ie:
        familyid indid father mother sex whatever whatever whatever
    into a pydigree pedigree object, with optional population to
    assign to pedigree members. If you don't provide a population
    you can't simulate genotypes!

    Arguements
    -----
    filename: The file to be read
    population: The population to assign individuals to
    delimiter: a string defining the field separator, default: any whitespace
    affected_labels: The labels that determine affection status.

    Returns: An object of class PedigreeCollection
    """
    p = Pedigree()
    pc = PedigreeCollection()
    sex_codes = {'1': 0, '2': 1, 'M': 0, 'F': 1}
    if not affected_labels:
        affected_labels = {'1': 0, '2': 1,
                           'A': 1, 'U': 0,
                           'X': None,
                           '-9': None}

    def getph(ph):
        try:
            return affected_labels[ph]
        except KeyError:
            return None

    with open(filename) as f:
        for line in f:
            fam, id, fa, mo, sex, aff = line.strip().split(delimiter)[0:6]
            # Give a special id for now, to prevent overwriting duplicated
            # ids between families
            id = (fam, id)
            p[id] = Individual(population, id, fa, mo, sex)
            p[id].phenotypes['affected'] = getph(aff)
            p[id].pedigree = p
        for ind in p:
            fam, id = ind.id
            # Actually make the references instead of just pointing at strings
            ind.father = p[(fam, ind.father)] if ind.father != '0' else None
            ind.mother = p[(fam, ind.mother)] if ind.mother != '0' else None
            ind.sex = sex_codes[ind.sex]
            ind.register_with_parents()
        for pedigree_label in set(ind.id[0] for ind in p):
            ped = Pedigree(label=pedigree_label)
            thisped = [x for x in p if x.id[0] == pedigree_label]
            for ind in thisped:
                ind.id = ind.id[1]
                ped[ind.id] = ind
                ind.population = ped
            pc[pedigree_label] = ped
    return pc


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


def read_gs_chromosome_template(templatef):
    with open(templatef) as f:
        f.readline()  # The label and
        f.readline()  # the number of markers, both of which we dont need.
        c = Chromosome()
        last_cm = 0
        for line in f:
            if line == '\n':
                continue
            label, majf, minf, cm, bp = line.strip().split()
            cm = float(cm)
            last_cm += cm
            c.add_genotype(float(minf), last_cm, label=label, bp=bp)
    return c


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


def read_phenotypes(pedigrees, csvfile, delimiter=',', missingcode='X'):
    """
    Reads a csv with header
    famid,ind,phen,phen,phen,phen etc etc

    Arguements
    ------
    Pedigrees:   An object of class PedigreeCollection
    csvfile:     the filename of the file containing phenotypes.
    delimiter:   the field delimiter for the file
    missingcode: the code for missing values

    Returns: Nothing
    """
    with open(csvfile) as f:
        header = f.readline().strip().split(delimiter)
        for line in f:
            # Match columns to their column name
            d = dict(zip(header, line.strip().split(delimiter)))
            for k, v in d.items():
                # Convert all phenotypes into floats
                try:
                    v = float(v)
                except ValueError:
                    if not v:
                        v = None
                if k in set(['famid', 'ind']):
                    continue
                fam, ind = d['famid'], d['ind']
                pedigrees[fam][ind].phenotypes[k] = v
