from itertools import chain

import numpy as np

from pydigree.io import smartopen
from pydigree.io.base import genotypes_from_sequential_alleles
from pydigree.genotypes import ChromosomeTemplate, ChromosomeSet
from pydigree.exceptions import FileFormatError
import collections
from functools import reduce


def create_pop_handler_func(chromosomes):

    def pop_handler(pop):
        for chrom in chromosomes:
            pop.add_chromosome(chrom)

    return pop_handler


def create_data_handler_func(chromosomes,
                             phenotypes, phenotype_indices,
                             genotypes, genotype_indices):

    def data_handler(ind, data):

        # Handle the phenotypes first
        phendata = data[phenotype_indices]
        for phenotype, value in zip(phenotypes, phendata):
            kind, label = phenotype
            if kind == 'A':
                value = True if value == '1' else False
            else:
                value = float(value) if '.' in value else int(value)
            ind.phenotypes[label] = value

        gts = data[genotype_indices]
        ind.genotypes = genotypes_from_sequential_alleles(chromsomes,
                                                          gts,
                                                          missing_code='X')


def read_dat(filename):
    with open(filename) as f:
        return [tuple(line.strip().split()) for line in f]


def read_merlin(prefix=None, pedfile=None, datfile=None, mapfile=None):
    field_types = read_dat(datfile)
    data_marker_labels = [lab for kind, lab in field_types if kind == 'M']

    # We only want the map positions for the genotypes in the dat file
    chromosomes = read_map(mapfile, only=data_marker_labels)

    pop_handler = create_pop_handler_func(chromosomes)

    maplabs = list(chain(x.labels for x in chromosomes))

    # Require that all markers in the data file be in the map
    if set(data_marker_labels) ^ set(maplabs):
        raise FileFormatError('Datfile and map file contain different markers')

    # Make sure that the order of markers in the dat file is the same as that
    # in the map file.
    for datlab, maplab in zip(data_marker_labels, maplabs):
        if datlab != maplab:
            raise FileFormatError('Datfile and mapfile in different order')

    phenidx = phenotype_indices([kind for kind, label in field_types])
    genidx = genotype_indices([kind for kind, label in field_types])
    pop_handler = create_pop_handler_func(chromosomes)
    data_handler = create_data_handler_func(chromosomes,
                                            phenotypes, phenidx,
                                            genotypes, gtidx)


def phenotype_indices(types):
    indices = []
    for property_type in types:
        if property_type in {'S2', 'M'}:
            # Genotypes (M) or skipped genotypes (S2) take up two columns
            indices.extend([False, False])
        elif property_type == 'S':
            indices.append(False)
        elif property_type in {'A', 'T', 'C'}:
            indices.append(True)
        else:
            raise FileFormatError(
                'Unknown MERLIN field type: {}'.format(property_type))
    return np.array(indices, dtype=np.bool)


def genotype_indices(types):
    indices = []
    for property_type in types:
        if property_type == 'M':
            # Genotypes (M) take up two columns
            indices.extend([True, True])
        elif property_type == 'S2':
            # Skipped genotypes take up two columns
            indices.extend([False, False])
        elif property_type in {'A', 'S', 'T', 'C'}:
            indices.append(False)
        else:
            raise FileFormatError(
                'Unknown MERLIN field type: {}'.format(property_type))
    return np.array(indices, dtype=np.bool)


def read_map(filename, only=None):
    with open(filename) as f:
        line = f.readline()
        header = line.strip().split()
        sex_specific_map = len(header) == 5

        lastchrom = None
        chromosomes = ChromosomeSet()
        for markernum, line in enumerate(f):
            l = line.strip().split()

            if sex_specific_map:
                chrom, marker, pos, male_pos, female_pos = l
            else:
                chrom, marker, pos = l

            if only and marker not in only:
                continue

            pos = float(pos)

            if chrom != lastchrom:
                if lastchrom:
                    c.finalize()
                c = ChromosomeTemplate(label=chrom)
                chromosomes.add_chromosome(c)
                lastchrom = chrom

            c.add_genotype(map_position=pos, label=marker)

        for chrom in chromosomes:
            chrom.finalize()
        return chromosomes


def write_phenotypes(pedigrees, prefix, phenotypes=None,
                     phenotype_types=None, predicate=None,
                     skip_all_missing=True, missing_code='X', trait=None):
    output_inds = pedigrees.individuals
    if isinstance(predicate, collections.Callable):
        output_inds = [x for x in output_inds if predicate(x)]

    if phenotypes is None:
        available_phenotypes = reduce(set.union,
                                      [set(x.phenotypes.keys())
                                       for x in pedigrees.individuals])
        available_phenotypes = sorted(available_phenotypes)
    else:
        available_phenotypes = phenotypes

    if not phenotype_types:
        # Default to covariate
        phenotype_types = [('T' if ph == trait else 'C')
                           for ph in available_phenotypes]

    if len(phenotype_types) != len(available_phenotypes):
        raise ValueError('Not enough types for specified phenotypes')

    def get_ph(ind, phen):
        if phen not in ind.phenotypes:
            return str(missing_code)

        p = ind.phenotypes[phen]

        if p is None:
            return str(missing_code)

        return str(ind.phenotypes[phen])

    with open(prefix + '.phenotypes.dat', 'w') as datf:
        for ptype, phen in zip(phenotype_types, available_phenotypes):
            datf.write('{0} {1}\n'.format(ptype, phen))

    with open(prefix + '.phenotypes.ped', 'w') as phenf:
        for ind in pedigrees.individuals:
            pedinfo = [str(x) for x in (ind.population.label,
                                        ind.label,
                                        ind.father.label if ind.father is not None else '0',
                                        ind.mother.label if ind.mother is not None else '0',
                                        1 if ind.sex == 0 else 2)]
            phenotypes = [get_ph(ind, ph) for ph in available_phenotypes]

            if all(ph == missing_code for ph in phenotypes):
                continue

            oline = ' '.join(str(x) for x in pedinfo + phenotypes)
            phenf.write(oline + '\n')
