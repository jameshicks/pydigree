import argparse

import pydigree as pyd

parser = argparse.ArgumentParser()
parser.add_argument('--ped', required=True, help='Plink formatted PED file')
parser.add_argument('--map', required=True, help='Plink formatted MAP file')
args = parser.parse_args()

peds = pyd.io.plink.read_plink(pedfile=args.ped, mapfile=args.map)

def formatted(*cells):
    return '\t'.join([str(x) for x in cells])
for chromidx, chromobj in enumerate(peds.chromosomes):
    for locidx, markername in enumerate(chromobj.labels):
        locus = chromidx, locidx
        freqs = peds.allele_frequencies(locus).items()
        freqs = sorted(freqs, key=lambda x: x[1], reverse=True)
        maj_allele = freqs[0][0]
        for min_allele, maf in freqs[1:]:
            maf_str = '{:5.4g}'.format(maf)
            print formatted(chromobj.label, chromobj.physical_map[locidx],
                            chromobj.labels[locidx], maj_allele, min_allele,
                            maf_str)
