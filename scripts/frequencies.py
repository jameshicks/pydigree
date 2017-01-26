"A script that calculates allele frequencies from plink data"

import argparse

import pydigree as pyd

parser = argparse.ArgumentParser()
parser.add_argument('--ped', required=True, help='Plink formatted PED file')
parser.add_argument('--map', required=True, help='Plink formatted MAP file')
parser.add_argument('--snps', required=None, nargs='*', metavar='SNP', 
                    default=None, 
                    help='Only calculate frequencies for the specified SNPS')
args = parser.parse_args()

peds = pyd.io.plink.read_plink(pedfile=args.ped, mapfile=args.map)

if args.snps is not None:
    onlysnps = set(args.snps)

def formatted(*cells):
    
    "Print values in tab-delimited format"
    
    return '\t'.join([str(x) for x in cells])

for chromidx, chromobj in enumerate(peds.chromosomes):
    for locidx, markername in enumerate(chromobj.labels):
        if args.snps is not None and markername not in onlysnps:
            continue

        locus = chromidx, locidx
        freqs = list(peds.allele_frequencies(locus).items())
        freqs = sorted(freqs, key=lambda x: x[1], reverse=True)
        maj_allele = freqs[0][0]
        for min_allele, maf in freqs[1:]:
            maf_str = '{:5.4g}'.format(maf)
            print(formatted(chromobj.label, chromobj.physical_map[locidx],
                            chromobj.labels[locidx], maj_allele, min_allele,
                            maf_str))

