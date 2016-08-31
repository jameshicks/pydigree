from libc.string cimport strsep, strdup, strcmp
from libc.stdlib cimport atoi, free
import numpy as np
cimport numpy as np
cimport cython

from pydigree.cydigree.datastructures import SparseArray

def vcf_allele_parser(datastr, int desired, int nallele):
    if not datastr:
        return None
    databytes = datastr.encode('utf8')
    cdef char* data = databytes

    outp = SparseArray(nallele, 0)
    cdef char* datadup = strdup(data)

    cdef char* delim = " \t"
    cdef char* subdelim = ':'
    cdef char* alleledelim = "/|" 
    
    cdef char* geno_tok
    cdef char* allele_tok

    cdef int allele = 0
    cdef int alleleidx = 0
    cdef int subtokidx = 0
    
    cdef char* token = strsep(&datadup, delim)
    while token:
        subtokidx = 0
        
        if token[0]:
            geno_tok = strsep(&token, subdelim)
            while subtokidx < desired:
                # Fast-forward to the genotype token
                geno_tok = strsep(&token, subdelim)
                subtokidx += 1


            if (strcmp(geno_tok, "./.") == 0) or (strcmp(geno_tok, ".|.") == 0):
                # Missing data
                outp.set_item(alleleidx, -1)
                outp.set_item(alleleidx + 1,-1)
                alleleidx += 2 

            elif (strcmp(geno_tok, "0/0") == 0) or (strcmp(geno_tok, "0|0") == 0):
                # Most common scenario
                alleleidx += 2

            else:
                # Allele 1
                allele_tok = strsep(&geno_tok, alleledelim)
                allele = atoi(allele_tok)
                
                if allele != 0:
                    outp.set_item(alleleidx, allele)

                alleleidx += 1

                # Allele 2
                allele_tok = strsep(&geno_tok, alleledelim)
                allele = atoi(allele_tok)

                if allele != 0:
                    outp.set_item(alleleidx, allele)

                alleleidx += 1
        
        token = strsep(&datadup, delim)

    free(datadup)
    return outp
