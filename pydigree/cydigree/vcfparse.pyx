from libc.string cimport strsep, strcmp, strlen
from libc.stdlib cimport atoi, free
from libc.stdint cimport int8_t

cimport cython
cimport pydigree.cydigree.datastructures as datastructures

def vcf_allele_parser(datastr, int desired):
    if not datastr:
        return None
    databytes = datastr.encode('utf8')
    cdef char* data = databytes

    cdef int ntok = 0
    cdef int stridx = 0

    while data[stridx] != '\0':
        if data[stridx] == ' ' or data[stridx] == '\t':
            ntok += 1
        stridx += 1
    ntok += 1 # Plus one more token

    cdef datastructures.SparseArray outp = datastructures.SparseArray(ntok * 2, 0)

    cdef char* delim = " \t"
    cdef char* subdelim = ':'
    cdef char* alleledelim = "/|" 
    
    cdef char* geno_tok
    cdef char* allele_tok

    cdef int allele = 0
    cdef int alleleidx = 0
    cdef int subtokidx = 0
    
    cdef char* token = strsep(&data, delim)
    while token:
        subtokidx = 0
        
        if token[0]:
            geno_tok = strsep(&token, subdelim)
            while subtokidx < desired:
                # Fast-forward to the genotype token
                geno_tok = strsep(&token, subdelim)
                subtokidx += 1

            if (strcmp(geno_tok, "0/0") == 0) or (strcmp(geno_tok, "0|0") == 0):
                # Most common scenario
                alleleidx += 2
            
            elif (strcmp(geno_tok, "./.") == 0) or (strcmp(geno_tok, ".|.") == 0):
                # Missing data
                outp.set_item(alleleidx, -1)
                outp.set_item(alleleidx + 1,-1)
                alleleidx += 2 

            elif strlen(geno_tok) == 3:
                allele = <int8_t>geno_tok[0] - 48 # '0' is ascii/utf8 48
                if allele != 0:
                    outp.set_item(alleleidx, allele)
                alleleidx += 1

                allele = <int8_t>geno_tok[2] - 48
                if allele != 0:
                    outp.set_item(alleleidx, allele)
                alleleidx += 1 

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
        
        token = strsep(&data, delim)

    return outp

@cython.boundscheck(False)
def assign_genorow(datastructures.SparseArray row, inds, int chromidx, int markidx):
    cdef datastructures.NodeStack s = row.container.to_stack()
    cdef datastructures.SparseArray spchrom
    cdef int hapidx, indidx

    cdef datastructures.IntTreeNode* denseval = s.pop()
    while denseval:
        indidx = denseval.key // 2
        hapidx = denseval.key % 2

        spchrom = inds[indidx].genotypes[chromidx][hapidx].container
        spchrom.set_item(markidx, denseval.value)

        denseval = s.pop()