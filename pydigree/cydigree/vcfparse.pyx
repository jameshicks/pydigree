from libc.string cimport strsep, strcmp, strlen
from libc.stdlib cimport atoi, free
from libc.stdint cimport int8_t, uint32_t
from cpython.mem cimport PyMem_Malloc, PyMem_Free

cimport cython
cimport pydigree.cydigree.datastructures as datastructures

cdef struct VariantCall:
    uint32_t alleleidx
    int8_t allele
    VariantCall* following 

cdef class VariantStack(object):
    cdef VariantCall* front
    def __cinit__(self):
        self.front = NULL

    cdef void push(self, uint32_t alleleidx, int8_t allele):
        cdef VariantCall* var = <VariantCall*>PyMem_Malloc(sizeof(VariantCall))
        var.alleleidx = alleleidx
        var.allele = allele
        var.following = self.front
        self.front = var

    cdef VariantCall* pop(self):
        if not self.front:
            return NULL
        cdef VariantCall* item = self.front
        self.front = item.following
        return item

    def tolist(self):
        cdef VariantCall* item = self.front

        outp = []
        while item:
            outv = (item.alleleidx, item.allele)
            outp.append(outv)
            item = item.following

        return outp

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

    cdef VariantStack outstack = VariantStack()

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
                outstack.push(alleleidx, -1)
                outstack.push(alleleidx + 1,-1)
                alleleidx += 2 

            elif strlen(geno_tok) == 3:
                allele = <int8_t>geno_tok[0] - 48 # '0' is ascii/utf8 48
                if allele != 0:
                    outstack.push(alleleidx, allele)
                alleleidx += 1

                allele = <int8_t>geno_tok[2] - 48
                if allele != 0:
                    outstack.push(alleleidx, allele)
                alleleidx += 1 

            else:
                # Allele 1
                allele_tok = strsep(&geno_tok, alleledelim)
                allele = atoi(allele_tok)
                
                if allele != 0:
                    outstack.push(alleleidx, allele)

                alleleidx += 1

                # Allele 2
                allele_tok = strsep(&geno_tok, alleledelim)
                allele = atoi(allele_tok)

                if allele != 0:
                    outstack.push(alleleidx, allele)

                alleleidx += 1
        
        token = strsep(&data, delim)

    return outstack

@cython.boundscheck(False)
def assign_genorow(VariantStack row, inds, int chromidx, int markidx):
    cdef datastructures.SparseArray spchrom
    cdef int hapidx, indidx

    cdef VariantCall* denseval = row.pop()
    while denseval:
        indidx = denseval.alleleidx // 2
        hapidx = denseval.alleleidx % 2

        spchrom = inds[indidx].genotypes[chromidx][hapidx].container
        spchrom.set_item(markidx, denseval.allele)

        PyMem_Free(denseval)
        denseval = row.pop()
