#include <Python.h>
#include <math.h>

/* Function Prototypes */ 

PyObject* haldane_interface(PyObject* self, PyObject* args);
PyObject* recombine_haldane(PyObject* A, PyObject* B, PyObject* map, Py_ssize_t l);

PyObject* sample_with_replacement_interface(PyObject* self, PyObject* args);
PyObject* sample_with_replacement(PyObject* population, long int n);

PyObject* random_pairs_interface(PyObject* self, PyObject* args);
PyObject* random_pairs(PyObject* pop, long int numpairs);

PyObject* choice_probs_interface(PyObject* self, PyObject* args);
PyObject* choice_probs(PyObject* choices, PyObject* probs);

PyObject* chromatid_delabeler_interface(PyObject* self, PyObject* args);
PyObject* chromatid_delabeler(PyObject* chromatid, Py_ssize_t chromidx);

PyObject* linkeq_chrom_interface(PyObject* self, PyObject* args);
PyObject* linkeq_chrom(PyObject *frequencies);

PyObject* sgs_shares(PyObject* affecteds, PyObject* shared, Py_ssize_t nmark); 
PyObject* sgs_shares_interface(PyObject* self, PyObject* args);

/* docstrings */
static char module_docstring[] =
  "This sub-module provides C implementations of time-intensive pydigree functions";

static char haldane_docstring[] = "Takes 3 python lists, chromA, chromB, and map. The map is a list of floats indicating\
 the map position";
static char sample_repl_docstring[] = "Randomly choses n individuals (with replacement) from the population";
static char pairs_docstring[] = "Returns a list of n pairs from sampled with replacement from population"; 
static char choice_prob_docstring[] = "Chooses a random item based on probabilities provided in probs"; 
static char chromatid_delabeler_docstring[] = "Replaces genotypes from label_genotypes with alleles from the ancestor"; 
static char linkeq_chrom_docstring[] = "Returns a chromosome with all markers in linkage equilibrium"; 
static char sgs_shares_docstring[] = "Returns the proportion of individual pairs IBD";
/* Python C API boilerplate */
static PyMethodDef module_methods[] = {
  {"recombine_haldane", haldane_interface, METH_VARARGS, haldane_docstring},
  {"sample_with_replacement", sample_with_replacement_interface, METH_VARARGS, sample_repl_docstring},
  {"random_pairs",random_pairs_interface,METH_VARARGS,pairs_docstring},
  {"choice_with_probs", choice_probs_interface, METH_VARARGS, choice_prob_docstring},
  {"chromatid_delabeler", chromatid_delabeler_interface, METH_VARARGS, chromatid_delabeler_docstring},
  {"linkeq_chrom", linkeq_chrom_interface, METH_VARARGS, linkeq_chrom_docstring},
  {"sgs_shares", sgs_shares_interface, METH_VARARGS, sgs_shares_docstring},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_pydigree(void)
{
  (void) Py_InitModule("pydigree._pydigree", module_methods);
  srand(time(NULL));
}


/* Common functions */ 
int rand_lim(int limit) {
  /* return a random number between 0 (inclusive) and limit (exclusive).
   */
  int divisor = RAND_MAX/(limit);
  int retval;

  do { 
    retval = lrand48() / divisor;
  } while (retval > limit);

  return retval;
}

int rand_int_range(int min, int max) {
  return rand_lim(max) + min; 
}

int rand_bool(double prob) { 
  return drand48() < prob; 
}

double rexp(double rate) { 
  return log(1-drand48()) / (-rate);
}



PyObject* haldane_interface(PyObject* self, PyObject* args) {
  PyObject *chromA,*chromB, *map;
  if (!PyArg_ParseTuple(args, "OOO", &chromA, &chromB, &map)) return NULL;
  Py_ssize_t length = PySequence_Size(chromA);
  PyObject* nc = recombine_haldane(chromA,chromB,map,length);
  return nc;

}

PyObject* recombine_haldane(PyObject* A, PyObject* B, PyObject* map, Py_ssize_t l) {
  Py_ssize_t i;
  unsigned char flipped = (lrand48() & 1) ? 1 : 0;
  double recombination_site = rexp(0.01);
  PyObject* newchrom = PyList_New(l);
  for (i = 0; i < l; i++) {

    /* Do we cross over? */ 
    double position = PyFloat_AsDouble( PyList_GetItem(map,i) );
    if (recombination_site < position) {
      flipped = !flipped;
      recombination_site = rexp(0.01);
    }
    
    /* Put it on the new chromatid */
    /*
       FIXME: check if I'm leaking memory here. PyList used to return a
       borrowed reference for GetItem, but PySequence returns a new reference!
    */

    PyObject* allele = PySequence_GetItem(flipped ? B : A, i);
    Py_INCREF(allele);
    PyList_SET_ITEM(newchrom,i,allele);

  }
  return newchrom;
}

/* sampling with replacement */

PyObject* sample_with_replacement_interface(PyObject* self, PyObject* args) {
  PyObject* pop;
  long int n;
  if (!PyArg_ParseTuple(args, "Ol", &pop, &n)) return NULL;
  PyObject* ns = sample_with_replacement(pop,n);
  return ns; 
}


PyObject* sample_with_replacement(PyObject* population, long int n) {
  /* TODO: figureout my long int situation */ 
  PyObject* sample = PyList_New(n); 
  Py_ssize_t population_size = PyList_Size(population);   
  Py_ssize_t i; 
  for (i=0;i < n; i++) {
    long int r = rand_int_range(0,(int)population_size-1); 
    PyObject* draw = PyList_GetItem(population,r);
    PyList_SET_ITEM(sample,i,draw);
    Py_INCREF(draw);
  }
  return sample;
}

/* Randomly pairing items (with replacement) */

PyObject* random_pairs_interface(PyObject* self, PyObject* args) {
  PyObject* pop; 
  long int numpairs; 
  if (!PyArg_ParseTuple(args, "Ol", &pop, &numpairs)) return NULL;
  PyObject* newobj = random_pairs(pop,numpairs);
  return newobj;
}

PyObject* random_pairs(PyObject* pop, long int numpairs) {
  PyObject* pairs = PyList_New(numpairs);
  Py_ssize_t i;
  for (i=0;i<numpairs;i++) {
    PyObject* pair = sample_with_replacement(pop,2);
    pair = PyList_AsTuple(pair);
    PyList_SET_ITEM(pairs,i,pair);
    Py_INCREF(pair);
  }
  return pairs;
}


/* Choice with probabilities */ 

PyObject* choice_probs_interface(PyObject* self, PyObject* args) {
  PyObject *choices, *probs; 
  if (!PyArg_ParseTuple(args, "OO", &choices, &probs)) return NULL;
  PyObject* newobj = choice_probs(choices,probs);
  return newobj; 
}

PyObject* choice_probs(PyObject* choices, PyObject* probs) {
  Py_ssize_t n_choices = PyList_Size(choices); 
  double cum_prob = 0;
  Py_ssize_t i;
  for (i=0;i<n_choices;i++) {
    cum_prob += PyFloat_AsDouble(PyList_GetItem(probs,i));
    if (drand48() < cum_prob) {
      PyObject* choice = PyList_GetItem(choices,i);
      Py_INCREF(choice);
      return choice;
    }
  }
  return NULL;
}

PyObject* chromatid_delabeler_interface(PyObject* self, PyObject* args) {
  PyObject *chromatid;
  Py_ssize_t chromidx;
  if (!PyArg_ParseTuple(args, "On", &chromatid, &chromidx)) return NULL;
  PyObject* newchromatid = chromatid_delabeler(chromatid, chromidx);
  return newchromatid;
}

PyObject* chromatid_delabeler(PyObject* chromatid, Py_ssize_t chromidx) {
  Py_ssize_t chromlength = PyList_Size(chromatid);
  PyObject* newchromatid = PyList_New(chromlength);

  Py_ssize_t i;
  PyObject* last_ancestor = NULL;
  Py_ssize_t last_ah = NULL;
  PyObject *anc_genotypes = NULL;
  PyObject *anc_chromosomes = NULL;
  PyObject *anc_chromosome = NULL;
  for (i=0; i < chromlength; i++) {
    PyObject* listitem = PyList_GET_ITEM(chromatid, i);
    PyObject* ancestor = PyTuple_GET_ITEM(listitem, 0);
    Py_ssize_t ancestral_hap = PyInt_AsSsize_t(PyTuple_GetItem(listitem, 1));
    
    if (!(ancestor == last_ancestor && ancestral_hap == last_ah)) {
      /* Get the genotypes of the ancestor, which are a list of pairs of lists */
      anc_genotypes = PyObject_GetAttrString(ancestor, "genotypes");
      anc_chromosomes = PySequence_GetItem(anc_genotypes, chromidx);
      anc_chromosome = PySequence_GetItem(anc_chromosomes, ancestral_hap);
    }

    PyObject* allele = PySequence_GetItem(anc_chromosome, i);
    PyList_SET_ITEM(newchromatid, i, allele);
    Py_INCREF(allele);

    last_ancestor = ancestor;
    last_ah = ancestral_hap;
  }
  return newchromatid;
}

PyObject* linkeq_chrom_interface(PyObject* self, PyObject* args) {
  PyObject* frequencies;
  if (!PyArg_ParseTuple(args, "O", &frequencies)) return NULL;
  PyObject* chrom = linkeq_chrom(frequencies);
  return chrom;
}

PyObject* linkeq_chrom(PyObject* frequencies) {
  Py_ssize_t chromlen = PySequence_Size(frequencies);
  PyObject* newchrom = PyList_New(chromlen);
  Py_ssize_t marker; 

  for (marker = 0; marker < chromlen; marker++) {
    PyObject* freq = PySequence_GetItem(frequencies, marker);
    double f = PyFloat_AS_DOUBLE(freq);
    PyObject* allele = PyInt_FromLong(rand_bool(f) ? 1 : 2);
    PyList_SET_ITEM(newchrom, marker, allele);
    Py_INCREF(allele);
  }

  return newchrom;
}

PyObject* sgs_shares_interface(PyObject* self, PyObject* args) {
  PyObject *shared, *affecteds;
  Py_ssize_t nmark;
  if (!PyArg_ParseTuple(args, "OOn", &affecteds, &shared, &nmark)) return NULL;
  return sgs_shares(affecteds, shared, nmark);

}

PyObject* sgs_shares(PyObject* affecteds, PyObject* shared, Py_ssize_t nmark) {
  Py_ssize_t i,j,k,m, nsegs, start, stop; 
  PyObject *a, *b, *pairtuple, *pair;
  PyObject *sharelocs, *loc, *nshareslist;

  Py_ssize_t naff = PyList_Size(affecteds);
  /* Initialize empty array */ 
  unsigned long int shares[nmark];
  for (i = 0; i < nmark; i++) {
    shares[i] = 0;
  }

  /* Start incrementing */
  for (i=0; i < naff; i++) {
    a = PyList_GET_ITEM(affecteds, i);
    for (j=i+1; j < naff; j++) {
      b = PyList_GET_ITEM(affecteds, j); 
      
      /* Get shares between individuals */
      pairtuple = PyTuple_Pack(2, a, b);
      pair = PyFrozenSet_New(pairtuple);

      sharelocs = PyDict_GetItem(shared, pair);
      if (sharelocs == NULL) {
	/* I'm expecting a KeyError here but in the case of anything
	   I'm just going to skip it anyway
	*/
	PyErr_Clear();
	goto finish;
      }

      /* For each share, increment the count inside it */
      nsegs = PyList_Size(sharelocs);
      for (k=0; k < nsegs; k++) {
	loc = PyList_GET_ITEM(sharelocs, k);
	start = PyInt_AsSsize_t(PyTuple_GET_ITEM(loc, 0));
	stop = PyInt_AsSsize_t(PyTuple_GET_ITEM(loc, 1));
	for (m=start; m <= stop; m++) {
	  shares[m] += 1;
	}
      }
    finish:
      Py_XDECREF(pair);
      Py_XDECREF(pairtuple);
    }
  }
  /* Make into a python list */ 
  nshareslist = PyList_New(nmark);
  for (i=0; i < nmark; i++) {
    PyList_SET_ITEM(nshareslist, i, PyLong_FromUnsignedLong(shares[i]));
  }
  return nshareslist;
}
