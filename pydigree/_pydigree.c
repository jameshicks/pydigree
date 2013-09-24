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


/* docstrings */
static char module_docstring[] =
  "This sub-module provides C implementations of time-intensive pydigree functions";

static char haldane_docstring[] = "Takes 3 python lists, chromA, chromB, and map. The map is a list of floats indicating\
 the map position";
static char sample_repl_docstring[] = "Randomly choses n individuals (with replacement) from the population";
static char pairs_docstring[] = "Returns a list of n pairs from sampled with replacement from population"; 
static char choice_prob_docstring[] = "Chooses a random item based on probabilities provided in probs"; 

/* Python C API boilerplate */
static PyMethodDef module_methods[] = {
  {"recombine_haldane", haldane_interface, METH_VARARGS, haldane_docstring},
  {"sample_with_replacement", sample_with_replacement_interface, METH_VARARGS, sample_repl_docstring},
  {"random_pairs",random_pairs_interface,METH_VARARGS,pairs_docstring},
  {"choice_with_probs", choice_probs_interface, METH_VARARGS, choice_prob_docstring},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_pydigree(void)
{
  (void) Py_InitModule("pydigree._pydigree", module_methods);
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
  Py_ssize_t length = PyList_Size(chromA);
  PyObject* nc = recombine_haldane(chromA,chromB,map,length);
  return nc;

}

PyObject* recombine_haldane(PyObject* A, PyObject* B, PyObject* map, Py_ssize_t l) {
  Py_ssize_t i;
  unsigned char flipped = 0;
  double recombination_site = rexp(0.01);
  PyObject* newchrom = PyList_New(l);
  for (i = 0; i < l; i++) {
    /* Get the alleles */
    PyObject *a, *b;
    a = PyList_GetItem(A,i);
    b = PyList_GetItem(B,i);

    /* Do we cross over? */ 
    double position = PyFloat_AsDouble( PyList_GetItem(map,i) );
    if (recombination_site < position) {
      flipped = !flipped;
      recombination_site = rexp(0.01);
    }
    
    /* Put it on the new chromatid */
    PyObject* allele = flipped ? b : a;
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

