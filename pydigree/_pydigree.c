#include <Python.h>
#include <math.h>

#include <numpy/arrayobject.h>
/* Function Prototypes */ 



PyObject* choice_probs_interface(PyObject* self, PyObject* args);
PyObject* choice_probs(PyObject* choices, PyObject* probs);


static inline PyObject* sgs_shares(PyObject* affecteds, PyObject* shared, Py_ssize_t nmark); 
PyObject* sgs_shares_interface(PyObject* self, PyObject* args);

/* docstrings */
static char module_docstring[] =
  "This sub-module provides C implementations of time-intensive pydigree functions";


static char choice_prob_docstring[] = "Chooses a random item based on probabilities provided in probs"; 
static char sgs_shares_docstring[] = "Returns the proportion of individual pairs IBD";
/* Python C API boilerplate */
static PyMethodDef module_methods[] = {
  {"choice_with_probs", choice_probs_interface, METH_VARARGS, choice_prob_docstring},
  {"sgs_shares", sgs_shares_interface, METH_VARARGS, sgs_shares_docstring},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_pydigree(void)
{
  (void) Py_InitModule("pydigree._pydigree", module_methods);
  srand(time(NULL));
  import_array();
}


/* Common functions */ 
static inline int rand_lim(int limit) {
  /* return a random number between 0 (inclusive) and limit (exclusive).
   */
  int divisor = RAND_MAX/(limit);
  int retval;

  do { 
    retval = lrand48() / divisor;
  } while (retval > limit);

  return retval;
}

static inline int rand_int_range(int min, int max) {
  return rand_lim(max) + min; 
}

static inline int rand_bool(double prob) { 
  return drand48() < prob; 
}

static inline double rexp(double rate) { 
  return log(1-drand48()) / (-rate);
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


PyObject* sgs_shares_interface(PyObject* self, PyObject* args) {
  PyObject *shared, *affecteds;
  Py_ssize_t nmark;
  if (!PyArg_ParseTuple(args, "OOn", &affecteds, &shared, &nmark)) return NULL;
  return sgs_shares(affecteds, shared, nmark);

}

static inline PyObject* sgs_shares(PyObject* affecteds, PyObject* shared, Py_ssize_t nmark) {
  Py_ssize_t i,j,k,m, nsegs, start, stop; 
  PyObject *a, *b, *pair;
  PyObject *sharelocs, *loc;
  PyObject *nsharesarray;

  Py_ssize_t naff = PyList_Size(affecteds);
  PyObject *inds[naff];
  for (i=0; i < naff; i++) {
    inds[i] = PyList_GET_ITEM(affecteds, i);
  }

  /* Initialize empty array */ 
  unsigned long int shares[nmark];
  for (i = 0; i < nmark; i++) {
    shares[i] = 0;
  }

  /* Start incrementing */
  for (i=0; i < naff; i++) {
    a = inds[i];
    for (j=i+1; j < naff; j++) {
      b = inds[j];

      /* Get shares between individuals */
      pair = PyFrozenSet_New(NULL);
      PySet_Add(pair, a);
      PySet_Add(pair, b);
      sharelocs = PyDict_GetItem(shared, pair);
      if (sharelocs == NULL) {
	/* 
	   I'm expecting a KeyError here but in the case of anything
	   I'm just going to skip it anyway
	*/
	PyErr_Clear();
	goto finish;
      }

      /* For each share, increment the count inside it */
      nsegs = PyList_GET_SIZE(sharelocs);
      for (k=0; k < nsegs; k++) {
	loc = PyList_GET_ITEM(sharelocs, k);
	start = PyInt_AsSsize_t(PyTuple_GET_ITEM(loc, 0));
	stop = PyInt_AsSsize_t(PyTuple_GET_ITEM(loc, 1));
	for (m=start; m <= stop; m++) {
	  shares[m]++;
	}
      }
    finish:
      Py_XDECREF(pair);
    }
  }

  /* Make results into a numpy array */
  npy_intp dims[1] = {nmark};
  nsharesarray = PyArray_SimpleNew(1, &dims, NPY_UINT16);

  for (i=0; i < nmark; i++) {
    void *itemPtr = PyArray_GETPTR1(nsharesarray, i);
    PyArray_SETITEM(nsharesarray, itemPtr, PyLong_FromUnsignedLong(shares[i]));
  }


  return PyArray_Return(nsharesarray);
}
