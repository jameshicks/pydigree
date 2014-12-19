# setup.py

import os.path

#from setuptools import Extension
#from distutils.core import setup
from setuptools import setup
from distutils.extension import Extension
from distutils.log import warn, error


from Cython.Build import cythonize
from Cython.Distutils import build_ext

import numpy 



c_ext = Extension(
    "pydigree._pydigree",
    ["pydigree/_pydigree.c"],
    include_dirs=[numpy.get_include()]
    )

cysources = ['pydigree/cyfuncs.pyx']
if not all(os.path.exists(x) for x in cysources):
    error('ERROR: Sources not found! Giving up.')
    exit(1)


cyext = [Extension('pydigree.cyfuncs', cysources, include_dirs=[numpy.get_include()])]

setup(
    packages=['pydigree'],
    ext_modules=[c_ext] + cythonize(cyext),
    requires=['numpy', 'scipy']
    )
