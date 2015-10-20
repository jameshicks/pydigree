# setup.py

import os.path

#from setuptools import Extension
#from distutils.core import setup
from setuptools import setup
from distutils.extension import Extension
from distutils.log import warn, error


from Cython.Build import cythonize
from Cython.Distutils import build_ext

cyprofile = False
if cyprofile:
    from Cython.Compiler.Options import directive_defaults
    directive_defaults['linetrace'] = True
    directive_defaults['binding'] = True

import numpy


c_ext = Extension(
    "pydigree._pydigree",
    ["pydigree/_pydigree.c"],
    include_dirs=[numpy.get_include()]
)

cysources = ['pydigree/cyfuncs.pyx']
if not all(os.path.exists(x) for x in cysources):
    error('ERROR: Cython sources not found! Giving up.')
    exit(1)

if cyprofile:
    macros = [('CYTHON_TRACE', '1')]
else:
    macros = None

cyext = [Extension('pydigree.cyfuncs', cysources,
                   include_dirs=[numpy.get_include()], define_macros=macros)]

setup(
    packages=['pydigree'],
    ext_modules=[c_ext] + cythonize(cyext),
    requires=['numpy', 'scipy', 'pandas']
)
