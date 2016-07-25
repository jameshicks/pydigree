# setup.py

import os.path

from setuptools import setup

from distutils.extension import Extension
from distutils.log import warn, error

from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy

cyprofile = False

if cyprofile:
    from Cython.Compiler.Options import directive_defaults
    directive_defaults['linetrace'] = True
    directive_defaults['binding'] = True


cysources = ['pydigree/cydigree/cyfuncs.pyx']
if not all(os.path.exists(x) for x in cysources):
    error('ERROR: Cython sources not found! Giving up.')
    exit(1)

if cyprofile:
    macros = [('CYTHON_TRACE', '1')]
else:
    macros = None

cyext = [Extension('pydigree.cyfuncs',
                   sources=cysources,
                   include_dirs=[numpy.get_include()],
                   extra_compile_args=['-Wno-unused-function'],
                   define_macros=macros)]

with open('LICENSE.txt') as f:
    license = f.read()

setup(
        name='pydigree',
        description='A package for operations on pedigree and genotype data',
        author='James Hicks',
        url='https://github.com/jameshicks/pydigree',
        license=license,
        packages=['pydigree'],
        ext_modules=cythonize(cyext),
        requires=['numpy', 'scipy', 'pandas', 'cython']
    )
