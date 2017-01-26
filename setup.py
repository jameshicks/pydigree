"Package management script"

import sys
import os
import os.path

if sys.version_info[0:2] < (3, 3):
    print("Pydigree requires Python 3.2 or higher", file=sys.stderr)
    exit(1)

from setuptools import setup, find_packages

from distutils.extension import Extension
from distutils.log import error

from Cython.Build import cythonize
import numpy

cyprofile = bool(os.getenv("PYDIGREEPROFILING", False))

if cyprofile:
    from Cython.Compiler.Options import directive_defaults
    directive_defaults['linetrace'] = True
    directive_defaults['binding'] = True
    cython_macros = [('CYTHON_TRACE', '1')]
else:
    cython_macros = None

pydigree_dir = os.path.dirname(__file__)
cydigree_dir = os.path.join(pydigree_dir, 'pydigree', 'cydigree')
cysources = ['pydigree/cydigree/cyfuncs.pyx',
             'pydigree/cydigree/datastructures.pyx',
             'pydigree/cydigree/vcfparse.pyx']

cysources = [os.path.join(pydigree_dir, x) for x in cysources]

if not all(os.path.exists(x) for x in cysources):
    error('ERROR: Cython sources not found! Giving up.')
    exit(1)

cyext = Extension('pydigree.cydigree.cyfuncs',
                  sources=[os.path.join(cydigree_dir, 'cyfuncs.pyx')],
                  include_dirs=[numpy.get_include()],
                  extra_compile_args=['-Wno-unused-function'],
                  define_macros=cython_macros)

dsext = Extension('pydigree.cydigree.datastructures',
                  sources=[os.path.join(cydigree_dir, 'datastructures.pyx')],
                  extra_compile_args=['-Wno-unused-function'],
                  define_macros=cython_macros)

vcfext = Extension('pydigree.cydigree.vcfparse',
                   sources=[os.path.join(cydigree_dir, 'vcfparse.pyx')],
                   extra_compile_args=['-Wno-unused-function'],
                   define_macros=cython_macros)

sparray2 = Extension('pydigree.cydigree.sparsearray',
                     language='c++',
                     sources=[os.path.join(cydigree_dir, 'sparsearray.pyx')],
                     extra_compile_args=['-Wno-unused-function'],
                     define_macros=cython_macros)


with open(os.path.join(pydigree_dir, 'classifiers.txt')) as f:
    lib_class = [x.strip() for x in f if x.strip()]

lib_requirements = ['numpy', 'scipy', 'pandas', 'cython']

setup(
    name='pydigree',
    version='0.8a',
    description='A package for operations on pedigree and genotype data',
    author='James Hicks',
    author_email='jhicks22@wustl.edu',
    url='https://github.com/jameshicks/pydigree',
    download_url="https://github.com/jameshicks/pydigree/tarball/0.8a",
    license="MIT",
    packages=find_packages(),
    ext_modules=cythonize([cyext, dsext, vcfext, sparray2]),
    requires=lib_requirements,
    classifiers=lib_class)
