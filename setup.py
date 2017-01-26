"Package management script"

import sys
import os.path

if sys.version_info[0:2] < (3, 3):
    raise ImportError('pydigree requires python >3.2')

from setuptools import setup, find_packages

from distutils.extension import Extension
from distutils.log import error

from Cython.Build import cythonize
import numpy

cyprofile = False

if cyprofile:
    from Cython.Compiler.Options import directive_defaults
    directive_defaults['linetrace'] = True
    directive_defaults['binding'] = True

pydigree_dir = os.path.dirname(__file__)
cydigree_dir = os.path.join(pydigree_dir, 'pydigree', 'cydigree')
cysources = ['pydigree/cydigree/cyfuncs.pyx',
             'pydigree/cydigree/datastructures.pyx',
             'pydigree/cydigree/vcfparse.pyx']

cysources = [os.path.join(pydigree_dir, x) for x in cysources]

if not all(os.path.exists(x) for x in cysources):
    error('ERROR: Cython sources not found! Giving up.')
    exit(1)

if cyprofile:
    macros = [('CYTHON_TRACE', '1')]
else:
    macros = None

cyext = Extension('pydigree.cydigree.cyfuncs',
                  sources=[os.path.join(cydigree_dir, 'cyfuncs.pyx')],
                  include_dirs=[numpy.get_include()],
                  extra_compile_args=['-Wno-unused-function'],
                  define_macros=macros)

dsext = Extension('pydigree.cydigree.datastructures',
                  sources=[os.path.join(cydigree_dir, 'datastructures.pyx')],
                  extra_compile_args=['-Wno-unused-function'],
                  define_macros=macros)

vtext = Extension('pydigree.cydigree.varianttree',
                  sources=[os.path.join(cydigree_dir, 'varianttree.pyx')],
                  extra_compile_args=['-Wno-unused-function'],
                  define_macros=macros)

vcfext = Extension('pydigree.cydigree.vcfparse',
                   sources=[os.path.join(cydigree_dir, 'vcfparse.pyx')],
                   extra_compile_args=['-Wno-unused-function'],
                   define_macros=macros)

sparray2 = Extension('pydigree.cydigree.sparsearray',
                     language='c++',
                     sources=[os.path.join(cydigree_dir, 'sparsearray.pyx')],
                     extra_compile_args=['-Wno-unused-function'],
                     define_macros=macros)


with open(os.path.join(pydigree_dir, 'LICENSE.txt')) as f:
    liblicense = f.read()

lib_class = ['Programming Language :: Python :: 3 :: Only',
             'Programming Language :: Cython',
             'Intended Audience :: Science/Research',
             'Topic :: Scientific/Engineering :: Bio-Informatics',
             'Topic :: Scientific/Engineering :: Mathematics',
             'Topic :: Sociology :: Genealogy',
             'Topic :: Software Development :: Libraries :: Python Modules',
             'Development Status :: 4 - Beta',
             'License :: OSI Approved :: Apache Software License']
lib_requirements = ['numpy', 'scipy', 'pandas', 'cython']

setup(
    name='pydigree',
    version='0.8a',
    description='A package for operations on pedigree and genotype data',
    author='James Hicks',
    author_email='jhicks22@wustl.edu',
    url='https://github.com/jameshicks/pydigree',
    download_url="https://github.com/jameshicks/pydigree/tarball/0.8a",
    license=liblicense,
    packages=find_packages(),
    ext_modules=cythonize([cyext, dsext, vtext, vcfext, sparray2]),
    requires=lib_requirements,
    classifiers=lib_class)
