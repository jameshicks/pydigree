# setup.py

import os.path

#from setuptools import Extension
#from distutils.core import setup
from setuptools import setup
from distutils.extension import Extension
from distutils.log import warn, error

try:
    import numpy
except ImportError:
    error('ERROR: Installation requires numpy')
    exit(1)

try:
    from Cython.Distutils import build_ext, cythonize
    cmdclass = {'build_ext': build_ext}
    use_cy = True
except ImportError:
    use_cy = False
    cmdclass = dict()
    warn('warning: Cython not found. Using C sources.')

if use_cy:
    fastio_sources = ['pydigree/fastio/fastio.pyx']
else:
    fastio_sources = ['pydigree/fastio/fastio.c']

if not all(os.path.exists(x) for x in fastio_sources):
    error('ERROR: Sources not found! Giving up.')
    exit(1)

c_ext = Extension(
    "pydigree._pydigree",
    ["pydigree/_pydigree.c"],
    include_dirs=[numpy.get_include()]
    )

fastio_ext = Extension('pydigree.io.fastio', fastio_sources)

setup(
    packages=['pydigree'],
    ext_modules=[c_ext, fastio_ext],
    requires=['numpy', 'scipy', 'bitarray']
    )
