# setup.py
import numpy
from setuptools import setup, Extension

c_ext = Extension(
    "pydigree._pydigree",
    ["pydigree/_pydigree.c"],
    include_dirs=[numpy.get_include()]
    )

setup(
    packages=['pydigree'],
    ext_modules=[c_ext],
    requires=['numpy', 'scipy']
    )
