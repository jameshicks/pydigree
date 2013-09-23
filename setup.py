# setup.py

from setuptools import setup,Extension

c_ext = Extension("pydigree._pydigree", ["pydigree/_pydigree.c"])

setup(
    packages=['pydigree'],
    ext_modules=[c_ext],
    requires=['numpy']
    )
    
