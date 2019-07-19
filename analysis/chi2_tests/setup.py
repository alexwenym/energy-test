
from distutils.core import setup
from Cython.Build import cythonize

setup(name='distance list', ext_modules=cythonize('distances.pyx'))
