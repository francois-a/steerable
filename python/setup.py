from distutils.core import setup, Extension
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy

setup(ext_modules = cythonize(Extension(
    "steerablefilter",
    sources=["steerablefilter.pyx", "../src/steerableDetector.cpp"],
    include_dirs=[numpy.get_include(), "/usr/local/include/", '../src/'],
    language="c++",
    extra_link_args=["-lgsl"],
)))
