from setuptools import Extension, setup
from Cython.Build import cythonize

setup(
    name="fastz",
    ext_modules=cythonize(
        [Extension("fastz.bloom", ["fastz/bloom.pyx"],
                   libraries=["m"],
                   include_dirs=["fastz/bloom"])])
)
