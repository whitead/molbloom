from gettext import install
from setuptools import Extension, setup
from Cython.Build import cythonize
from ctypes.util import find_library
import os

exec(open("molbloom/version.py").read())

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name="molbloom",
    ext_modules=cythonize(
        [
            Extension(
                "molbloom.bloom",
                ["molbloom/bloom.pyx"],
                libraries=[find_library("m") if os.name == "nt" else "m"],
                include_dirs=["molbloom/bloom"],
            )
        ]
    ),
    include_dirs=["molbloom/bloom"],
    version=__version__,
    description="Purchaseable SMILES filter",
    author="Andrew White",
    author_email="andrew.white@rochester.edu",
    url="https://whitead.github.io/molbloom/",
    license="MIT",
    package_data={"molbloom": ["data/*.bloom"]},
    install_requires=["importlib_resources"],
    setup_requires=["cython"],
    packages=["molbloom"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Typing :: Typed",
    ],
)
