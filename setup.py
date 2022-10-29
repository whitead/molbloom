from setuptools import Extension, setup
from Cython.Build import cythonize

exec(open("fastz/version.py").read())

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name="fastz",
    ext_modules=cythonize(
        [Extension("fastz.bloom", ["fastz/bloom.pyx"],
                   libraries=["m"],
                   include_dirs=["fastz/bloom"])]),
    version=__version__,
    description="Purchaseable SMILES filter",
    author="Andrew White",
    author_email="andrew.white@rochester.edu",
    url="https://whitead.github.io/fastz/",
    license="MIT",
    packages=["fastz"],
    package_data={"fastz": ["filters/*.bloom"]},
    install_requires=[
        "importlib-resources",
    ],
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
