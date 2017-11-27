#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup
except:
    from distutils.core import setup

long_description = None
with open('README.rst') as file:
    long_description = file.read()

setup(
    name = "pymultinest",
    version = "2.6",
    description = "Access modules for the MultiNest and Cuba integration libraries",
    author = "Johannes Buchner",
    author_email = "johannes.buchner.acad@gmx.com",
    maintainer = "Johannes Buchner",
    maintainer_email = "johannes.buchner.acad@gmx.com",
    url = "http://johannesbuchner.github.com/PyMultiNest/",
    license = "GPLv3",
    packages = ["pymultinest", "pycuba", "pypolychord"],
    provides = ["pymultinest", "pycuba", "pypolychord"],
    requires = ["ctypesGsl", "numpy (>=1.5)", "matplotlib", "scipy"],
    scripts=['multinest_marginals.py'],
    long_description=long_description,
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)

