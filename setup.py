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
    version = "1.4",
    description = "Access modules for MultiNest and APEMoST",
    author = "Johannes Buchner",
    author_email = "johannes.buchner.acad [@t] gmx.com",
    maintainer = "Johannes Buchner",
    maintainer_email = "johannes.buchner.acad [@t] gmx.com",
    url = "http://johannesbuchner.github.com/PyMultiNest/",
    license = "GPLv3",
    packages = ["pymultinest", "pyapemost", "pycuba"],
    provides = ["pymultinest", "pyapemost", "pycuba"],
    requires = ["ctypesGsl", "numpy (>=1.5)", "matplotlib", "scipy"],
    scripts=['multinest_marginals.py'],
    long_description=long_description,
)

