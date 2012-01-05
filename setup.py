#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup
except:
    from distutils.core import setup

setup(
    name = "pymultinest",
    version = "0.1",
    description = "Access modules for MultiNest and APEMoST",
    author = "Johannes Buchner",
    author_email = "johannes.buchner.acad [ät] gmx.com",
    url = "http://johannesbuchner.github.com/PyMultiNest/",
    license = "GPLv3",
    packages = ["pymultinest", "pyapemost"]
)
