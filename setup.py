try:
    from setuptools import setup
except:
    from distutils.core import setup

setup(
    name = "pymultinest",
    version = "0.1",
    packages = ["pymultinest", "pyapemost"]
)
