from setuptools import find_packages, setup


with open("README.md") as fh:
    long_description = fh.read()

setup(
    name="mbe",
    version="0.1.0",
    description="Tools for generating many-body expansion expressions",
    long_description=long_description,
    author="Eric Berquist",
    author_email="eric.berquist@gmail.com",
    url="https://github.com/berquist/mbe",
    packages=find_packages(),
    install_requires=["numpy", "sympy", "toolz"],
)
