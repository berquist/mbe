from setuptools import setup


with open('README.md') as fh:
    long_description = fh.read()


setup(
    name='MBE',
    version='0.1.0',
    description='',
    long_description=long_description,
    author='Eric Berquist',
    author_email='erb74@pitt.edu',
    url='github or bitbucket',
    packages=['mbe'],
    install_requires=['six'],
    classifiers=[
    ],
)
