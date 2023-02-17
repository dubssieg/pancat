#!/usr/bin/env python3
from sys import version_info, stderr
from setuptools import setup, find_packages


CURRENT_PYTHON = version_info[:2]
REQUIRED_PYTHON = (3, 10)

if CURRENT_PYTHON < REQUIRED_PYTHON:
    stderr.write(
        f"Pangraphs requires Python 3.10 or higher and your current version is {CURRENT_PYTHON}.")
    exit(1)


setup(
    name='pangraphs',
    version='0.0.1',
    description='Abstraction layer for GFA file format',
    url='https://github.com/Tharos-ux/gfatypes',
    author='Tharos',
    author_email='dubois.siegfried@gmail.com',
    packages=find_packages(),
    zip_safe=False,
    license="LICENSE",
    long_description=open("README.md", encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    install_requires=['networkx', 'tharos-pytools', 'pyvis',
                      'mycolorpy', 'levenshtein', 'gfatypes', 'gfagraphs', 'Bio'],
    entry_points={'console_scripts': ['pangraphs=scripts.pangraphs:main']}
)
