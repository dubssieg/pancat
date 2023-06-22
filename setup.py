#!/usr/bin/env python3
from sys import version_info, stderr
from setuptools import setup, find_packages

NAME = 'pangraphs'
CURRENT_PYTHON = version_info[:2]
REQUIRED_PYTHON = (3, 10)

if CURRENT_PYTHON < REQUIRED_PYTHON:
    stderr.write(
        f"{NAME} requires Python 3.10 or higher and your current version is {CURRENT_PYTHON}.")
    exit(1)

setup(
    name=NAME,
    version='0.1.0',
    description='Tools for manipulating and visualising GFA file format',
    url='https://github.com/Tharos-ux/gfatypes',
    author='Tharos',
    author_email='dubois.siegfried@gmail.com',
    packages=find_packages(),
    package_data={'': ['template.html']},
    include_package_data=True,
    zip_safe=False,
    license="LICENSE",
    long_description=open("README.md", encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    install_requires=['matplotlib', 'seaborn', 'numpy', 'pyvis', 'networkx', 'mycolorpy', 'levenshtein',
                      'gfagraphs', 'tharos-pytools', 'BubbleGun', 'statsmodels', 'resource', 'Bio', 'rich'],
    entry_points={'console_scripts': ['pangraphs=workspace.main:main']}
)
