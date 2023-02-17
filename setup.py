#!/usr/bin/env python3
from sys import version_info, stderr
from setuptools import setup, find_packages


CURRENT_PYTHON = version_info[:2]
REQUIRED_PYTHON = (3, 10)

if CURRENT_PYTHON < REQUIRED_PYTHON:
    stderr.write(
        f"Pangraphs requires Python 3.10 or higher and your current version is {CURRENT_PYTHON}.")
    exit(1)

with open("requirements.txt", "r", encoding='utf-8') as fh:
    requirements = [line.strip() for line in fh]

setup(
    name='pangraphs',
    version='0.0.1',
    description='Tools for manipulating and visualising GFA file format',
    url='https://github.com/Tharos-ux/gfatypes',
    author='Tharos',
    author_email='dubois.siegfried@gmail.com',
    packages=find_packages(),
    zip_safe=False,
    license="LICENSE",
    long_description=open("README.md", encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    install_requires=requirements,
    entry_points={'console_scripts': ['pangraphs=scripts.main:main']}
)
