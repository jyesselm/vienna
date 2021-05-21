#!/usr/bin/env python

import os
import sys
import shutil

def does_program_exist(prog_name):
    if shutil.which(prog_name) is None:
        return False
    else:
        return True

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if not does_program_exist("RNAfold"):
    raise ValueError("cannot install without rnafold being present")

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.rst').read()
doclink = """
Documentation
-------------
"""
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='vienna',
    version='0.1.0',
    description='a short python wrapper for vienna tools ',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='Joe Yesselman',
    author_email='jyesselm@unl.edu',
    url='https://github.com/jyesselm/vienna',
    packages=[
        'vienna',
    ],
    package_dir={'vienna': 'vienna'},
    py_modules=[
        'vienna/vienna'
    ],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='vienna',
    classifiers=[
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
    entry_points = {
        'console_scripts' : [
        ]
    }
)
