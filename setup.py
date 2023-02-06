"""
Setup script for installing vienna module
"""

# !/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

with open("README.md", "r", encoding="UTF-8") as f:
    readme = f.read()

with open("requirements.txt", "r", encoding="UTF-8") as f:
    requirements = f.read().splitlines()

setup(
    name="vienna",
    version="0.6.0",
    description="a short python wrapper for vienna tools ",
    long_description=readme + "\n\n",
    long_description_content_type="text/markdown",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    url="https://github.com/jyesselm/vienna",
    packages=[
        "vienna",
    ],
    package_dir={"vienna": "vienna"},
    py_modules=["vienna/vienna"],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords="vienna",
    classifiers=[
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
    ],
    entry_points={"console_scripts": []},
)
