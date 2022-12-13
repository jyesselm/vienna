# vienna

[![PYPI package](https://badge.fury.io/py/vienna.png)](http://badge.fury.io/py/vienna)
[![Build status](https://travis-ci.org/jyesselm/vienna.png?branch=main)](https://travis-ci.org/jyesselm/vienna)
[![linting: pylint](https://img.shields.io/badge/linting-pylint-yellowgreen)](https://github.com/PyCQA/pylint)
[![formatting: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

a short python wrapper for vienna
tools (https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/install.html) Note I did not develop this
software just built this short wrapper.

## how to install

### install anaconda

highly recommended to install anaconda first and create a py3 environment

```shell
# to install conda on linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh

conda create --name py3 python=3.7 --yes
conda activate py3

# to install conda on mac
# https://docs.anaconda.com/anaconda/install/mac-os/
bash ~/Downloads/Anaconda3-2020.02-MacOSX-x86_64.sh
```

### install vienna

Note this is just a wrapper so you must install the vienna fold code

```shell
# this can be accomplished using conda 
conda install -c bioconda viennarna

# if you do not wish to use conda 
# you can install it via many options from the website 
https://www.tbi.univie.ac.at/RNA/
```

### install vienna python package

```shell
# vienna python package is avialible on PYPI 
pip install vienna

# you can also install it from github 
git clone https://github.com/jyesselm/vienna
cd vienna
python setup.py install 
```

## how to use

```shell
>>> import vienna 
>>> fr = vienna.fold('GGGGAAAACCCC')
>>> print(fr)
FoldResults(dot_bracket='((((....))))', mfe=-5.4, ens_defect=1.62)
```