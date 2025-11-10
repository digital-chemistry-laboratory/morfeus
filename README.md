![Logo Light](docs/_static/logo-text-light.svg#gh-light-mode-only)
![Logo Dark](docs/_static/logo-text-dark.svg#gh-dark-mode-only)

![PyPI - License](https://img.shields.io/pypi/l/morfeus-ml)
[![PyPI](https://img.shields.io/pypi/v/morfeus-ml)](https://pypi.org/project/morfeus-ml/)
[![Conda (channel only)](https://img.shields.io/conda/vn/conda-forge/morfeus-ml)](https://anaconda.org/conda-forge/morfeus-ml)
![Conda](https://img.shields.io/conda/pn/conda-forge/morfeus-ml)
![Python requires](https://img.shields.io/badge/dynamic/json?query=info.requires_python&label=python&url=https%3A%2F%2Fpypi.org%2Fpypi%2Fmorfeus-ml%2Fjson)
[![Testing](https://github.com/digital-chemistry-laboratory/morfeus/actions/workflows/test.yml/badge.svg)](https://github.com/digital-chemistry-laboratory/morfeus/actions/workflows/test.yml)
[![DOI](https://zenodo.org/badge/291745112.svg)](https://zenodo.org/badge/latestdoi/291745112)

A Python package for calculating molecular features.

# Installation

## pip

```shell
$ pip install morfeus-ml
```

## conda

```shell
$ conda install -c conda-forge morfeus-ml
```

# Usage

ᴍᴏʀғᴇᴜs can be imported as a Python module that is easily integrated into
workflows. Here is an example for calculating the exact ligand cone angle.

```shell
>>> from morfeus import ConeAngle, read_xyz
>>> elements, coordinates = read_xyz("PdPMe3.xyz")
>>> cone_angle = ConeAngle(elements, coordinates, 1)
>>> print(cone_angle.cone_angle)
117.11012922937587
```

It can also be used from the command line.

```console
$ morfeus cone_angle PdPMe3.xyz - 1 - cone_angle
117.11012922937587
```
For further information, see the separate [documentation](https://digital-chemistry-laboratory.github.io/morfeus/).

# Features

* Bite angle
* Buried volume
* Conformer tools
* Dispersion descriptor
* Exact ligand cone angle
* Ligand solid angle
* Local force constant
* Pyramidalization
* Solvent accessible surface area
* Sterimol parameters
* XTB electronic descriptors

# Acknowledgements

ᴍᴏʀғᴇᴜs was started by Kjell Jorner as a post doc at AstraZeneca UK in
collaboration with the groups of Alán Aspuru-Guzik at the University of
Toronto, Matthew Sigman at the University of Utah and Tobias Gensch at TU
Berlin. The package was further developed and maintained in the group of Kjell Jorner at ETH Zurich.

In particular, the following people (in alphabetical order) have contributed significantly to
developing its functionality:

* Gabriel dos Passos Gomes
* Kjell Jorner
* Lauriane Jacot-Descombes
* Pascal Friedrich
* Tobias Gensch
