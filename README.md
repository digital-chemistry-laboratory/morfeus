![Logo](docs/_static/logo-light.svg)
![PyPI - License](https://img.shields.io/pypi/l/morfeus-ml)
[![PyPI](https://img.shields.io/pypi/v/morfeus-ml)](https://pypi.org/project/morfeus-ml/)
![Python requires](https://img.shields.io/badge/dynamic/json?query=info.requires_python&label=python&url=https%3A%2F%2Fpypi.org%2Fpypi%2Fmorfeus-ml%2Fjson)
[![Testing](https://github.com/kjelljorner/morfeus/actions/workflows/test.yml/badge.svg)](https://github.com/kjelljorner/morfeus/actions/workflows/test.yml)

A Python package for calculating molecular features.

# Installation

## pip

```console
$ pip install morfeus-ml
```

## conda

```console
$ conda install -c conda-forge morfeus-ml
```

# Usage

ᴍᴏʀғᴇᴜs can be imported as a Python module that is easily integrated into
workflows. Here is an example for calculating the exact ligand cone angle.

```python
>>> from morfeus import ConeAngle, read_xyz
>>> elements, coordinates = read_xyz("PdPMe3.xyz")
>>> cone_angle = ConeAngle(elements, coordinates, 1)
>>> print(cone_angle.cone_angle)
117.11012922937584 
```

It can also be used from the command line.

```console
$ morfeus cone_angle PdPMe3.xyz - 1 - cone_angle
117.11012922937584   
```
For further information, see the separate [documentation](https://kjelljorner.github.io/morfeus/).

# Features

* Buried volume
* Conformer tools
* Dispersion descriptor
* Exact ligand cone angle
* Local force constant
* Pyramidalization
* Solvent accessible surface area
* Sterimol parameters
* XTB electronic descriptors

# Acknowledgements

ᴍᴏʀғᴇᴜs was started by Kjell Jorner as a post doc at AstraZeneca UK in
collaboration with the groups of Alán Aspuru-Guzik at the University of
Toronto, Matthew Sigman at the University of Utah and Tobias Gensch at TU
Berlin. In particular, the following people have contributed significantly to
developing its functionality:

* Gabriel dos Passos Gomes
* Pascal Friedrich
* Tobias Gensch
