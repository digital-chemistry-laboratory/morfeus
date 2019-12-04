=========
Morfeus
=========

A Python package for calculation of molecular features. Morfeus can be run
either as a command line script or imported as a module.

************
Installation
************

###
pip
###

Installation of core functionality with pip

.. code-block:: console

  pip install morfeus

Installation with extra graphics functionality

.. code-block:: console

  pip install morfeus[extras]


#####
conda
#####

Conda currently only supports installation with extra graphics functionality.

.. code-block:: console

  conda install -c conda-forge morfeus

*************
Documentation
*************

https://morfeus.readthedocs.io/en/latest/index.html

************
Dependencies
************

Core dependencies:

* numpy
* scipy

Optional depedencies:

* matplotlib
* pyvista
* vtk

The optional dependencies are used for 3D visualization and in the Dispersion
descriptor calculations.