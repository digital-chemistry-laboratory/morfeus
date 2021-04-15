=========
Morfeus
=========
.. image:: docs/images/logo.svg
.. image:: https://github.com/kjelljorner/morfeus/actions/workflows/ci-pip.yml/badge.svg

A Python package for calculation of molecular features. Morfeus can be run
either as a command line script or imported as a module.

************
Installation
************

Clone the repository and install with

.. code-block:: console

  pip install .

Installation with extra graphics functionality

.. code-block:: console

  pip install .[extras]

*************
Documentation
*************

To be added later.

************
Dependencies
************

Core dependencies:

* numpy
* scipy

Optional depedencies when installed with [extras]:

* matplotlib
* pyvista
* vtk

The optional dependencies are used for 3D visualization and in the Dispersion
descriptor calculations.