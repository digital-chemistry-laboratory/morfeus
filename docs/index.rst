.. steriplus documentation master file, created by
   sphinx-quickstart on Sat Mar 23 22:40:04 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========
Steriplus
=========

.. figure:: images/logo.svg
  :align: left

Steriplus calculates molecular descriptors with a focus on sterics. It can be
used as a Python library or with console scripts.

********
Examples
********

************
Installation
************

.. code-block:: console
  :caption: pip

  pip install steriplus[extras]

.. code-block:: console
  :caption: conda

  conda install -c conda-forge steriplus

********
Features
********

* Buried volume
* Dispersion descriptor
* Ligand cone angles
* Local force constants
* Solvent accessible surface area
* Sterimol parameters

.. toctree::
   :maxdepth: 2
   :hidden:

   installation
   notes
   buried_volume
   cone_angle
   dispersion
   local_force
   sasa
   sterimol
   api
