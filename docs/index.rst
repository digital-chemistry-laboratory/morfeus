.. steriplus documentation master file, created by
   sphinx-quickstart on Sat Mar 23 22:40:04 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========
Steriplus
=========

.. figure:: images/logo.svg

Steriplus calculates molecular descriptors from 3D structures with a focus on 
sterics. It can be used as a Python library or through console scripts.

********
Examples
********

Steriplus can be imported as a Python module that is easily integrated into
workflows.

.. code-block:: python
  :caption: Cone angle calculation

  >>> from steriplus import ConeAngle, read_xyz
  >>> elements, coordinates = read_xyz("phosphines/PdPMe3.xyz")
  >>> cone_angle = ConeAngle(elements, coordinates, 1)
  >>> print(cone_angle.cone_angle)
  117.11012922937584  

Command-line scripts are also available for convenience.

.. code-block:: console
  :caption: Sterimol calculation
  
  $ steriplus_sterimol tBu.xyz 1 2
  L         B_1       B_5
  4.21      2.87      3.27

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

*****
About
*****

Steriplus was developed by Kjell Jorner working as a Post Doc at AstraZeneca
in collaboration with the Aspuru-Guzik group at the University of Toronto.

*******
License
*******

Steriplus is released under the `MIT license`_ and is thus completely free for 
both academic and commercial use.

.. toctree::
   :maxdepth: 3
   :caption: Getting started
   :hidden:

   installation
   notes

.. toctree::
   :maxdepth: 3
   :caption: Methods
   :hidden:

   buried_volume
   cone_angle
   dispersion
   local_force
   sasa
   sterimol
   api

.. _MIT license: https://en.wikipedia.org/wiki/MIT_License