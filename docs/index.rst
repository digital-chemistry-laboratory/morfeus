.. figure:: images/logo.svg

Morfeus calculates molecular features from 3D structures with a focus on 
steric descriptors. It can be used as a Python library or through console
scripts.

********
Examples
********

Morfeus can be imported as a Python module that is easily integrated into
workflows.

.. code-block:: python
  :caption: Cone angle calculation

  >>> from morfeus import ConeAngle, read_xyz
  >>> elements, coordinates = read_xyz("phosphines/PdPMe3.xyz")
  >>> cone_angle = ConeAngle(elements, coordinates, 1)
  >>> print(cone_angle.cone_angle)
  117.11012922937584  

Command-line scripts are also available for convenience.

.. code-block:: console
  :caption: Sterimol calculation
  
  $ morfeus_sterimol tBu.xyz 1 2
  L         B_1       B_5
  4.21      2.87      3.27

************
Installation
************

Clone the repository from GitHub and install with pip:

.. code-block:: console
  :caption: pip

  pip install .

To install with extra graphics functionality:

.. code-block:: console
  :caption: pip

  pip install .[extras]

********
Features
********

* Buried volume
* Dispersion descriptor
* Ligand cone angle
* Local force constant
* Pyramidalization
* Solvent accessible surface area
* Sterimol parameters

*****
About
*****

Morfeus is developed by Kjell Jorner working as a Post Doc at AstraZeneca
in collaboration with the Aspuru-Guzik group at the University of Toronto.

*******
License
*******

Morfeus is released under the `MIT license`_ and is thus completely free for 
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
   pyramidalization
   sasa
   sterimol
   api

.. _MIT license: https://en.wikipedia.org/wiki/MIT_License