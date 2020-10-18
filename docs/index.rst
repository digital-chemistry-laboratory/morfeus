.. figure:: images/logo.svg

ᴍᴏʀғᴇᴜs calculates molecular features from 3D structures with a focus on 
steric descriptors. It can be used as a Python library or through console
scripts.

********
Examples
********

ᴍᴏʀғᴇᴜs can be imported as a Python module that is easily integrated into
workflows. Here is an example for calculating the exact ligand cone angle:

.. code-block:: python
  :caption: Cone angle calculation

  >>> from morfeus import ConeAngle, read_xyz
  >>> elements, coordinates = read_xyz("phosphines/PdPMe3.xyz")
  >>> cone_angle = ConeAngle(elements, coordinates, 1)
  >>> print(cone_angle.cone_angle)
  117.11012922937584  

The metal atom here has index 1. Command-line scripts are also available for
convenience. Here is an example for a Sterimol calculation:

.. code-block:: console
  :caption: Sterimol calculation
  
  $ morfeus_sterimol tBu.xyz 1 2
  L         B_1       B_5
  4.21      2.87      3.27

The dummy atom here has index 1, while the connected atom in the substituent
has index 2.

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
* XTB electronic descriptors

*****
About
*****

ᴍᴏʀғᴇᴜs is developed by Kjell Jorner working as a Post Doc at AstraZeneca
in collaboration with the Aspuru-Guzik group at the University of Toronto.

*******
License
*******

ᴍᴏʀғᴇᴜs is released under the `MIT license`_ and is thus completely free for 
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
   xtb

.. toctree::
   :maxdepth: 3
   :caption: API
   :hidden:   
   
   api/morfeus

..
  api

.. _MIT license: https://en.wikipedia.org/wiki/MIT_License
