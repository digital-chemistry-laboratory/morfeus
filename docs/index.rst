.. raw:: html

  <picture>
      <source srcset="_static/logo-dark.svg" media="(prefers-color-scheme: dark)">
      <img src="_static/logo-light.svg?text=Light+mode" /> 
  </picture>

ᴍᴏʀғᴇᴜs calculates molecular features from 3D structures with a focus on steric
descriptors. It can be used as a Python library or through command line
scripts.

********
Examples
********

ᴍᴏʀғᴇᴜs can be imported as a Python module that is easily integrated into
workflows. Here is an example for calculating the exact ligand cone angle:

.. tab:: Module
  
  .. code-block:: python

    >>> from morfeus import ConeAngle, read_xyz
    >>> elements, coordinates = read_xyz("PdPMe3.xyz")
    >>> cone_angle = ConeAngle(elements, coordinates, 1)
    >>> print(cone_angle.cone_angle)
    117.11012922937584 

.. tab:: Command line

  .. code-block:: console

    $ morfeus cone_angle PdPMe3.xyz - 1 - cone_angle
    117.11012922937584    

************
Installation
************

Clone the repository from GitHub and install with pip:

.. code-block:: console
  :caption: pip installation

  pip install .

********
Features
********

* Buried volume
* Conformer tools
* Dispersion descriptor
* Exact ligand cone angle
* Local force constant
* Pyramidalization
* Solvent accessible surface area
* Sterimol parameters
* XTB electronic descriptors

*****
About
*****

ᴍᴏʀғᴇᴜs was developed by Kjell Jorner working as a postdoctoral fellow at
AstraZeneca UK.

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
   cli

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

.. toctree::
  :maxdepth: 3
  :caption: Source
  :hidden:

  GitHub <https://github.com/kjelljorner/morfeus>

.. _MIT license: https://en.wikipedia.org/wiki/MIT_License
