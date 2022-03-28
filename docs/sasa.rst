===============================
Solvent accessible surface area
===============================
Solvent accessible surface area (SASA), atomic SASA and volumes under the
SASA can be calculated.

******
Module
******

The SASA class calculates and stores the total and atomic SASA as well as the
volume.

.. code-block:: python
  :caption: Example

  >>> from morfeus import SASA, read_xyz
  >>> elements, coordinates = read_xyz("n-heptane.xyz")
  >>> sasa = SASA(elements, coordinates)
  >>> print(sasa.atom_areas[1])
  18.380429455791376
  >>> print(sasa.area)
  331.5607124071378
  >>> print(sasa.volume)
  475.5699458352845

The ``atom_areas`` dictionary contains the atomic SASAs indexed from 1. Type of
radii can be changed with the keyword argument ``radii=<str>``  and custom
radii can be supplied with ``radii=<list>``. The probe radius is changed with
``probe_radius=<float>``.

For more information, use ``help(SASA)`` or consult the API:
:py:class:`SASA <morfeus.morfeus.SASA>`

*******************
Command line script
*******************

The command-line script outputs total SASA and volume as well as SASA per atom.

.. code-block:: console
  :caption: Example

  $ morfeus sasa PdPMe3.xyz - - print_report
  Probe radius (Å): 1.4
  Solvent accessible surface area (Å²): 288.3
  Volume inside solvent accessible surface (Å³): 410.7
  $ morfeus sasa PdPMe3.xyz - - print_report --verbose=True
  Probe radius (Å): 1.4
  Solvent accessible surface area (Å²): 288.3
  Volume inside solvent accessible surface (Å³): 410.7
  Symbol    Index     Area (Å²)
  Pd        1         91.8
  P         2         0.0
  C         3         13.4
  H         4         18.2
  H         5         15.6
  H         6         18.2
  C         7         13.5
  H         8         18.2
  H         9         15.6
  H         10        18.2
  C         11        13.5
  H         12        18.2
  H         13        18.2
  H         14        15.6

**********
Background
**********

Solvent accessible surface area is a measure of how much of the area of a
molecule is available to the solvent. The atomic SASA can be used as a measure
of the steric availability of an atom. ᴍᴏʀғᴇᴜs uses a modified version of the
method of Shrake and Rupley :footcite:`shrake_environment_1973` where a
constant surface density of points is used instead of a fixed number of points
regardless of the atom area. The atomic SASA and volumes are computed as
described by Eisenhaber *et al.* :footcite:`eisenhaber_double_1995`. ᴍᴏʀғᴇᴜs is
not optimized for larger molecules and other programs are recommended for,
*e.g.*, proteins.

.. footbibliography::
