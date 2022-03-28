==============
Notes on usage
==============

************
Atom indices
************

The ᴍᴏʀғᴇᴜs API uses atom indices starting from 1, so called one-based
indexing. This means that the first atom in a molecule has index 1. In
contrast, Python generally uses zero-based indexing for lists and arrays.

.. code-block:: none
  :caption: One-based indexing from xyz file
  :emphasize-lines: 3

  14
  PdPMe3
  Pd         -0.743397 -0.000992 -1.518770                # Atom with index 1
  P           1.005896  0.002102 -0.163618                # Atom with index 2
  ...

We believe that one-based indexing is more intuitive for the user when, *e.g.*,
getting indices visually from molecular viewers. To accomodate for this
difference, atom properties in morfeus are given as dictionaries, with keys
corresponding to the atom indices. These can easily be turned into regular
lists and used with zero-indexing if needed:

.. code-block:: python
  :caption: Accessing atom properties

  >>> from morfeus import SASA, read_geometry
  >>> elements, coordinates = read_geometry("Et.gjf")
  >>> sasa = SASA(elements, coordinates)
  >>> sasa.atom_areas
  {1: 22.489947408719903,
  2: 19.180448148100027,
  3: 22.529947315182724,
  ...
  >>> list(sasa.atom_areas.values())
  [22.489947408719903,
  19.180448148100027,
  22.529947315182724,
  ...

**************
Geometry files
**************

Many features of ᴍᴏʀғᴇᴜs makes use of atomic numbers/symbols and coordinates.
These can be read from geometry files, and currently the files formats ``gjf``
(Gaussian input file) and ``xyz`` (XMOL )are supported. These files can be
read with two different functions. There is also the
:py:func:`read_geometry <morfeus.io.read_geometry>`
function that will try to guess the file type based on its suffix. If the
cclib__ package is installed, it can be used to read many file formats.

====== =============================================
Format Function
====== =============================================
gjf    :py:func:`read_gjf <morfeus.io.read_gjf>`
xyz    :py:func:`read_xyz <morfeus.io.read_gjf>`
\*     :py:func:`read_cclib <morfeus.io.read_cclib>`
====== =============================================

.. __: https://github.com/cclib/cclib

*****
Radii
*****

ᴍᴏʀғᴇᴜs makes extensive use of atomic radii in calculating the different
descriptors. There are mainly two types of radii: vdW and covalent.

#######
Options
#######

=============================================== =================================================
vdW                                             Covalent
=============================================== =================================================
:py:data:`alvarez <morfeus.data.radii_alvarez>` :py:data:`pyykko <morfeus.data.cov_radii_pyykko>`
:py:data:`bondi <morfeus.data.radii_bondi>`
:py:data:`crc <morfeus.data.radii_crc>`
:py:data:`rahm <morfeus.data.radii_rahm>`
:py:data:`truhlar <morfeus.data.radii_truhlar>`
=============================================== =================================================

- Default vdW radii are from the curated collection in the CRC handbook
- The only option for covalent radii at the moment is those from pyykko
- Rahm radii are a special type derived from electron density calculation. They
  are the default for dispersion descriptor calculations.

Changing the type of radii is done by specifying the ``--radii_type`` keyword
argument.

#######
Sources
#######

The data is extracted from the excellent mendeleev__ package.

.. __: https://github.com/lmmentel/mendeleev/

***************
Without display
***************

Matplotlib may cause problems on computers without displays, e.g., nodes on
computer clusters. This can be solved by changing the plotting backend to
"agg". This could be done at the top of the Python script using ᴍᴏʀғᴇᴜs:

.. code-block:: python

  import matplotlib
  matplotlib.use('Agg')

Alternatively, set the environment variable in the shell before launching the
script (Linux example):

.. code-block:: console

  export MPLBACKEND="agg"
