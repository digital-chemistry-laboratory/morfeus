##########
Cone angle
##########

Exact cone angles are implemented as described by Allen and co-workers
:footcite:`bilbrey_exact_2013`.

******
Module
******

The ConeAngle class is provided to calculate and store the cone angles.

.. code-block:: python
  :caption: Example

  >>> from morfeus import ConeAngle, read_xyz
  >>> elements, coordinates = read_xyz("phosphines/PdPMe3.xyz")
  >>> cone_angle = ConeAngle(elements, coordinates, 1)
  >>> print(cone_angle.cone_angle)
  117.11012922937584
  >>> print(cone_angle.tangent_atoms)
  [5, 9, 12]
  >>> cone_angle.print_report()
  Cone angle: 117.1
  No. tangent atoms: 3
  >>> cone_angle.plot_3D()

The Bondi vdW radii are used in reference :footcite:`bilbrey_exact_2013`, but
radii from the CRC Handbook is the default here. It can be changed with
``radii_type=<str>`` with either ``crc`` or ``bondi``. Custom radii can passed
with ``radii=<list>``.

The default setting is ``method="libconeangle"``, which uses the fast
libconeangle__ package as a backend. If it is not installed, an internal
algorithm will be used (``method="internal"``), printing a warning method.

For more detailed information, use ``help(ConeAngle)`` or see the API:
:py:class:`ConeAngle <morfeus.cone_angle.ConeAngle>`.

.. __: https://github.com/kjelljorner/libconeangle

*******************
Command line script
*******************

The command line script provides access to the basic functionality through the
terminal.

.. code-block:: console
  :caption: Example

  $ morfeus cone_angle PdPMe3.xyz - 1 - print_report
  Cone angle: 117.1
  No. tangent atoms: 3
  Tangent to: H6 H10 H13

**********
Background
**********

Cone angles is a method invented by Tolman for assessing the steric size of
ligands :footcite:`tolman_steric_1977`. The original Tolman cone angles for
phosphines have problems with asymmetric ligands and are not implemented in
this package. Instead, the exact cone angles :footcite:`bilbrey_exact_2013` are
used. These are also defined for multidentate ligands.

The method implemented in ᴍᴏʀғᴇᴜs is taken directly from the article by Allen
:footcite:`bilbrey_exact_2013`. The results have been benchmarked against the
original article and agree within numerical accuracy.

.. footbibliography::
