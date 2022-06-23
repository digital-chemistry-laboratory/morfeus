###########
Solid angle
###########

Ligand solid angle and solid cone angles are implemented with a numerical recipe as
described by Guzei and Wendt :footcite:`guzei_2006`.

******
Module
******

The SolidAngle class is provided to calculate and store the solid angles.

.. code-block:: python
  :caption: Example

  >>> from morfeus import SolidAngle, read_xyz
  >>> elements, coordinates = read_xyz("phosphines/PdPMe3.xyz")
  >>> solid_angle = SolidAngle(elements, coordinates, 1)
  >>> print(solid_angle.cone_angle)
  112.53614024551955
  >>> print(solid_angle.solid_angle)
  2.794082404625142
  >>> print(solid_angle.G)
  22.234601305109027
  >>> solid_angle.print_report()
  Solid angle (sr): 2.794
  Cone angle (°): 112.536
  G: 22.235

The Bondi vdW radii are used in reference :footcite:`bilbrey_solid_2013`, but
radii from the CRC Handbook is the default here. It can be changed with
``radii_type=<str>``. Custom radii can passed with ``radii=<list>``. The units
for the solid cone angle are in degrees, steradians for the solid angle and the
G parameter is in percent.

For more detailed information, use ``help(SolidAngle)`` or see the API:
:py:class:`SolidAngle <morfeus.solid_angle.SolidAngle>`.

*******************
Command line script
*******************

The command line script provides access to the basic functionality through the
terminal.

.. code-block:: console
  :caption: Example

  $ morfeus solid_angle PdPMe3.xyz - 1 - print_report
  Solid angle (sr): 2.794
  Cone angle (°): 112.536
  G: 22.235

**********
Background
**********

The ligand solid angle, Ω is the solid angle of the complete shadow cast by a
ligand when hypothetically illuminated from the metal center.
:footcite:`bilbrey_solid_2013` From Ω, a solid "cone" angle Θ can also be
calculated, which is analogous to Tolman's cone angle. The G parameter is a
measure of how much of the metal coordination sphere is shielded by the
ligands. :footcite:`guzei_2006`

ᴍᴏʀғᴇᴜs uses a numerical method, in a similar way to the Solid-G__ program, but
with different radii. The results have been benchmarked against the exact solid
angle program by Allen and co-workers, :footcite:`bilbrey_solid_2013` and agree
within numerical accuracy.

.. __: https://xray.chem.wisc.edu/solid-g/

.. footbibliography::
