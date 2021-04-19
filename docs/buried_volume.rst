=============
Buried volume
=============

Buried volumes are implemented as described by Cavallo and co-workers
:footcite:`falivene_sambvca_2016`. 

******
Module
******

The BuriedVolume class calculates the buried volume. Steric maps can also be
plotted.

.. code-block:: python
  :caption: Example

  >>> from morfeus import BuriedVolume, read_xyz
  >>> elements, coordinates = read_xyz("1.xyz")
  >>> bv = BuriedVolume(elements, coordinates, 1, excluded_atoms=[1, 2, 3, 4, 5, 6, 7])
  >>> print(bv.buried_volume)
  0.2962110976518822
  >>> bv.print_report()
  V_bur (%): 29.6
  >>> bv.plot_steric_map([14])

.. image:: images/steric_map.png

Plots can be saved by passing the ``filename=<str>`` keyword argument.
:py:meth:`draw_3D <morfeus.buried_volume.BuriedVolume.draw_3D>` gives a
three-dimensional representation of the buried volume.


By default, hydrogen atoms are excluded in the calculation. They can be added
by giving the keyword argument ``indclude_hs=True``. The default sphere radius
is 3.5 Å, but can be changed with ``radius=<float>``. Default radii type is
Bondi which are scaled by a factor of 1.17. This can be changed with
``radii_type=<str>`` and ``radii_scale=<float>``. Custom radii can be supplied
as a list with ``radii=<list>``.

For more information, use ``help(BuriedVolume)`` or see the API:
:py:class:`BuriedVolume <morfeus.buried_volume.BuriedVolume>`

*******************
Command line script
*******************

The basic functionality is available through the command line script.

.. code-block:: console
  :caption: Example

  $ morfeus buried_volume 1.xyz - 1 --excluded_atoms='[1,2,3,4,5,6,7]' - print_report
  V_bur (%): 29.6

**********
Background
**********

The percent of buried volume is a measure of the steric hindrance induced by a
ligand of a transition metal complex :footcite:`falivene_sambvca_2016`. A web tool
to calculate buried volumes, SambVca, was made available for scientific
purposes by Cavallo and co-workers in 2009 :footcite:`poater_sambvca_2009` with
version 2 in 2016 :footcite:`falivene_sambvca_2016`. .

The approach of ᴍᴏʀғᴇᴜs differs somewhat from that in ref.
:footcite:`falivene_sambvca_2016` in that points are generated uniformly in the
test sphere rather than considering voxels. The numerical results with standard
settings are the same though as shown by benchmarks on complexes 1-18 from ref.
:footcite:`falivene_sambvca_2016`. Steric maps also match those in ref.
:footcite:`falivene_sambvca_2016`.

.. footbibliography::