=============
Buried volume
=============

Buried volumes are implemented as described by Cavallo and co-workers
:cite:`falivene_sambvca_2016`. 

*******************
Command line script
*******************

The basic functionality is available through the command line script.

.. code-block:: console
  :caption: Example

  $ morfeus_buried_volume 1.xyz 1 --exclude 1 2 3 4 5 6 7
  V_bur (%): 29.6

The first argument is the structure file, and the second argument the atom
index of the metal atom. The list of excluded atoms is also needed (must
contain at least the metal atom).

--density <float>
  Density of sphere grid (default 0.001)
--exclude <list>
  List of atoms to exclude from the calculation
--include_hs <True/False>
  Include H atoms (default: False)
--radii <str>  
  Type of radii, "crc" or "bondi" (default)
--radii_scale <float>
  Radii scale factor (default: 1.17)
--radius <float>
  Radius of sphere (default: 3.5)
--steric_map <list>
  Draw steric map with specified atoms to define z axis.
  Map is saved as "steric_map.png"

More information can be found with ``morfeus_buried_volume --help``

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

``bv.plot_3D()`` gives a three-dimensional representation of the buried volume.
Plots can be saved by passing the ``filename=<str>`` keyword argument.

By default, hydrogen atoms are excluded in the calculation. They can be added
by giving the keyword argument ``indclude_hs=True``. The default sphere radius
is 3.5 Å, but can be changed with ``radius=<float>``. Default radii type is
Bondi which are scaled by a factor of 1.17. This can be changed with
``radii_type=<str>``, choosing etiher ``crc`` or ``bondi`` and
``radii_scale=<float>``. Custom radii can be supplied as a list with
``radii=<list>``.

For more information, use ``help(BuriedVolume)`` or see the API:
:py:class:`morfeus.morfeus.BuriedVolume`

**********
Background
**********

The percent of buried volume is a measure of the steric hindrance induced by a 
ligand of a transition metal complex :cite:`falivene_sambvca_2016`. A web tool
to calculate buried volumes, SambVca, was made available for scientific
purposes by Cavallo and co-workers in 2009 :cite:`poater_sambvca_2009` with
version 2 in 2016 :cite:`falivene_sambvca_2016`. .

The approach of ᴍᴏʀғᴇᴜs differs somewhat from that in ref.
:cite:`falivene_sambvca_2016` in that points are generated uniformly in the
test sphere rather than considering voxels. The numerical results with standard
settings are the same though as shown by benchmarks on complexes 1-18 from ref.
:cite:`falivene_sambvca_2016`. Steric maps also match those in ref.
:cite:`falivene_sambvca_2016`.

.. todo::
  Correlation to BVs of ref [X]
  figure

**********
References
**********

.. bibliography:: refs.bib
  :style: unsrt
  :filter: docname in docnames