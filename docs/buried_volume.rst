=============
Buried volume
=============

Buried volumes are implemented as described by Cavallo and co-workers [1]_. 

*******************
Command line script
*******************

The basic functionality is available through the command line script.

.. code-block:: console
  :caption: Example

  $ morfeus_buried_volume 1.xyz 1 --exclude 1 2 3 4 5 6 7
  V_bur (%): 29.6

--include_hs <True/False>
  Include H atoms (default: ``False``)
--radii <str>  
  Type of radii, ``crc`` or ``bondi`` (default)
--radii_scale <float>
  Scaling of the radii (default: ``1.17``)
-r <float>
  Radius of sphere (default: ``3.5``)
--steric_map <list>
  Save steric map with supplied atoms to define z axis.
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
ligand of a transition metal complex [1]_. A web tool to calculate buried 
volumes, SambVca, was made available for scientific purposes by Cavallo and 
co-workers in 2009 [2]_ with version 2 in 2016 [1]_.

The approach of ᴍᴏʀғᴇᴜs differs somewhat from that in ref. [1]_ in that points
are generated uniformly in the test sphere rather than considering voxels. The 
numerical results with standard settings are the same though as shown by
benchmarks on complexes 1-18 from ref. [1]_. Steric maps also match those in 
ref. [1]_.

.. figure:: benchmarks/buried_volume/correlation.png



.. todo::
  Correlation to BVs of ref [X]
  .. figure:: 

**********
References
**********

.. [1] Falivene, L.; Credendino, R.; Poater, A.; Petta, A.; Serra, L.;
       Oliva, R.; Scarano, V.; Cavallo, L. *Organometallics* **2016**, *35*,
       2286.
.. [2] Poater, A.; Cosenza, B.; Correa, A.; Giudice, S.; Ragone, F.;
       Scarano, V.; Cavallo, L. Eur. J. *Inorg. Chem.* **2009**, *2009*, 1759.