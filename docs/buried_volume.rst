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

  $ steriplus_buried_volume 1.xyz 1 --exclude 1 2 3 4 5 6 7
  V_bur (%): 29.6

--include_hs  Whether to include H atoms in the calculations (default: False)
--radii  Type of radii, "crc" or "bondi" (default: "crc")
--radii_scale  Scaling of the vdW radii (default: 1.17)
-r  Radius of sphere (default: 3.5)
--steric_map  Plots a steric map with the supplied atoms to define the z axis.
  Map is saved as "steric_map.png"

More information can be found with ``steriplus_buried_volume --help``

******
Module
******

The BuriedVolume class calculates the BV. Steric maps can also be plotted.

.. code-block:: python
  :caption: Example

  >>> from steriplus import BuriedVolume, read_xyz
  >>> elements, coordinates = read_xyz("1.xyz")
  >>> bv = BuriedVolume(elements, coordinates, 1, exclude_list=[1, 2, 3, 4, 5, 6, 7])
  >>> print(bv.buried_volume)
  0.2962110976518822
  >>> bv.print_report()
  V_bur (%): 29.6
  >>> bv.plot_steric_map([14])

.. image:: images/steric_map.png

``bv.plot_3D()`` gives a three-dimensional representation of the buried volume.
Plots can be saved by passing the ``filename={filename}`` keyword argument.

By default, hydrogen atoms are excluded in the calculation. They can be added
by giving ``indclude_hs=True``. The default radius sphere is 3.5 Å, but can be
changed with ``radius={radius}``. Default radii type is Bondi which are scaled
by a factor of 1.17. This can be changed with ``radii_type={"crc", "bondi"}`` 
and ``radii_scale={radii_scale}``. Custom radii can be supplied as a list with
``radii={radii_list}``.

**********
Background
**********

The percent of buried volume is a measure of the steric hindrance induced by a 
ligand of a transtion metal complex [1]_. A web tool to calculate buried 
volumes, SambVca was made available for scientific purposes by Cavallo and 
co-workers in 2009 [2]_ with version 2 in 2016 [1]_.

The approach of Steriplus differs somewhat from that in ref. [1]_ in that points
are generated uniformly in the test sphere rather than considering voxels. The 
numerical results with standard settings are the same though as shown by
benchmarks on complexes 1-18 from ref. [1]_. Steric maps also match those in 
ref. [1]_.

.. figure:: benchmarks/buried_volume/correlation.png

  Correlation to BVs of ref [X]_

.. figure:: 
  
  


**********
References
**********

.. [1] Falivene, L.; Credendino, R.; Poater, A.; Petta, A.; Serra, L.;
       Oliva, R.; Scarano, V.; Cavallo, L. Organometallics 2016, 35, 2286.
.. [2] Poater, A.; Cosenza, B.; Correa, A.; Giudice, S.; Ragone, F.;
       Scarano, V.; Cavallo, L. Eur. J. Inorg. Chem. 2009, 2009, 1759.