================
Pyramidalization
================

Pyramidalization can be calculated as described by Agranat and Radhakrishnan
:footcite:`radhakrishnan_measures_1991`. 

.. note::
  Documentation to be completed.

*******************
Command line script
*******************

******
Module
******

**********
Background
**********

Pyramidalization can be calculated for any tetracoordinate atom as described in
:footcite:`radhakrishnan_measures_1991`:. An important alteration to the
original recipe is made for extreme cases of pyramidalization. As suggested by
Tobias Gensch [ref], the sign of the α angle is adjusted to be negative
(:numref:`fig_pyramidalization`). When the angle between vector c and the sum
of vectors a and b is acute, the α angle will be taken as negative.
Correspondingly, the pyramidalization value *P* is taken as 2 - *P* and can
therefore be larger than 1.

.. figure:: images/pyramidalization/pyramidalization.svg
  :name: fig_pyramidalization
  :width: 75%
  
  Definition of α angle as negative for extreme pyramidalization.

.. footbibliography::