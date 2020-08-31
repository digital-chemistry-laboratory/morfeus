================
Pyramidalization
================

Pyramidalization can be calculated as described by Agranat and Radhakrishnan
[1]_.

*******************
Command line script
*******************

******
Module
******

**********
Background
**********

Pyramidalization can be calculated for any tetracoordinate atom as described
in [1]_. An important alteration to the original recipe is made for extreme
cases of pyramidalization. As suggested by Tobias Hensch [ref], the sign of
the α angle is adjusted to be negative (:numref:`pyramidalization`). When the
angle between vector c and the sum of vectors a and b is acute, the α angle
will be taken as negative. Correspondingly, the pyramidalization value *P* is
taken as 2 - *P* and can therefore be larger than 1.

.. figure:: images/pyramidalization/pyramidalization.png
  :name: pyramidalization
  :width: 50%
  
  Definition of α angle as negative for extreme pyramidalization.


**********
References
**********

.. [1] Radhakrishnan, T. P.; Agranat, *I. Struct. Chem.* **1991**, *2*, 107