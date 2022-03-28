===
XTB
===

Electronic parameters can be calculated at the GFN1-xTB or GFN2-xTB levels
using the xtb-python__ interface (which needs to be installed).

******
Module
******

The XTB class is used to calculate electronic properties. Simple quantities
such as the ionization potential, electron affinity, HOMO, LUMO energies and
dipole moment are available as well as the atomic charges and bond orders.

.. code-block:: python
  :caption: Example

  >>> from morfeus import read_xyz
  >>> from morfeus import XTB
  >>> elements, coordinates = read_xyz("ammonia.xyz")
  >>> xtb = XTB(elements, coordinates)
  >>> xtb.get_ip()
  16.49875404054677
  >>> xtb.get_ip(corrected=True)
  11.652754040546771
  >>> xtb.get_ea()
  -5.736519991532061
  >>> xtb.get_homo()
  -0.39012010481047976
  >>> xtb.get_charges()
  {1: -0.42539265,
   2: 0.14180091,
   3: 0.14179421,
   4: 0.14179754}
  >>> xtb.get_bond_order(1, 2)
  0.9786781554103448
  >>> xtb.get_dipole()
  array([0.48187417, 0.06877519, 0.55556546])

In addition, global and local descriptors from conceptual density functional
theory can also be calculated.

.. code-block:: python
  :caption: Example

  >>> from morfeus import read_xyz
  >>> from morfeus import XTB
  >>> elements, coordinates = read_xyz("ammonia.xyz")
  >>> xtb = XTB(elements, coordinates)
  >>> xtb.get_global_descriptor("electrophilicity", corrected=True)
  0.00643909828825333
  >>> xtb.get_global_descriptor("nucleophilicity", corrected=True)
  -11.652754040546771
  >>> xtb.get_fukui("electrophilicity")
  {1: -0.20661935,
   2: -0.26445605,
   3: -0.26448747,
   4:-0.26443713}
  >>> xtb.get_fukui("nucleophilicity")
  {1: -0.42294271,
   2: -0.19234974,
   3: -0.19235729,
   4:-0.19235026}

The version of GFNX-xTB can be set using ``version=<int>`` with versions 1 and
2 currently supported. A correction term can be applied for the calculation of
the ionization potential and electron affinity using ``corrected=True``, which
also affects some of the global and local descriptors. For a full list of
descriptors and their definitions, see the Background_.

For more information, use ``help(XTB)`` or consult the API:
:py:class:`XTB <morfeus.xtb.XTB>`.

*******************
Command line script
*******************

The command line script gives access to the basic functionality from the
terminal.

.. tab:: Get charge

  .. code-block:: console

    $ morfeus xtb Et.gjf - - get_charges - 1
    0.03220302786615441

.. tab:: Change version

  .. code-block:: console

    $ morfeus xtb Et.gjf - --version='"1"' - get_charges - 1
    0.02564834649261168

**********
Background
**********

ᴍᴏʀғᴇᴜs can compute both simple electronic parameters such as charges, HOMO
and LUMO energies and bond orders, as well as descriptors from conceptual
density functional theory :footcite:`domingo_applications_2016`.
The following global descriptors are available:

* Electrophilicity: :math:`\omega`
* Nucleophilicity: :math:`N`
* Electrofugality: :math:`\nu_{electrofugality}`
* Nucleofugality: :math:`\nu_{nucleofugality}`

They are calculated according to the following definitions
:footcite:`domingo_applications_2016,ayers_indices_2005`:

.. math::

  \omega &= \frac{(IP + EA)^2}{8(IP - EA)} = \frac{\mu^2}{2\eta}

  N &= -IP

  \nu_{electrofugality} &= \frac{(IP - 3EA)^2}{8(IP - EA)} = -EA + \omega

  \nu_{nucleofugality} &= \frac{(3IP - EA)^2}{8(IP - EA)} = IP + \omega

Where :math:`IP` is the ionization potential, :math:`EA` is the electron
affinity, :math:`\mu` is the chemical potential and :math:`\eta` is the
hardness given by

.. math::

  \mu &= - \frac{IP + EA}{2}

  \eta &= IP - EA

The Fukui coefficients are calculated calculated via the finite differences
approach using the atomic charges from *xtb*. These include:

* Electron removal: :math:`f^-`
* Electron addition: :math:`f^+`
* Radical attack: :math:`f`
* Dual descriptor: :math:`f^{(2)}`

Which are calculated as follows.

.. math::

  f^- &= q_{N-1} - q_{N}

  f^+ &= q_{N} - q_{N+1}

  f &= (q_{N-1} - q_{N+1}) / 2

  f^{(2)} &= f^+ - f^- = 2 q_{N} - q_{N+1} - q_{N-1}

The Fukui coefficient for electron removal is also called the coefficient for
electrophilic attack and is a measure of nucleophilicity. The coefficient for
electron addition is also called the coefficient for nucleophilic attack and is
a measure of electrophilicity. The somewhat unintuitive names is due to the
notion that *another* molecule would attack as a nucleophile/electrophile. The
coefficient for radical attack is often used for radical reactivity. In
addition, the local electrophilicity (:math:`l_{\omega}`) and nucleophilicity
(:math:`l_N`) are also available and calculated as
:footcite:`oller_global_2018`.

.. math::

  l_{\omega} &= - \frac{\mu}{\eta}f + \frac{1}{2}(\frac{\mu}{\eta})^2 f^{(2)}

  l_N &= f^-

The ionization potentials and electron affinities calculated with *xtb* can be
corrected using the empirical terms given by Grimme and co-workers
:footcite:`neugebauer_benchmark_2020`.

.. footbibliography::

.. __: https://github.com/grimme-lab/xtb-python/
