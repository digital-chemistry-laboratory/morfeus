===
XTB
===

Electronic parameters can be calculated at the GFN1-xTB, GFN2-xTB or PTB levels
using the xtb__ package (which needs to be installed).

******
Module
******

The XTB class is used to calculate electronic properties, such as, among others, bond orders, 
atomic charges, total, HOMO and LUMO energies, dipole moments, ionization potential, electron affinity.

.. code-block:: python
  :caption: Example

  >>> from morfeus import XTB, read_xyz
  >>> elements, coordinates = read_xyz("ammonia.xyz")
  >>> xtb = XTB(elements, coordinates)
  >>> xtb.get_bond_order(1, 2)
  0.978270297857
  >>> xtb.get_charges()
  {1: -0.4344091, 2: 0.14479795, 3: 0.14480895, 4: 0.1448022}
  >>> xtb.get_energy()
  -4.42584832
  >>> xtb.get_homo()
  -0.3847857
  >>> xtb.get_dipole()
  array([ 0.38617738,  0.48077917, -0.31171428])
  >>> xtb.get_ip()
  11.6047
  >>> xtb.get_ip(corrected=False)
  16.4502
  >>> xtb.get_ea()
  -11.121

In addition, global and local descriptors from conceptual density functional
theory can also be calculated.

.. code-block:: python
  :caption: Example

  >>> xtb.get_global_descriptor("electrophilicity")
  0.0012869003485041107
  >>> xtb.get_global_descriptor("nucleophilicity")
  -11.6047
  >>> xtb.get_fukui("electrophilicity")
  {1: 0.205, 2: 0.265, 3: 0.265, 4: 0.265}
  >>> xtb.get_fukui("nucleophilicity")
  {1: 0.428, 2: 0.191, 3: 0.191, 4: 0.191}

Properties related to the solvation are also available with the XTB class.

.. code-block:: python
  :caption: Example

  >>> xtb = XTB(elements, coordinates, solvent="water")
  >>> xtb.get_solvation_energy(unit="kcal/mol")
  -4.644374184273
  >>> xtb.get_solvation_h_bond_correction(unit="kcal/mol")
  -2.808641396951

The version of GFNn-xTB or the PTB method can be set using ``method=<int>|<str>`` with versions 1 and
2 and the additional method "ptb" currently supported. The molecular charge, number of unpaired electrons 
and solvent are also set while initializing the class, with ``charge=<int>``, ``n_unpaired=<int>`` and 
``solvent=<str>`` respectively.

A empirical correction term, given by Grimme and co-workers :footcite:`neugebauer_benchmark_2020`, 
is applied for the calculation of the ionization potential and electron affinity, which also 
affects some of the global and local descriptors. It can be disabled with ``corrected=False``. 

For more information and a list of all available electronic properties, 
use ``help(XTB)`` or consult the API:
:py:class:`XTB <morfeus.xtb.XTB>`.

*******************
Command line script
*******************

The command line script gives access to the basic functionality from the
terminal.

.. tab:: Get charge

  .. code-block:: shell

    $ morfeus xtb ammonia.xyz - - get_charges - 1
    -0.4344091

.. tab:: Change version

  .. code-block:: shell

    $ morfeus xtb ammonia.xyz - --method='"1"' - get_charges - 1
    -0.57697413

**********
Background
**********

###############################
Descriptors from conceptual DFT
###############################

ᴍᴏʀғᴇᴜs can compute both simple electronic parameters such as charges, HOMO
and LUMO energies, and bond orders, as well as descriptors from conceptual
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

  \nu_{electrofugality} &= \frac{(3IP - EA)^2}{8(IP - EA)} = IP + \omega

  \nu_{nucleofugality} &= \frac{(IP - 3EA)^2}{8(IP - EA)} = -EA + \omega

Where :math:`IP` is the ionization potential, :math:`EA` is the electron
affinity, :math:`\mu` is the chemical potential and :math:`\eta` is the
hardness given by

.. math::

  \mu &= - \frac{IP + EA}{2}

  \eta &= IP - EA

The Fukui coefficients are calculated via the finite differences
approach using the atomic charges from *xtb*. These include:

* Electron removal (nucleophilicity): :math:`f^-`
* Electron addition (electrophilicity): :math:`f^+`
* Radical attack: :math:`f`
* Dual descriptor: :math:`f^{(2)}`

which are calculated as follows:

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
:footcite:`oller_global_2018`:

.. math::

  l_{\omega} &= - \frac{\mu}{\eta}f + \frac{1}{2}(\frac{\mu}{\eta})^2 f^{(2)}

  l_N &= f^-

####################
Solvation properties
####################

The solvation free energy :math:`\Delta G_{solv}` is calculated with *xtb* 
using the ALPB model :footcite:`ehlert_robust_2021`. Because the solvent is 
described implicitly as a dielectric medium, the hydrogen bonding (HB) between solute 
and solvent needs to be accounted for separately. The HB correction 
is calculated with: 

.. math::

  \Delta G^{HB} = \sum_{A}^{N} \Delta G_{A}^{HB}

where :math:`\Delta G_{A}^{HB}` is the atomic contribution to the HB correction
and approximated in the ALPB model as :footcite:`ehlert_robust_2021`:

.. math::

  \Delta G_{A}^{HB} ≈ - g_{A}^{HB} q_{A}^{2} \frac{\sigma_{A}}{4\pi(R_{A}^{surf})^{2}}

where :math:`g_{A}^{HB}` is the HB strength of the atom, :math:`q_{A}` its charge, 
:math:`\sigma_{A}` its SASA, and :math:`R_{A}^{surf}` the atom vdW radius
combined with the solvent probe radius.


.. footbibliography::

.. __: https://github.com/grimme-lab/xtb
