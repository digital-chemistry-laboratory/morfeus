===============================
Local force constants
===============================

Local force constants can be calculated with the local modes method [1]_ or the
compliance matrix method [2]_. A Gaussian log file is currently required.
Increased accuracy might be achieved by supplying a fchk file with
high-precision normal modes and normal mode force constants.

************************
Preparing Gaussian input
************************

Depending on if you want to use the local modes method (recommended) or the 
compliance matrix method, different inputs for Gaussian are needed.

###########
Local modes
###########

The local modes method needs normal modes decomposed into contributions from
internal coordinates as well as normal mode force constants.

1. Normal mode decomposition in Gaussian
  This is the most simple way. Only the log file is required.
    ``#p freq(intmodes) iop(7/75=-1) iop(1/33=3)``
2. Normal mode decomposition from high precision normal modes
  This might increase the accuracy somewhat. Only log file is required.
    ``#p freq(intmodes,hpmodes) iop(1/33=3)``
3. Normal mode decomposition from very high precision normal modes and force
constants
  This is the highest accuracy mode. Both the log file and fchk file are
  required. The fchk file is generated with the Gaussian formchk program.
    ``#p freq(intmodes) iop(1/33=3)``

#################
Compliance matrix
#################

The compliance matrix methods needs the force constant matrix (Hessian).

1. Force constant matrix from the log file
    ``#p freq(intmodes) iop(1/33=3)``
2. Higher-accuracy force constant matrix from the fchk file
    ``#p freq(intmodes) iop(1/33=3)``
3. Highest-accuarcy force constant matrix from PES file.
    ``#p freq(intmodes) iop(1/33=3) iop(7/32=3)``

*******************
Command line script
*******************

Some dummy text here.

******
Module
******

The LocalForce class is provided to calculate and store the local force
constants.

.. code-block:: python
  :caption: Example with local modes method

  >>> from steriplus import LocalForce
  >>> lf = LocalForce("freq-lm.log")
  >>> fc = lf.get_local_force_constant([1, 2])
  >>> print(fc)
  5.364289643211871

.. code-block:: python
  :caption: Example with compliance matrix method

  >>> from steriplus import LocalForce
  >>> lf = LocalForce("freq-lm.log", fchk_file="freq.fchk", method="compliance")
  >>> fc = lf.get_local_force_constant([1, 2])
  >>> print(fc)
  5.364476039405804

For the local modes method, projection of imaginary frequencies can be
controlled with the ``project_imag=<bool>``. The cutoff for low-freqency modes
can be controlled with ``cutoff=<float>``. Choice of method is controlled with
``method=<str>`` using either ``local`` (default) or ``compliance``. File names
of any fchk file and PES are specified with the ``fchk_file=<str>`` and
``pes_file=<str>`` keywords.

**********
Background
**********

Local force constants describe the bond strength based on vibrational
frequencies. In the literature, there are two approaches to this, the local
modes method of Cremer [1]_ and the compliance matrix method championed by
Grunenberg [2]. They have been shown to be equivalent within numerical accuracy
[3]_. Steriplus can use either method, and they give almost identical results
for most bonds. The exception is when imaginary or very small vibrational
frequencies exist. In this case, the numerical stability of the local modes
approach can be improved by two methods: (1) projecting out normal modes with
imaginary frequencies and (2) raising the force constants of the low-frequency
modes to a cutoff value. Steriplus does this projection by default and uses a 
cutoff of 0.001 mDyne/Å for low-frequency modes. We therefore recommend the
local modes with default setting for this case. Expert users can turn off this
projection and alter the cutoff value.

Note that interactions involving imaginary modes (such as breaking/forming
bonds in transition states) cannot be assessed by the local force constants.

**********
References
**********

.. [1] Konkoli, Z.; Cremer, D. Int. J. Quantum Chem. 1998, 67, 1.
.. [2] Brandhorst, K.; Grunenberg, J. Chem. Soc. Rev. 2008, 37, 1558.
.. [3] Zou, W.; Kalescky, R.; Kraka, E.; Cremer, D. J. Chem. Phys. 2012, 137, 84114.
