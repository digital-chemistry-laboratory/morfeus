=====================
Local force constants
=====================

Local force constants can be calculated with the local modes method [1]_ or
the compliance matrix method [2]_. Steriplus can use the output of Gaussian_,
xtb_, or UniMoVib_ programs.

===============
Preparing input
===============

The LocalForce class needs input from quantum-chemical frequency calculations.
Therefore, the instructions are a bit more involved. This input can either be
read using the ``load_file`` method, and Steriplus can also provide by 
built-in normal mode analysis and internal coordinate codes.

***********
Local modes
***********

The local modes approach needs the normal modes in terms of internal
coordinates, here called *internal modes*, as well as the normal mode force
constants. These can either be read from a file or calculated with Steriplus.
For example, the internal modes can be computed from the normal modes and the
Wilson B matrix. For local frequencies, additional input is needed.

.. figure:: images/local_force/diagram_local.svg
  
  Input needed for the local modes approach.

*****************
Compliance matrix
*****************

The compliance matrix method needs the Hessian matrix and the Wilson B matrix.
For local frequencies, the elements are needed (to get the atomic masses).

.. figure:: images/local_force/diagram_compliance.svg
  
  Input needed for the compliance matrix method.

*****************
Recommended input
*****************

The table below summarizes the recommended input to Steriplus from the
supported programs. More details are given below.

.. table:: Recommended input.
  :widths: auto
  :align: center

  =========== ======== ========= =======
  Method      Gaussian UniMoVib  xtb
  =========== ======== ========= =======
  Local modes log      log       hessian
  Compliance  fchk     umv/local hessian
  =========== ======== ========= =======

**********************
Geometry optimizations
**********************

Local force constants and frequencies are based on the harmonic approximation
that is valid only at stationary points (minima and transition states) on the
potential energy surface. Therefore, the geometry optimizations need to be of
good quality. For difficult examples on flat surfaces, this can mean
increasing the convergence criteria of the quantum-chemical program.
Vibrations with small imaginary frequencies should be eliminated as much as
possible. The local modes method is somewhat robust to the presence of these
vibrations with the standards settings, while the they can result in artifacts
with the compliance matrix method. Therefore, the local modes method is
recommended in these cases.

For transition states, the imaginary mode corresponding to the reaction is
projected out. This also means that forces involving theseatoms are
meaningless and should not be used. Only the local modes method can be used
for transition states.

.. warning:: Transition states

  Transition states can only be treate with the local modes method. Force
  constants and frequencies of the atoms corresponding to the imaginary mode
  should not be used.

********
Gaussian
********

Gaussian can provide the input for both the local modes and the compliance
matrix methods directly. For the local modes method, it is recommended to
use the log file, while for the the compliance matrix method, the fchk file
is neeed.

==== =========================================
log  ``freq=intmodes iop(7/75=1) iop(1/33=3)``
fchk ``freq=intmodes``
==== =========================================

``iop(7/75=1)`` makes Gaussian print the internal modes with full accuracy,
while ``iop(1/33=3)`` triggers printing of the Wilson B matrix. The table
below lists all the quantities read by each approach.

.. table:: Quantities from Gaussian files.
  :widths: auto
  :align: center

  ==================== === ===
  Quantity             log fck
  ==================== === ===
  Elements              x   x
  Coordinates           x   x
  Internal modes        x   
  Internal coordinates  x   x
  B Matrix              x   
  Normal modes              x
  Hessian                   x
  Force constants       x   x
  ==================== === ===

.. note:: Additional accuracy
  
  Additional accuracy for local modes calculations be achieved by loading
  also the fchk file which contains the normal mode force constants to
  higher accuracy.

********
UniMoVib
********

The UniMoVib program can do vibrational analysis for a range of different
quantum-chemical programs. The user needs to define internal coordinates in
addition to the information that can be read from the files. For the local
modes method, the *log* file is recommended. For the compliance matrix method,
either the *umv* or the *local* files are recommended.


.. table:: Quantities from UniMoVib files.
  :widths: auto
  :align: center

  ==================== ===== === ===
  Quantity             local log umv
  ==================== ===== === === 
  Elements             x     x   x
  Coordinates          x     x   x
  Internal modes          
  Internal coordinates 
  B Matrix                
  Normal modes         x     x
  Hessian              x         x  
  Force constants            x
  ==================== ===== === === 

The *log* file is the standard output of the program. The *umw* file and the
*local* file can be generated by specifying them in the input file

.. code-block::
  :emphasize-lines: 5, 6
  :caption: Example UniMoVib input file

  a test job

  $contrl
    qcprog="gaussian"
    iflocal=.t.
    ifsave=.t.
  $end

  $qcdata
    fchk="freq.fchk"
  $end

***
xtb
***

The xtb program can provide the Hessian and the normal modes and normal mode
force constants. The files *hessian* and *xtb_normalmodes* are generated by
the xtb program as a results of a frequency calculation. The recommended
approach for both the local modes method and the compliance method matrix is
to read the *hessian* file.

.. table:: Quantities from xtb files.
  :widths: auto
  :align: center

  ==================== ======= ===========
  Quantity             hessian normalmodes
  ==================== ======= ===========
  Elements             
  Coordinates          
  Internal modes          
  Internal coordinates 
  B Matrix                
  Normal modes                 x
  Hessian              x       
  Force constants              x
  ==================== ======= =========== 

.. warning:: Linear molecules

  xtb 6.2 has a bug which gives the wrong number of normal modes for linear
  molecules in the *xtb_normalmodes* file. Therefore, the approach of reading
  the Hessian and doing a normal mode analysis with Steriplus is recommended.  

*******************
Command line script
*******************

The command line script provides access to the basic functionality through
the terminal.

.. code-block:: console
  :caption: Example of single internal coordinate
  
  $ steriplus_local_force freq-hp.log -a 1 2
  5.364

.. code-block:: console
  :caption: Example of report
  
  $ steriplus_local_force freq-hp.log
    Atom_1    Atom_2    Atom_3    Atom_4                           Force constant(mDyne/Å)                       Frequency (cm^-1)
         1         2                                                                 5.364                                    3252
         1         3                                                                 5.364                                    3252
         1         4                                                                 5.364                                    3252
         1         5                                                                 5.364                                    3252

-a <list>
  List of atoms in the bond/internal coordinate.
--cutoff <float>
  Cutoff value for low-frequency modes (default:0.001)
--fchk_file <str>
  Name of Gaussian fchk file
--method <str>
  Method: "local" (default) or "compliance"
--no_project_imaginary
  Flag to disable projection of imaginary modes
--pes_file <str>
  Name of Gaussian PES file
  
More information is given with ``steriplus_local_force --help``

******
Module
******

The LocalForce class is provided to calculate and store the local force
constants.

.. code-block:: python
  :caption: Example with Gaussian and local modes method.

  >>> from steriplus import LocalForce
  >>> lf = LocalForce()
  >>> lf.load_file("freq-lm.log", "gaussian", "log")
  >>> lf.compute_local()
  >>> lf.compute_frequencies()
  >>> fc = lf.get_local_force_constant([1, 2])
  >>> print(fc)
  5.364289643211871
  >>> freq = lf.get_local_frequency([1, 2])
  >>> print(freq)
  3129.3126301763527
  
.. code-block:: python
  :caption: Example with Gaussian and compliance matrix method.

  >>> from steriplus import LocalForce
  >>> lf = LocalForce()
  >>> lf.load_file("freq.fchk", "gaussian", "fchk")
  >>> lf.compute_local()
  >>> lf.compute_frequencies()
  >>> fc = lf.get_local_force_constant([1, 2])
  >>> print(fc)
  5.364398642985929
  >>> freq = lf.get_local_frequency([1, 2])
  >>> print(freq)
  3129.352986019491

.. code-block:: python
  :caption: Example with xtb and local modes method.

  >>> from steriplus import LocalForce, read_xyz
  >>> elements, coordinates = read_xyz("xtbopt.xyz")
  >>> lf = LocalForce(elements, coordinates)
  >>> lf.load_file("hessian", "xtb", "hessian")
  >>> lf.normal_mode_analysis()
  >>> lf.detect_bonds()
  >>> print(lf.internal_coordinates)
  [Bond(1, 4), Bond(1, 3), Bond(1, 2), Bond(1, 5)]
  >>> lf.compute_local()
  >>> lf.compute_frequencies()
  >>> fc = lf.get_local_force_constant([1, 2])
  >>> print(fc)
  5.190222259808879
  >>> freq = lf.get_local_frequency([1, 2])
  >>> print(freq)
  3078.130379468432

.. code-block:: python
  :caption: Example with UniMoVib and the local modes method.

  >>> from steriplus import LocalForce
  >>> lf = LocalForce()
  >>> lf.load_file("job.out", "unimovib", "log")
  >>> lf.detect_bonds()
  >>> lf.compute_local()
  >>> lf.compute_frequencies()
  >>> fc = lf.get_local_force_constant([1, 2])
  >>> print(fc)
  5.364347084281302
  >>> freq = lf.get_local_frequency([1, 2])
  >>> print(freq)
  3129.337947449028

.. code-block:: python
  :caption: Example with adding internal coordinates manually
  :emphasize-lines: 4-7

  >>> from steriplus import LocalForce
  >>> lf = LocalForce()
  >>> lf.load_file("job.out", "unimovib", "log")
  >>> lf.add_internal_coordinate([1, 2])
  >>> lf.add_internal_coordinate([1, 2, 3])
  >>> print(lf.internal_coordinates)
  [Bond(1, 2), Angle(1, 2, 3)]
  >>> lf.compute_local()
  >>> lf.compute_frequencies()
  >>> fc = lf.get_local_force_constant([1, 2])
  >>> print(fc)
  5.364347084281298
  >>> freq = lf.get_local_frequency([1, 2])
  >>> print(freq)
  3129.337947449028
  >>> lf.print_report(angles=True, angle_units=True)
  Coordinate                            Force constant (mDyne/Å, mDyne Å rad^(-2))             Frequency (cm^-1)
  Bond(1, 2)                                                                 5.364                          3129
  Angle(1, 2, 3)                                                             2.416                          1687

For the local modes method, projection of imaginary frequencies can be
controlled with the ``project_imag=<bool>`` keyword to the ``compute_local``
method. The cutoff for low-freqency modes can be controlled with 
``cutoff=<float>``. Internal coordinates can be added with 

For more detailed information, use ``help(LocalForce)`` or see the API:
:py:class:`steriplus.steriplus.LocalForce`

**********
Background
**********

Local force constants describe the bond strength based on vibrational
frequencies. There are two approachces in the literature, the local modes
method of Cremer [1]_ and the compliance matrix method of Grunenberg [2]_.
They have been shown to be equivalent within numerical accuracy [3]_.
Steriplus can use either method, and they give almost identical results for
most cases. The exception is when there are modes with imaginary or very small
frequencies exist. In this case, the numerical stability of the local modes
approach can be improved by two methods: (1) projecting out normal modes with
imaginary frequencies and (2) raising the force constants of low-frequency
modes to a cutoff value. Steriplus does this projection by default and uses a 
cutoff of 0.001 mDyne/Å for low-frequency modes. We therefore recommend local
modes with default settings as the most robust method in these cases.
Expert users can turn off the projection and alter the cutoff value.

Note that interactions involving imaginary modes (such as breaking/forming
bonds in transition states) cannot be assessed by the local force constants.

Steriplus has been benchmarked against the local force constants and
frequencies for small organic molecules given by Cremer [3]_. 

.. figure:: benchmarks/local_force/benchmark.png
  
  Benchmark of local force constants and frequencies against data from Table 1
  of ref. [3]_. Data obtained with Gaussian log file and the local modes
  method.

**********
References
**********

.. [1] Konkoli, Z.; Cremer, D. Int. J. Quantum Chem. 1998, 67, 1.
.. [2] Brandhorst, K.; Grunenberg, J. Chem. Soc. Rev. 2008, 37, 1558.
.. [3] Zou, W.; Kalescky, R.; Kraka, E.; Cremer, D. J. Chem. Phys. 2012, 137, 84114.

.. _Gaussian: https://gaussian.com/
.. _UniMoVib: https://github.com/zorkzou/UniMoVib
.. _xtb: https://xtb-docs.readthedocs.io/en/latest/contents.html
