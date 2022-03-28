===============
Conformer tools
===============

Conformer ensembles can be generated and handled with the conformer tools in
Mᴏʀғᴇᴜs. Generation of new conformers from SMILES strings relies on interfaces
to the RDKit__ and OpenBabel__.

.. __: https://www.rdkit.org
.. __: http://openbabel.org/wiki/Main_Page

OpenBabel and RDKit are avaible for install via, *e.g.*, the conda-forge:

.. code-block:: console

  conda install -c conda-forge rdkit
  conda install -c conda-forge openbabel

******
Module
******

Below is an example where the RDKit is used to generate conformers for
*n*-butane.

.. code-block:: python
  :caption: Example

  >>> ce = ConformerEnsemble.from_rdkit("CCCC", optimize="MMFF94")
  >>> ce
  ConformerEnsemble(50 conformers)
  >>> ce.prune_rmsd()
  ConformerEnsemble(3 conformers)
  >>> ce.sort()
  >>> ce.get_relative_energies()
  array([0.        , 0.78220907, 0.78220907])
  >>> ce.write_xyz("conformers.xyz")

Here, conformers are first generated with the RDKit and optimized with the
MMFF94 force field. Then a pruning based on RMSD brings the number of
conformers down to 3, which are are sorted based on energy. For full
information, use ``help(ConformerEnsemble)`` or consult the API:
:py:class:`ConformerEnsemble <morfeus.conformer.ConformerEnsemble>`.

####################
Conformer generation
####################

There are three different methods for generating conformers with Mᴏʀғᴇᴜs. The
available arguments are listed for the standalone functions.

.. list-table::
  :header-rows: 1

  * - Method
    - Standalone function
    - Class method
  * - RDKit
    - :py:func:`conformers_from_rdkit <morfeus.conformer.conformers_from_rdkit>`
    - :py:func:`ConformerEnsemble.from_rdkit <morfeus.conformer.ConformerEnsemble.from_rdkit>`
  * - Openbabel FF
    - :py:func:`conformers_from_ob_ff <morfeus.conformer.conformers_from_ob_ff>`
    - :py:func:`conformers_from_ob_ff <morfeus.conformer.conformers_from_ob_ff>`
  * - Openbabel GA
    - :py:func:`conformers_from_ob_ga <morfeus.conformer.conformers_from_ob_ga>`
    - :py:func:`conformers_from_ob_ga <morfeus.conformer.conformers_from_ob_ga>`

The conformers from OpenBabel are generally pruned with respect to RSMD while
those from RDKit require post-generation pruning. Here is an example of using
OpenBabel.

.. code-block:: python
  :caption: Conformers from OpenBabel

  >>> ce = ConformerEnsemble.from_ob_ga("CCCC")
  >>> ce
  ConformerEnsemble(3 conformers)

#####
CREST
#####

Output from the CREST__ program can be parsed. Mᴏʀғᴇᴜs does not currently have
any functionality for running CREST. The following files must be present in
the output folder: *crest_conformers.xyz*, *cre_members* and *crest.eneriges*.

.. tab:: Mᴏʀғᴇᴜs

  .. code-block:: python
    :caption: Conformer ensemble from CREST output

    >>> ce = ce = ConformerEnsemble.from_crest("crest_output_folder")
    >>> ce
    ConformerEnsemble(8 conformers)
    >>> ce.sort()
    >>> ce.get_relative_energies()
    array([0.   , 0.468, 0.829, 0.832, 0.834])
    >>> ce.get_degeneracies()
    array([ 4, 27,  8,  2,  1])
    >>> ce.boltzmann_weights()
    array([0.45642055, 0.20716615, 0.11264294, 0.11207402, 0.11169634])

.. tab:: CREST output

  .. code-block:: none

    Erel/kcal        Etot weight/tot  conformer     set   degen     origin
    1    0.000    -2.40317    0.05283    0.21084       1       4     mtd3
    2    0.000    -2.40317    0.05281                                mtd4
    3    0.002    -2.40317    0.05262                                mtd3
    4    0.003    -2.40316    0.05258                                mtd4
    5    0.468    -2.40242    0.02400    0.64605       2      27     mtd2
    6    0.468    -2.40242    0.02399                                mtd4
    7    0.468    -2.40242    0.02398                                mtd6
    8    0.468    -2.40242    0.02398                                mtd1
    9    0.469    -2.40242    0.02398                                mtd1
    10   0.469    -2.40242    0.02398                                mtd1
    11   0.469    -2.40242    0.02397                                mtd2
    12   0.469    -2.40242    0.02397                                mtd1
    13   0.469    -2.40242    0.02397                                mtd3
    14   0.469    -2.40242    0.02396                                mtd1
    15   0.469    -2.40242    0.02395                                mtd1
    16   0.469    -2.40242    0.02395                                mtd4
    17   0.469    -2.40242    0.02394                                mtd4
    18   0.469    -2.40242    0.02394                                mtd5
    19   0.469    -2.40242    0.02394                                mtd3
    20   0.470    -2.40242    0.02393                                mtd1
    21   0.470    -2.40242    0.02393                                mtd1
    22   0.470    -2.40242    0.02391                                mtd3
    23   0.470    -2.40242    0.02391                                mtd5
    24   0.470    -2.40242    0.02390                                mtd5
    25   0.471    -2.40242    0.02389                                mtd2
    26   0.471    -2.40242    0.02388                                mtd1
    27   0.471    -2.40242    0.02387                                mtd3
    28   0.472    -2.40242    0.02386                                mtd3
    29   0.472    -2.40242    0.02385                                mtd5
    30   0.473    -2.40242    0.02381                                mtd3
    31   0.473    -2.40242    0.02381                                mtd6
    32   0.829    -2.40185    0.01306    0.10429       3       8     mtd3
    33   0.829    -2.40185    0.01305                                mtd2
    34   0.829    -2.40185    0.01305                                mtd6
    35   0.829    -2.40185    0.01305                                mtd2
    36   0.829    -2.40185    0.01305                                mtd2
    37   0.830    -2.40185    0.01304                                mtd5
    38   0.832    -2.40184    0.01300                                mtd5
    39   0.832    -2.40184    0.01299                                mtd5
    40   0.832    -2.40184    0.01299    0.02588       4       2     mtd3
    41   0.837    -2.40184    0.01289                                mtd3
    42   0.834    -2.40184    0.01295    0.01295       5       1     mtd5

.. __: https://github.com/grimme-lab/crest

#############################
Boltzmann-weighted properties
#############################

A convenient way to calculate properties for the ensemble is to iterate over
the ``ConformerEnsemble`` object:

.. code-block:: python
  :caption: Boltzmann-weighted SASA

  >>> ce = ConformerEnsemble.from_rdkit("CCCO", optimize="MMFF94")
  >>> ce.prune_rmsd()
  >>> ce.sort()
  ConformerEnsemble(3 conformers)
  >>> for conformer in ce:
  >>> ... sasa = SASA(ce.elements, conformer.coordinates)
  >>> ... conformer.properties["sasa"] = sasa.area
  >>> ce.get_properties()
  {'sasa': array([221.26889622, 217.38905484, 216.53891818])}
  >>> ce.boltzmann_weights()
  array([0.56040173, 0.28260296, 0.15699531])
  >>> ce.boltzmann_statistic("sasa")
  219.4298571958332

The default of the function
:py:meth:`ConformerEnsemble.boltzmann_statistic <morfeus.conformer.ConformerEnsemble.boltzmann_statistic>`
is to calculate the Boltzmann average at 298.15 K, but this can be changed with
``temperature=<float>`` and ``statistic=<str>``, where "var" or "std" are
available. The temperature derivative of the Boltzmann average can also be
calculated with the method
:py:meth:`ConformerEnsemble.boltzmann_average_dT <morfeus.conformer.ConformerEnsemble.boltzmann_average_dT>`

############
RMSD pruning
############

Conformers are usually pruned on root mean square deviation in terms of (heavy)
atom coordinates to remove redundant structures which correspond to essentially
the same conformation. In Mᴏʀғᴇᴜs, this is achieved with the
:py:meth:`ConformerEnsemble.prune_rmsd <morfeus.conformer.ConformerEnsemble.prune_rmsd>`
method. By default, the ``AllChem.AlignMolConformers`` function from RDKit is
used to calculate the RMSD, but this can be changed with the keyword argument
``method=<str>``. The following options are available. spyrmsd__ needs to be
installed for that option to work.

.. list-table::
  :header-rows: 1

  * - method
    - symmetry
    - include_hs
  * - obrms-batch
    - Always
    - Never
  * - obrms-iter
    - Always
    - Never
  * - openbabel
    - Optional
    - Optional
  * - rdkit
    - Never
    - Never
  * - spyrmsd
    - Optional
    - Optional

The distinguishing factors are whether symmetry and non-heavy atoms are
considered when calculating the RMSD. For the ``method="openbabel"`` and
``method="spyrmsd"``, the keyword arguments ``symmetry=<bool>``and
``include_hs=<bool>`` are used to control the behavior. For the rest of the
methods, these arguments will be ignored. Pruning out conformers that are the
same by symmetry can lower the computational cost, but might also lead to
errors in Boltzmann weighting as degeneracy is not taken into account.

.. warning::

    ``include_hs`` and ``symmetry`` are ignored unless ``method`` is
    "openbabel" or "spyrmsd".

.. __: https://github.com/RMeli/spyrmsd

##############
Energy pruning
##############

Conformers are often pruned based on energy. The Boltzman weight for conformers
above 3 kcal/mol are expected to contribute in a neglible fashion to the
properties at room temperature. Therefore, the default of the
:py:meth:`ConformerEnsemble.prune_energy <morfeus.conformer.ConformerEnsemble.prune_energy>`
is to prune out all conformers above this energy. This can be changed with the
keyword argument ``threshold=<float>``.

.. code-block:: python
  :caption: Energy pruning

  >>> ce = ConformerEnsemble.from_rdkit("C1CCCCC1", optimize="MMFF94")
  >>> ce.prune_rmsd(method="obrms-batch")
  >>> ce.sort()
  >>> ce.get_relative_energies()
  array([0.       , 5.9297458])
  >>> ce.prune_energy()
  >>> ce.get_relative_energies()
  array([0.])

########################
Optimization and ranking
########################

Mᴏʀғᴇᴜs has an interface to QCEngine__ that allows calculation of single-point
energies and geometry optimizations for conformers. The ``program``, ``model``,
``keywords`` and ``local_options`` keyword arguments are all passed on to
QCEngine and more information can be found in the documentation__. Here is an
example, optimizing a conformer ensemble with GFN-FF and doing single points
with GFN2-xTB.

.. code-block:: python
  :caption: Optimization and ranking

  >>> ce = ConformerEnsemble.from_rdkit("CCCC", optimize="MMFF94")
  >>> ce.prune_rmsd()
  >>> ce.sort()
  >>> ce.get_relative_energies()
  array([0.        , 0.78220907, 0.78220907])
  # Optimize with GFN-FF
  >>> model={"method": "GFN-FF"}
  >>> ce.optimize_qc_engine(program="xtb", model=model, procedure="geometric")
  >>> ce.get_relative_energies()
  array([0.        , 0.45271867, 0.45271867])
  # Do single points with GFN2-xTB
  >>> model={"method": "GFN2-xTB"}
  >>> ce.sp_qc_engine(program="xtb", model=model)
  >>> ce.get_relative_energies()
  array([0.        , 0.63087431, 0.6308743 ])

.. note::

  Optimization of many molecules with many conformers through the QCEngine
  interface is not efficient. The CREST_ program is recommended in these
  cases.

.. __: https://github.com/MolSSI/QCEngine
.. __: http://docs.qcarchive.molssi.org/projects/QCEngine/

#######################
Enantiomeric conformers
#######################

Mᴏʀғᴇᴜs can handle degeneracy stemming from enatiomeric conformations to some
extent. This can be exemplified fro diethyl ether, which has three types of
conformers *tt* (1), *tg* (4) and *gg* (2), where *t* stands for *trans*, *g*
for *gauche* and the number in parenthesis is the total degeneracy of that
type. :footcite:`merrill_solvent_2020`

.. code-block:: python
  :caption: Degeneracies

  >>> ce = ConformerEnsemble.from_rdkit("CCOCC", optimize="MMFF94")
  >>> ce
  ConformerEnsemble(50 conformers)
  >>> ce.prune_rmsd()
  ConformerEnsemble(7 conformers)
  >>> ce.add_inverted()  # Invert conformers to add potentially missing enantiomers
  ConformerEnsemble(14 conformers)
  >>> ce.prune_rmsd()  # Prune out all duplicates added in previous step
  ConformerEnsemble(7 conformers)
  >>> ce.condense_enantiomeric()  # Condense enantiomers to single conformer
  ConformerEnsemble(4 conformers)
  >>> ce.sort()
  >>> ce.get_degeneracies()
  array([1, 2, 2, 2])
  >>> ce.get_relative_energies()
  array([0.        , 1.51922335, 1.51922335, 3.03207887])

We can see that Mᴏʀғᴇᴜs has indeed obtained three different types of conformers
(as seen by the energies) with degeneracies 1, 4 (2 + 2) and 2.

###########
Enantiomers
###########

Another situation is when there are enantiomers not due to conformation but to
configuration. In these cases, the methods above are not safe to use. If the
stereochemistry is given in the SMILES, RDKit will generate only that
enantiomer. However, if the sterochemistry is not given, both enantiomers will
be generated. Mᴏʀғᴇᴜs can prune out one of the enantiomers in this case to save
computational time if there is no interest in the actual stereochemistry. Here
is an example for alanine.

.. code-block:: python
  :caption: Generating specific enantiomer

  >>> ce = ConformerEnsemble.from_rdkit("C[C@@H](C(=O)O)N", optimize="MMFF94")
  >>> ce.prune_rmsd()
  >>> ce.sort()
  >>> ce.get_cip_labels()
  [('', 'S', '', '', '', '', '', '', '', '', '', '', ''),
   ('', 'S', '', '', '', '', '', '', '', '', '', '', ''),
   ('', 'S', '', '', '', '', '', '', '', '', '', '', '')]

.. code-block:: python
  :caption: Generating both enantiomers

  >>> ce = ConformerEnsemble.from_rdkit("CC(N)C(O)=O", optimize="MMFF94")
  >>> ce.prune_rmsd()
  >>> ce.sort()
  >>> ce.get_cip_labels()
  [('', 'R', '', '', '', '', '', '', '', '', '', '', ''),
   ('', 'S', '', '', '', '', '', '', '', '', '', '', ''),
   ('', 'S', '', '', '', '', '', '', '', '', '', '', ''),
   ('', 'R', '', '', '', '', '', '', '', '', '', '', ''),
   ('', 'S', '', '', '', '', '', '', '', '', '', '', ''),
   ('', 'R', '', '', '', '', '', '', '', '', '', '', ''),
   ('', 'S', '', '', '', '', '', '', '', '', '', '', ''),
   ('', 'R', '', '', '', '', '', '', '', '', '', '', '')]
  >>> ce.prune_enantiomers()
  >>> ce.get_cip_labels()
  [('', 'R', '', '', '', '', '', '', '', '', '', '', ''),
   ('', 'R', '', '', '', '', '', '', '', '', '', '', ''),
   ('', 'R', '', '', '', '', '', '', '', '', '', '', ''),
   ('', 'R', '', '', '', '', '', '', '', '', '', '', '')]

For more information, see the
:py:meth:`ConformerEnsemble.prune_enantiomers <morfeus.conformer.ConformerEnsemble.prune_enantiomers>`
method.

*******************
Command line script
*******************


.. footbibliography::
