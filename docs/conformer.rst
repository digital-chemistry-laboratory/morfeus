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
conformers down to 3, which are are sorted based on energy. For full information, 
use ``help(ConformerEnsemble)`` or consult the API:
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
  interface is not efficient. The CREST__ program is recommended in these
  cases.

.. __: https://github.com/MolSSI/QCEngine
.. __: http://docs.qcarchive.molssi.org/projects/QCEngine/
.. __: https://github.com/grimme-lab/crest

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