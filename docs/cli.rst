======================
Command line interface
======================

The command line interface (CLI) for ᴍᴏʀғᴇᴜs is powered by `Python Fire`__. In
this way, allmost all of the functionality of the module system is made
available also in the command line.

.. __: https://github.com/google/python-fire

************
Introduction
************

The entry point to the CLI is the command line script ``morfeus``. Running this
command in the terminal gives a help page with all the possible subcommands.

.. tab:: Input 1

  .. code-block:: console

    $ morfeus

.. tab:: Output 1

  .. code-block:: none

    NAME
      morfeus

    SYNOPSIS
      morfeus COMMAND

    COMMANDS
      COMMAND is one of the following:

       buried_volume
         CLI for buried volume.

       cone_angle
         CLI for cone angle.

       conformer
         CLI for cone angle.

       dispersion
         CLI for dispersion descriptor.

       local_force
         CLI for local force.

       pyramidalization
         CLI for pyramidalization.

       sasa
         CLI for solvent accessible surface area.

       sterimol
         CLI for Sterimol.

       visible_volume
         CLI for visible volume.

       xtb
         CLI for XTB.

There is also a help page for each subcommand:

.. tab:: Input 2

  .. code-block:: console

    $ morfeus sterimol -- --help

.. tab:: Output 2

  .. code-block:: none

    NAME
        morfeus sterimol - CLI for Sterimol.

    SYNOPSIS
        morfeus sterimol FILE

    DESCRIPTION
        CLI for Sterimol.

    POSITIONAL ARGUMENTS
        FILE
            Type: str
            Geometry file

    NOTES
        You can also use flags syntax for POSITIONAL ARGUMENTS

All the subcommands serve as entry points to the corresponding Python classes
used in the module system. For example, ``morfeus sterimol tBu.xyz -`` creates
a partially instantiated :py:class:`Sterimol <morfeus.sterimol.Sterimol>`
object where the ``elements`` and ``coordinates`` arguments have been set from
the geometry file. To find out what other argument can be given to this object,
we run

.. tab:: Input 3

  .. code-block:: console

    $ morfeus sterimol tBu.xyz - -- --help

.. tab:: Output 3

  .. code-block:: none

    NAME
        morfeus sterimol tBu.xyz - partial(func, *args, **keywords) - new function with partial application of the given arguments and keywords.

    SYNOPSIS
        morfeus sterimol tBu.xyz - GROUP | COMMAND | --dummy_index=DUMMY_INDEX --attached_index=ATTACHED_INDEX <flags>

    DESCRIPTION
        partial(func, *args, **keywords) - new function with partial application of the given arguments and keywords.

    ARGUMENTS
        DUMMY_INDEX
            Type: int
        ATTACHED_INDEX
            Type: typing.Union[int, ...

    FLAGS
        --radii=RADII
            Type: Optional[typing.Union[int, float, complex, str, bytes, numpy.ge...
            Default: None
        --radii_type=RADII_TYPE
            Type: str
            Default: 'crc'
        --n_rot_vectors=N_ROT_VECTORS
            Type: int
            Default: 3600
        --excluded_atoms=EXCLUDED_ATOMS
            Type: Optional[t...
            Default: None
        --calculate=CALCULATE
            Type: bool
            Default: True

    GROUPS
        GROUP is one of the following:

         args

         keywords

    COMMANDS
        COMMAND is one of the following:

         func
           Performs and stores results of Sterimol calculation.

We learn that we need to supply the two arguments: ``dummy_index`` and
``attached_index``. Optional keyword arguments can also be given. After
supplying the additional arguments, we get back a fully instantiated object.
Getting the help from this object returns the available commands and what
values are available.

.. tab:: Input 4

  .. code-block:: console

      $ morfeus sterimol tBu.xyz - 1 2 --radii_type=bondi - -- --help

.. tab:: Output 4

  .. code-block:: none

    NAME
        morfeus sterimol tBu.xyz 1 2 --radii_type=bondi - Performs and stores results of Sterimol calculation.

    SYNOPSIS
        morfeus sterimol tBu.xyz - 1 2 --radii_type=bondi - COMMAND | VALUE

    DESCRIPTION
        Performs and stores results of Sterimol calculation.

    COMMANDS
        COMMAND is one of the following:

         bury
           Do a Buried Sterimol calculation.

         calculate
           Calculate Sterimol parameters.

         draw_3D
           Draw a 3D representation of the molecule with the Sterimol vectors.

         print_report
           Prints the values of the Sterimol parameters.

         set_points
           Set points for calculation of Sterimol.

         surface_from_radii
           Create surface points from vdW surface.

    VALUES
        VALUE is one of the following:

         B_1

         B_1_value

         B_5

         B_5_value

         L

         L_value

         L_value_uncorrected

         bond_length

We could for example access the ``B_1_value`` attribute or run the method
``print_report`` with the ``verbose=True`` keyword argument.

.. code-block:: console
  :caption: Example of Sterimol CLI

  $ morfeus sterimol tBu.xyz - 1 2 --radii_type=bondi - B_1_value
  2.964748534441907
  $ morfeus sterimol tBu.xyz - 1 2 --radii_type=bondi - print_report
  L         B_1       B_5
  4.31      2.96      3.37

The last command correspond to the following Python code using the module
system.

.. code-block:: python
  :caption: Example of Sterimol with module code

  from morfeus import read_geometry, Sterimol
  elements, coordinates = read_geometry("tBu.xyz)
  sterimol = Sterimol(elements, coordinates, 1, 2, radii_type="bondi")
  sterimol.print_report(verbose=True)

************
Detailed use
************

Positional arguments are passed in sequence separated by spaces, for example:

.. code-block:: console

  $ morfeus sterimol tBu.xyz - 1 2

Keyword arguments are passed with or without an equals sign, so both of these
commands give the same result:

.. code-block:: console

  $ morfeus sterimol tBu.xyz - 1 2 --radii_type=bondi
  $ morfeus sterimol tBu.xyz - 1 2 --radii_type bondi

A single ``-`` is used to indicate that all arguments have been provided and
the function/class should be evaluated. For example, the ``-`` before
``print_report`` in the example below tells Fire to instantiate the Sterimol
class, as we don't want to give more keyword arguments. Then the
``print_report`` method is excuted.

.. code-block:: console

  $ morfeus sterimol tBu.xyz - 1 2 --radii_type bondi - print_report

Arguments following the ``--`` separator go directly to the Fire program. For
example, the ``--`` in the line below makes sure that ``--help`` is sent to
Python Fire instead of the ``Sterimol`` object.

.. code-block:: console

  $ morfeus sterimol tBu.xyz - -- --help

**********************
Partial initialization
**********************

ᴍᴏʀғᴇᴜs uses Fire together with functools.partial__ to partially initialize
classes using the information in the geometry file. When running ``morfeus
sterimol tBu.xyz -``, this is what goes on behind the scenes:

.. code-block:: python
  :caption: cli function for Sterimol

  def cli(file: str) -> Any:
      elements, coordinates = read_geometry(file)
      return functools.partial(Sterimol, elements, coordinates)

The geometry file is read, and the elements and coordinates are used to
partially instantiate a Sterimol object. By specifiying the final hyphen ``-``,
we are telling Fire to evaluate the ``cli`` function and return this partially
instantiated object. We can can then supply the remaining arguments required to
fully instantiate ``Sterimol`` with ``morfeus sterimol tBu.xyz - 1 2 -``.

.. __: https://docs.python.org/3/library/functools.html#functools.partial

********
Chaining
********

A very powerful feature of Fire is chaining, which allows a series of commands
to be run on the same object. The commands (and their arguments) are separated
by hyphens, and the chain should end with a command that gives some output.
Here is one example where we create an electron density isosurface from a cube
file for a dispersion descriptor calculation.

.. code-block:: console

  $ morfeus dispersion corannulene.xyz - --point_surface=False - surface_from_cube corannulene.cub - compute_p_int - print_report
  Surface area (Å²): 248.0
  Surface volume (Å³): 247.8
  P_int (kcal¹ᐟ² mol⁻¹ᐟ²): 25.8

This works because methods that modify
:py:class:`Dispersion <morfeus.dispersion.Dispersion>`
return the object itself:

.. code-block:: python

  def surface_from_cube(self, ...):
      ...
      return self

We can break down the chaining in detail:

``morfeus dispersion corannulene.xyz -``
  Partially instantiate the ``Dispersion`` object

``--point_surface=False -``
  Give the keyword argument ``point_surface=True`` and fully instantiate

``surface_from_cube corannulene.cub -``
  Run ``surface_from_cube`` method with cube filename as argument and return
  same object

``compute_p_int -``
  Run ``compute_p_int`` method and return same object

``print_report``
  Run print_report method and return output

****************
Interactive mode
****************

Another powerful feature of Fire is that it can send the result to an
interactive Python session, where it can be manipulated further. This is
triggered with the Fire argument ``--interactive`` and could be used to,
*e.g.*, access the 3D drawing capabilities of ᴍᴏʀғᴇᴜs.

.. code-block::

  $ morfeus sterimol tBu.xyz - 1 2 - -- --interactive

  Fire is starting a Python REPL with the following objects:
  Modules: fire
  Objects: component, main, morfeus, result, trace

  Python 3.9.2 | packaged by conda-forge | (default, Feb 21 2021, 05:02:20)
  Type 'copyright', 'credits' or 'license' for more information
  IPython 7.22.0 -- An enhanced Interactive Python. Type '?' for help.

  In [1]: result
  Out[1]: Sterimol(14 atoms)

  In [2]: result.bond_length
  Out[2]: 1.1

  In [3]: result.print_report()
  L         B_1       B_5
  4.21      2.86      3.27

  In [4]: result.draw_3D()
  ...

*********
Arguments
*********

Fire attemps to correctly guess the correct type of the arguments given in the
command line. The following style is recommended as it works across different
operating systems. Note that integers being intended as strings need double
quotes.

===== ==================
Type  Recommendation
===== ==================
str   bondi, '"1"'
int   1
float 1.0
list  '[1,2]'
dict  '{"key": "value"}'
bool  True
===== ==================

*********
Resources
*********

More detailed information on using the CLI can be found in the Fire
documentation:

- `The Python Fire guide`__
- `Using the CLI`__
- `Specifying arguments`__

.. __: https://google.github.io/python-fire/guide/
.. __: https://google.github.io/python-fire/using-cli/
.. __: https://google.github.io/python-fire/guide/#argument-parsing

