==============
Notes on usage
==============

***************
Without display
***************

Matplotlib may cause problems on computers without displays, e.g., nodes on
computer clusters. This can be solved by changing the plotting backend to "agg".
This could be done at the top of the Python script using morfeus:

.. code-block:: python

  import matplotlib
  matplotlib.use('Agg')

Alternatively, set the environment variable in the shell before launching the
script (Linux example):

.. code-block:: console

  export MPLBACKEND="agg"