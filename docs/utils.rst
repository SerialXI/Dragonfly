Utilities
==========

The utilities directory contains Python utilities, command-line tools, and Cython
extensions for data processing and simulation.

Command-line Tools
------------------

The package provides several command-line tools:

* ``dragonfly.init`` - Initialize a reconstruction with an interactive terminal wizard powered by ``prompt_toolkit`` and ``rich``; supports ``--copy-aux`` for self-contained setups

**Utilities to simulate data:**

* ``dragonfly.utils.make_densities`` - Generate density maps
* ``dragonfly.utils.make_intensities`` - Generate intensity maps
* ``dragonfly.utils.make_detector`` - Create detector files
* ``dragonfly.utils.make_data`` - Generate simulated data
* ``dragonfly.utils.sim_setup`` - Wrapper to run all utils above

**Graphical programs:**

* ``dragonfly.autoplot`` - Progress viewer (:doc:`user-guides/autoplot`)
* ``dragonfly.frameviewer`` - Frame viewer (:doc:`user-guides/frameviewer`)

Simulation Utilities
--------------------

.. automodule:: dragonfly.utils.sim_setup
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: dragonfly.utils.make_densities
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: dragonfly.utils.make_intensities

.. automodule:: dragonfly.utils.make_detector
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: dragonfly.utils.make_data
   :members:
   :undoc-members:
   :show-inheritance:

Configuration Parser
--------------------

.. automodule:: dragonfly.utils.py_src.read_config
   :members:
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 1

   utils/py_utils
   utils/slices
   utils/calc_cc
