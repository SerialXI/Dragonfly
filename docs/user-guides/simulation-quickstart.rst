Quick Start with Simulations
============================

This guide walks through the simulation workflow for single-particle imaging.

.. image:: /_static/images/emc_sim.png
   :alt: Simulation pipeline flowchart
   :width: 60%

For installation instructions, see :doc:`installation`.

Spawn a Reconstruction Directory
--------------------------------

Go to a directory with plenty of space and start a reconstruction session::

    dragonfly.init

This creates ``recon_0001`` with ``data/``, ``images/``, and ``logs/`` directories, a template ``config.ini``, and an ``aux`` link to the installed auxiliary data. If that directory already exists, the next available number is used. The command-line tools were installed with the package; initialization does not compile or copy executables.

For customization options::

    dragonfly.init -h

Configure Experiment
--------------------

Go to your newly created reconstruction directory::

    cd recon_0001

Edit ``config.ini`` to customize parameters:

* ``pdb_code``: PDB identifier, or replace it with ``in_pdb_file`` for a local PDB
* Detector setup (in ``[parameters]`` section):

  * ``detd``: Detector distance (mm)
  * ``lambda``: Photon wavelength (Angstrom)
  * ``detsize``: Linear detector size (pixels)
  * ``pixsize``: Pixel size (mm)
  * ``stoprad``: Beamstop radius (pixels)

* ``num_data``: Number of diffraction patterns
* ``fluence``: Incident beam intensity (photons/um^2)
* ``log_file``: Name of log file

Generate simulated data::

    dragonfly.utils.sim_setup

More details in :doc:`simulator` for a description of the constituent utilities.

Run EMC Reconstruction
----------------------

**On a single node (only OpenMP)**::

    dragonfly.emc -t <num_threads> 10

Starts reconstruction with default ``config.ini`` config file for 10 iterations using `<num_threads>` threads.

**Continue for 20 more iterations**::

    dragonfly.emc -t 4 -r 20

Resuming requires the log and the output file from the last completed iteration. Without ``-r``, a new run starts and will overwrite existing outputs.

**MPI reconstruction using different config file**::

    mpirun -np 4 dragonfly.emc -c config.ini -t 4 50

This uses four MPI processes with four OpenMP threads each (16 CPU cores).

Often these larger reconstructions are submitted through a job scheduler like SLURM, PBS etc.

Monitor Progress
----------------

In the reconstruction directory::

    dragonfly.autoplot -c <config_file>

This GUI monitors reconstruction progress and can automatically look for new reconstructed volumes. See :doc:`autoplot` for full documentation.

.. image:: /_static/images/Screenshot_for_quickstart.png
   :alt: Progress Viewer after simulation reconstruction
   :width: 100%

Example Workflow
----------------

Complete example with 8 cores::

    dragonfly.init
    cd recon_0001
    # Edit config.ini if needed before generating data.
    dragonfly.utils.sim_setup
    mpirun -np 2 dragonfly.emc -t 4 30
    dragonfly.autoplot

Additional Resources
--------------------

* See :doc:`config` for current 2D and 3D configuration examples
* Check the configured ``log_file`` (``logs/EMC.log`` in the template) and
  ``simdata.log`` for reconstruction and simulation progress
* Refer to :doc:`config` for description of configuration parameters
* Refer to :doc:`simulator` for data stream simulator details
