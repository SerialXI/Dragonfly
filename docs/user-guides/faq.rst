Frequently Asked Questions
==========================

Installation
------------

**How do I install dragonfly on an HPC cluster?**

Load your MPI module first, then install::

    module load mpi/openmpi-4.1
    python -m pip install mpi4py
    python -m pip install dragonfly-spi

Also load the site's compiler, HDF5, and GSL modules if needed. See :doc:`installation` for build prerequisites.

**Why does mpi4py fail to install or import?**

``mpi4py`` needs an available MPI runtime and, when built from source, its development headers and compiler wrapper. It must use an MPI library compatible with Dragonfly and the ``mpirun`` used to launch it. Loading the intended MPI module before installation helps select the correct environment, but matching the linked libraries matters more than installation order alone.

Check the runtime with::

    python -c "from mpi4py import MPI; print(MPI.Get_library_version())"

**Can I use conda instead of pip?**

Yes, but keep a consistent compiler and MPI environment. For a conda-managed OpenMPI runtime, one possible starting point is::

    conda install mpi4py openmpi

This does not replace the other Dragonfly build prerequisites. Avoid mixing that MPI runtime with an incompatible cluster MPI module.

Data Processing
---------------

**What file formats does dragonfly support?**

* Binary photon files (``.emc``)
* Sparse EMC-HDF5 photon files (usually ``.h5``)
* Detector files (``.h5``)
* PDB files for simulations

Binary EMC includes sparse counts, dense integer counts, and dense double data.  ``EMCReader`` can also browse dense datasets in other HDF5 layouts; those layouts need conversion before reconstruction. See :doc:`emc-format` and :doc:`detector-format`.

**How do I convert my data to dragonfly format?**

Read frames with ``h5py`` and write them with ``EMCWriter``. The :doc:`experimental-quickstart` tutorial includes a short streaming conversion script. Ensure that photon values, pixel ordering, and detector geometry match.

**Can I analyze data without doing a 3D SPI reconstruction?**

Yes. Use frame inspection, powder sums, and 2D classification to explore serial diffractive imaging data. The sample need not be described by a single reproducible 3D object to make those tools useful. In a 2D run, ``num_modes`` sets the number of classes and ``num_rot`` sets in-plane rotation sampling.  See :doc:`config` for an example and :doc:`autoplot` for analysis tools.

**Are 2D classes distinct structures?**

Not necessarily. Classes describe recurring diffraction patterns and can reflect orientation, sample heterogeneity, backgrounds, or detector artifacts.  Inspect individual member frames and class occupancies, and compare independent starts before assigning physical meaning to them.

**Can I use background-subtracted or analog detector intensities as photon counts?**

Not by silently rounding or clipping them. Sparse EMC stores nonnegative integer photon counts and uses a Poisson likelihood. Dense double EMC has a separate Gaussian path with a required positive ``gaussian_sigma`` and without per-frame scaling in its likelihood. The appropriate noise model depends on your calibration and preprocessing; see :doc:`config`.

**How do I choose the orientation sampling density (num_div)?**

This is a function of your dimensionless resolution, counted as the ratio of the size of your particle to your resolution, which is also around twice the highest order speckle/fringe.

For 3D without point-group reduction, the total number of orientations is ``10*(5*num_div^3 + num_div)``:

========  ============
num_div   Orientations
========  ============
2         420
4         3240
6         10860
8         25680
========  ============

The average angle between neighboring orientations is around ``0.91/num_div`` radians. A rule of thumb is to use ``num_div`` sampling for reconstruction up to a radius of ``num_div`` speckles/fringes.

For 2D classification, use ``num_rot`` instead. With ``friedel_sym = 0``, ``num_rot`` samples a full 360-degree turn; with ``friedel_sym = 1``, it samples 180 degrees. Do not impose that symmetry unless it is appropriate for the data.

Reconstruction
--------------

**How do I resume an interrupted reconstruction?**

Use the ``-r`` flag::

    dragonfly.emc -c config.ini -r 20

This resumes for 20 more iterations using the last completed iteration in the configured log and its ``output_XXX.h5`` file. Keep both files. If the log is missing, the command starts a new run instead, so check the startup messages.

**What do the output volumes represent?**

For 3D runs, volumes contain the intensity distribution of the sample in reciprocal space. The voxel size can be calculated as::

    1 / (lambda * ewald_rad)

where ``ewald_rad`` is the Ewald sphere radius in voxels, specified in the detector file.

For 2D runs, ``intens`` contains a stack of learned class intensities, not a 3D volume or a real-space density. See :doc:`autoplot` for browsing classes.

**The GUI doesn't display properly.**

The autoplot GUI requires a display. On remote systems:

* Use X11 forwarding: ``ssh -X``

Check that ``$DISPLAY`` is set and that your local machine has an X server.  The computation itself does not require a graphical display.

**How do I visualize the reconstructed volumes?**

Dragonfly outputs h5 files for every iteration. Use:

* The autoplot GUI
* Custom Python code to read the ``intens`` dataset. It has shape ``(num_modes, size, size, size)`` for 3D and ``(num_modes, size, size)`` for 2D.

Performance
-----------

**How many cores should I use?**

Optimal configuration depends on your system:

* Total cores = MPI processes × OpenMP threads

For a 64-core node::

    mpirun -np 8 dragonfly.emc -t 8 20

This runs 20 iterations with eight processes and eight threads per process.  Use a scheduler allocation that matches the requested total, and specify ``-t`` explicitly to avoid oversubscribing a node.

**Why are the earlier iterations slower?**

The first iterations typically have broader probability distributions, leading to slower updated tomogram calculations. At the end, frames are often assigned to single orientations, which makes the ``update`` step faster.

**Can I run on GPUs?**

The reconstruction implementation in this package uses CPUs through MPI and OpenMP. These documentation commands do not enable GPU computation.

Getting Help
------------

For additional support:

* GitHub Issues: https://github.com/SerialXI/Dragonfly/issues
* Documentation: https://dragonfly-spi.readthedocs.io
