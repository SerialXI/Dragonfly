Quick Start with Experimental Data
==================================

This guide walks through processing real experimental single-particle imaging (SPI) data.
We recommend completing at least one simulated reconstruction first using the
:doc:`simulation-quickstart` guide before proceeding.

.. image:: /_static/images/emc_exp.png
   :alt: Experimental pipeline flowchart
   :width: 60%

Setting Things Up
-----------------

As in the simulation case, create a reconstruction directory to keep things organized::

    dragonfly.init -t spi --legacy

For the AMO 86615 example, start from the sample configuration instead of using
the interactive TUI to build an experimental config from scratch. Copy the
sample config and its companion mask file from ``sample_configs/amo86615/`` into
your reconstruction directory and rename the config to ``config.ini``::

    cp sample_configs/amo86615/config.ini <recon_dir>/config.ini
    cp sample_configs/amo86615/mask_pnccd_back_260_257.byt <recon_dir>/

Data Source
-----------

This tutorial uses data from the AMO end-station of the Linac Coherent Light Source (LCLS),
published as Reddy, et al. *Scientific Data* **4**, 170079 (2017).

Download single hits from `CXIDB <http://cxidb.org/id-58.html>`_. Follow the HDF5 link to
download the data files. Each file contains the dataset ``photonConverter/pnccdBack/photonCount``
with photon converted data from the 4x4 down-sampled pnCCD detector.

Experimental Geometry
---------------------

Add the ``[parameters]`` section to your ``config.ini``::

    [parameters]
    detd = 586
    lambda = 7.75
    detsize = 260 257
    pixsize = 0.3
    stoprad = 40
    ewald_rad = 650.
    polarization = x

For units and other details, see :doc:`config`.
The ``ewald_rad`` parameter is the radius of curvature of the Ewald sphere in voxels.
The 3D grid size is determined by the highest resolution detector pixel distance.
A value of 650 gives a 125x125x125 voxel volume.

Detector File
-------------

Add the ``[make_detector]`` section to create the detector geometry file::

    [make_detector]
    in_mask_file = mask_pnccd_back_260_257.byt
    out_detector_file = data/det_pnccd_back.h5

Generate the detector file by running::

    dragonfly.utils.make_detector

Data Conversion
---------------

Convert HDF5 data to EMC format for each file::

    dragonfly.utils.convert.h5toemc -d photonConverter/pnccdBack/photonCount <HDF5_file>

This creates an EMC file in the ``data`` folder. For other options, run::

    dragonfly.utils.convert.h5toemc -h

The EMC file header contains the total pixel count from the configuration file,
so the ``[parameters]`` section must be complete before converting data.

View your data using the frame viewer::

    dragonfly.frameviewer

Click 'Random' a few times to inspect the patterns.

EMC Parameters
--------------

Add the ``[emc]`` section to ``config.ini``::

    [emc]
    in_photons_list = experiment_files.txt
    in_detector_file = make_detector:::out_detector_file
    output_folder = data/
    log_file = EMC.log
    num_div = 10
    need_scaling = 1
    beta = 0.001
    beta_schedule = 1.41421356 10

Create the photon list file (e.g., ``experiment_files.txt``) with paths to all EMC files::

    data/run_001.emc
    data/run_002.emc
    ...

For this experimental workflow, ``beta = 0.001`` sets a uniform initial
``beta_start`` across frames. The ``beta_schedule`` then specifies the
deterministic annealing schedule where that starting value is effectively scaled
by sqrt(2) every 10 iterations. These aid convergence with high signal data.

Running and Monitoring
----------------------

**Serial execution**::

    dragonfly.emc 100

**MPI parallel execution**::

    mpirun -np 8 dragonfly.emc 100

On a cluster, expect around 800 seconds per iteration at the start,
decreasing to ~240 seconds by iteration 100.

Monitor reconstruction progress with the autoplot GUI::

    dragonfly.autoplot

This displays likelihood values, best volume, and correlation coefficients.

Check Your Work
---------------

After running, verify your reconstruction by viewing the output volumes.
Compare with known results from similar experiments to ensure quality.

See Also
--------

* :doc:`simulation-quickstart` - Simulated data tutorial
* :doc:`faq` - Common troubleshooting questions
