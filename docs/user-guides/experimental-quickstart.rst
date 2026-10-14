3D Reconstruction from Experimental Data
========================================

This tutorial walks through a 3D reconstruction using publicly available single-particle imaging (SPI) data. It illustrates a particular sample and detector, rather than prescribing geometry or annealing settings for every experiment. Install Dragonfly first using :doc:`installation`; running a simulation beforehand is recommended, but optional.

3D reconstruction assumes that the selected patterns can be described by a common reciprocal-space intensity distribution in unknown orientations. If your aim is to explore heterogeneous serial diffraction data without that assumption, start with frame inspection or 2D classification instead; see :doc:`config` and :doc:`autoplot`.

.. image:: /_static/images/emc_exp.png
   :alt: Experimental pipeline flowchart
   :width: 60%

Setting Things Up
-----------------

Create a reconstruction directory to keep inputs, outputs, and logs organized::

    dragonfly.init -t spi
    cd spi_0001

Initialization creates ``data/``, ``images/``, ``logs/``, a template config, and an ``aux`` link to the installed auxiliary data. If ``spi_0001`` already exists, use the next available directory reported by the command. Executables are installed with the package, not compiled during initialization.

Data Source
-----------

This tutorial uses data from the AMO end-station of the Linac Coherent Light Source (LCLS), published as Reddy, et al. *Scientific Data* **4**, 170079 (2017).

Download single hits from `CXIDB <http://cxidb.org/id-58.html>`_. Follow the HDF5 link to download the data files. Each file contains the dataset ``photonConverter/pnccdBack/photonCount`` with photon converted data from the assembled, 4x4 binned pnCCD detector.

Experimental Geometry
---------------------

Replace the template's ``[parameters]`` section in ``config.ini`` with::

    [parameters]
    detd = 586
    lambda = 7.75
    detsize = 260 257
    pixsize = 0.3
    stoprad = 40
    ewald_rad = 650.
    polarization = x

For units and other details, see :doc:`config`.

The ``ewald_rad`` parameter is the radius of curvature of the Ewald sphere in voxels.  The 3D grid size is determined from the largest reciprocal-space radius of the unmasked detector pixels. It depends on the geometry and mask, not on ``ewald_rad`` alone. See :doc:`detector-format` for voxel units and conversion to physical reciprocal-space coordinates.

Detector File
-------------

Replace the template's ``[make_detector]`` section with::

    [make_detector]
    in_mask_file = aux/mask_pnccd_back_260_257.byt
    out_detector_file = data/det_pnccd_back.h5

Generate the detector file by running::

    dragonfly.utils.make_detector

See :doc:`detector-format` for the HDF5 geometry fields, units, and pixel masks.  The mask above is distributed with Dragonfly's auxiliary data. It is specific to this experiment.

Data Conversion
---------------

Save the following as ``convert_photons.py`` in the reconstruction directory.  It reads one frame at a time and writes sparse binary EMC data without loading the entire HDF5 dataset into memory:

.. code-block:: python

    import sys

    import h5py
    import numpy as np
    from dragonfly import Detector, EMCWriter

    input_fname, output_fname = sys.argv[1:]
    det = Detector('data/det_pnccd_back.h5')
    dataset_name = 'photonConverter/pnccdBack/photonCount'

    with h5py.File(input_fname, 'r') as source:
        frames = source[dataset_name]
        with EMCWriter(output_fname, det.num_pix, hdf5=False) as wemc:
            for frame_num in range(frames.shape[0]):
                frame = frames[frame_num].ravel()
                wemc.write_frame(frame.astype('i4'))

Run it for each downloaded file, using a different output name each time::

    python convert_photons.py amo86615_182_PR772_single.h5 data/photons_182.emc

See :doc:`emc-format` for the photon data layout and header fields.

Inspect the first converted file before configuring reconstruction::

    dragonfly.frameviewer -e data/photons_182.emc -d data/det_pnccd_back.h5 -M

Click 'Random' to inspect several patterns. Check the beam center, mask, detector layout, and count values. A powder sum is also useful::

    dragonfly.frameviewer -e data/run_001.emc -d data/det_pnccd_back.h5 -M -P

EMC Parameters
--------------

Replace the template's entire ``[emc]`` section with the following. In particular, remove its simulated ``in_photons_file`` entry rather than leaving both file and list inputs enabled::

    [emc]
    in_photons_list = experiment_files.txt
    in_detector_file = make_detector:::out_detector_file
    recon_type = 3d
    output_folder = data/
    log_file = EMC.log
    num_div = 10
    need_scaling = 1
    beta = auto
    beta_schedule = 1.41421356 10

Create ``experiment_files.txt`` with one photon-file path per line, in the
order you want frames concatenated. Prefer absolute paths so that the same list
works in both reconstruction and GUI tools::

    /absolute/path/to/spi_0001/data/run_001.emc
    /absolute/path/to/spi_0001/data/run_002.emc

``beta = auto`` estimates a starting value per frame. The schedule multiplies those starting values by approximately sqrt(2) every 10 iterations, up to a maximum effective beta of one. This controls the sharpness of orientation probabilities. 
Running and Monitoring
----------------------

First test the input configuration with a short run on one MPI process and four OpenMP threads::

    dragonfly.emc -t 4 5

Check ``EMC.log`` and ``data/output_005.h5`` before committing to a long run. If setup and early outputs look sensible, continue for another 95 iterations::

    dragonfly.emc -t 4 -r 95

Alternatively, resume with MPI parallelism (eight processes, four threads each)::

    mpirun -np 8 dragonfly.emc -t 4 -r 95

Choose one continuation command, and allocate the corresponding CPU cores.  Runtime depends on hardware, data volume, orientation sampling, and probability sparsity; measure it on your own setup. Keep the log and iteration output together for resuming. A fresh run without ``-r`` can overwrite existing results.

Monitor reconstruction progress with the autoplot GUI::

    dragonfly.autoplot

This displays model slices, RMS model change, mutual information, likelihood, and best-match assignments across iterations. See :doc:`autoplot` for controls and measured/model comparisons.

Check Your Work
---------------

Inspect model slices and compare measured frames against their best-match model slices. Look for persistent fringes and physically plausible structure rather than interpreting a decreasing RMS change alone as proof of correctness.  Check stability across independent starts and, where possible, independent data subsets. Results are reciprocal-space intensities, not electron density; phase retrieval is a separate downstream step.

See Also
--------

* :doc:`simulation-quickstart` - Simulated data tutorial
* :doc:`faq` - Common troubleshooting questions
