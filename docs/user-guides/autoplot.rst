Progress Viewer (autoplot)
==========================

The Progress Viewer is a GUI for monitoring EMC reconstruction progress.
It displays volume slices (3D mode) or class averages (2D mode) alongside
plots of reconstruction metrics as a function of iteration.

.. image:: /_static/images/Screenshot_autoplot.png
   :alt: Progress Viewer screenshot
   :width: 100%

Launch
------

From within a reconstruction directory::

    dragonfly.autoplot

With a specific config file::

    dragonfly.autoplot -c config.ini

To view slices of a particular volume file without tracking a running reconstruction::

    dragonfly.autoplot -f data/output_050.h5

Command-line Options
--------------------

.. code-block:: text

    usage: dragonfly.autoplot [-c CONFIG_FILE] [-f VOLUME_FILE]

    -c, --config_file   Path to config file (default: config.ini)
    -f, --volume_file   Show slices of a particular file instead of tracking output

Display Areas
-------------

The window is split into three panels:

**Volume slices** (top left)
    In 3D mode, three orthogonal slices through the reconstructed volume at the
    current layer number. The slice planes can be rotated by Ctrl+dragging or
    Alt+dragging on a slice panel, and the normal vectors can be set explicitly
    via right-click.

    In 2D mode, a grid of all class averages is shown on the left with an
    enlarged view of the selected class on the right. Click a class thumbnail
    to select it.

**Metrics plots** (bottom left)
    Plots of reconstruction metrics parsed from the log file:

    - **RMS change** of the 3D volume between iterations
    - **Mutual information** between frames and orientations
    - **Average log-likelihood**
    - **Orientation convergence** showing the most likely orientation for each
      frame as a function of iteration (sorted and colored by the last iteration)

    White dashed lines mark changes in ``beta`` (annealing schedule).
    Orange dashed lines mark changes in the number of rotation samples.

    Use **Analysis > Frame Metrics Plot** for frame-wise scatter plots, 2D
    histograms, or hexbins from the current HDF5 output file. Available metrics
    include frame index, ``likelihood``, ``mutual_info``, ``scale`` when present,
    ``orientation``, ``mode`` for multi-mode runs, ``blacklist``, and per-mode
    ``occupancy_*`` columns when occupancies are saved.

**Options panel** (right)
    Controls for log file path, volume file path, color scale range, exponent
    (gamma), iteration/layer/mode sliders, and action buttons.

Controls
--------

**Iteration slider**
    Select which iteration's output to display.

**Layer slider** (3D mode)
    Select the layer number for the orthogonal slices.

**Mode slider** (multi-mode reconstructions)
    Select which mode to display. In 2D mode, clicking a class thumbnail
    also selects the mode.

**Check / Keep checking**
    ``Check`` reads the log file to detect new completed iterations and updates
    all plots. ``Keep checking`` polls the log file every 5 seconds for
    automatic updates during a running reconstruction.

**Plot / Reparse**
    ``Plot`` re-renders with current settings. ``Reparse`` forces re-reading of
    the volume file from disk.

**VRange / Exp**
    Set the minimum and maximum of the color scale, and the exponent (gamma)
    for the power-law normalization. Enter ``log`` for symmetric logarithmic
    normalization.

Keyboard Shortcuts
------------------

===================  ==================================
Shortcut             Action
===================  ==================================
``Enter``            Plot / refresh
``Ctrl+S``           Save slices image
``Ctrl+K``           Check for new iterations
``Ctrl+Q``           Quit
===================  ==================================

Menu Reference
--------------

**File**

- **Load Volume** -- Open a volume file (HDF5 or binary) for viewing.
- **Load Config** -- Reload the GUI with a different config file.
- **Quit** -- Close the application.

**Image**

- **Save Slices Image** -- Save the current volume slices plot as a PNG.
- **Save Log Plot** -- Save the metrics panel as a PNG.
- **Save Layer Movie** -- Save an MP4 animation sweeping through all layers.
- **Save Iteration Movie** -- Save an MP4 animation sweeping through all iterations.
- **Color Map** -- Choose from ``coolwarm``, ``cubehelix``, ``CMRmap``, ``gray``,
  ``gray_r``, ``jet``.

**Analysis**

- **Open Frameviewer** -- Launch the :doc:`frameviewer` showing frames assigned to
  the current mode (in multi-mode reconstructions).
- **Subtract radial minimum** -- Remove the radial minimum from the current volume
  to highlight anisotropic features.
- **Normalize high q** (2D only) -- Normalize the outer region of all classes to 1.
- **Align models** (2D only) -- Align principal axes of all class averages.
- **Mode selection** -- Toggle mode selection for generating a blacklist file.
  Click class thumbnails to mark good classes, then save a blacklist that flags
  frames not belonging to the selected classes.
- **Open CLPCA** -- Launch the Common-Line PCA analysis window.
- **Frame Metrics Plot** -- Launch a pop-out window for frame-wise scatter plots,
  2D histograms, or hexbins from the selected HDF5 output file.
- **2D Class Phaser** (2D only) -- Launch 2D phasing of class averages.
