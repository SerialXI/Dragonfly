Progress Viewer (autoplot)
==========================

The Progress Viewer is both a monitor for EMC runs and an interactive analysis workspace for their results. In 2D classification, use it to explore class averages, inspect their measured member frames, examine per-frame metrics, and select subsets for further analysis. This may also be performed without a 3D reconstruction goal for inherently heterogeneous datasets.

In 3D mode it displays volume slices, assignments, and iteration metrics. The image below shows that layout; the walkthrough that follows concentrates on the larger 2D analysis feature set.

.. image:: /_static/images/Screenshot_autoplot.png
   :alt: Progress Viewer screenshot
   :width: 100%

Launch
------

From within a reconstruction directory::

    dragonfly.autoplot

With a specific config file::

    dragonfly.autoplot -c config.ini

To open a particular iteration output, supply the matching config as well::

    dragonfly.autoplot -c config.ini -f data/output_050.h5

The config determines the reconstruction type, class count, input files, and log path. Use the config belonging to the output, especially for 2D results.  Run the GUI on a machine with a graphical display, or use X11 forwarding from a remote system.

A 2D Classification Analysis Walkthrough
----------------------------------------

Follow one run through the questions below: **Have classes emerged? What do they contain? Which frames belong to them? Are the assignments convincing?  Which frames should I keep?**

For a simulated example, first generate data using :doc:`simulation-quickstart`. In that reconstruction directory, save a separate ``config_2d.ini`` with the following section, removing the other config sections::

    [emc]
    in_photons_file = data/photons.emc
    in_detector_file = data/det_sim.h5
    recon_type = 2d
    num_modes = 18
    num_rot = 180
    friedel_sym = 1
    need_scaling = 1
    beta_schedule = 2.0 10
    output_folder = data/2d/
    log_file = logs/EMC_2d.log

Run classification and open its viewer::

    dragonfly.emc -c config_2d.ini -t 4 30
    dragonfly.autoplot -c config_2d.ini

These settings are illustrative; adjust the class count and sampling to your data. A simulated single-particle dataset produces classes representing recurring views, not necessarily different particles or structures.

The screenshots below show an illustrative eight-class run with 90 angular samples and 20 completed iterations, following class 4 through several analysis views. The class count and output paths differ from the starting configuration above.

1. Find the Classes and Track Their Development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The upper display contains a grid of class thumbnails and an enlarged image of the selected class. Click a thumbnail or use the mode slider to select it.  Adjust **VRange** and **Exp** to reveal weak fringes without saturating the central scattering; ``Exp = log`` selects a symmetric logarithmic scale.

Use the iteration slider to compare an early completed iteration with a later one. Look for reproducible features emerging from the initial noise, and for duplicate or poorly supported classes. Class labels are arbitrary across independent runs, so compare their patterns rather than class numbers.

The lower plots show RMS model change, mutual information, likelihood, and the evolution of best-match assignments. White dashed lines mark beta changes; orange lines mark rotation-sampling changes. These plots help diagnose progress, but convergence of a scalar metric alone does not establish scientific validity.  Use **Check** or **Keep checking** to follow an active run.

.. figure:: /_static/images/autoplot_2d_overview.png
   :alt: 2D progress viewer showing eight class thumbnails, an enlarged class 4, iteration metrics, and run controls at iteration 20
   :width: 100%

   The 2D workspace at iteration 20: class thumbnails and the selected class above, iteration metrics below, and display controls and the log on the right.

.. figure:: /_static/images/autoplot_2d_progress.png
   :alt: The same eight-class run displaying class 4 at iteration 2, with less resolved fringes than at iteration 20
   :width: 100%

   Class 4 at iteration 2, using the same display limits as the iteration-20 view above. The fringes are less resolved at this early stage; the metrics panel still shows the full recorded iteration history.

2. Inspect the Measured Frames Behind a Class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Select a populated class, then choose **Analysis > Open Frameviewer**. This opens a masked **comparison-mode** viewer, restricted to frames whose best match is the selected class. It initially uses autoplot's selected iteration.  Selecting a different class in autoplot updates the class-specific frame list.

Browse several frames. Do their measured features agree with the best-match model slice, allowing for photon noise? Check for detector artifacts, residual background, or systematic displacement that a class average could conceal.  The comparison shows a model slice, not a residual map or a complete background-and-exposure-corrected prediction; do not judge agreement solely by the absolute brightness of the two panels.

.. figure:: /_static/images/autoplot_2d_compare.png
   :alt: Comparison frameviewer showing measured frame 8316 on the left and its best-match model slice from class 4 on the right
   :width: 100%

   A measured member of class 4 (left) beside its best-match model slice (right). The viewer retains the frame index and class-specific browsing controls; the central detector mask is visible in both panels.

3. Examine Assignment Confidence and Acquisition Trends
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Choose **Analysis > Frame Metrics Plot** to plot frame-wise values from the selected iteration output. Set **X**, **Y**, and optional **Color** metrics; choose a scatter plot, 2D histogram, or hexbin. Useful starting views include:

* ``frame`` versus ``scale``, colored by ``mode``, to inspect exposure trends.
* ``likelihood`` versus ``mutual_info``, colored by ``mode``, to compare classes and identify frames worth inspecting individually.
* ``occupancy_0`` versus ``occupancy_1`` for a pair of similar classes, to examine how posterior class membership is shared between them.

Use **Current mode only** to restrict the plot to the class selected in autoplot.  The metrics window reloads when the selected iteration changes. Available columns depend on the output: scales require scaling, and occupancy columns require a multi-mode run. Excluded frames are filtered using any saved iteration blacklist and the config's blacklist.

An occupancy is a probability summed over orientations, not just a hard class label. The best joint class/orientation assignment used for member browsing need not be the class with the largest marginal occupancy. Likelihood and mutual information also depend on count level and model assumptions; they are diagnostics, not universal thresholds for good data.

.. figure:: /_static/images/autoplot_2d_metrics.png
   :alt: Frame metrics scatter plot of mutual information versus scale with a logarithmic scale axis and Current mode only enabled
   :width: 85%

   Mutual information versus per-frame scale for 5,344 frames assigned to the selected class, with a logarithmic x axis and **Current mode only** enabled. The trend illustrates why assignment diagnostics should be considered alongside exposure variation.

For a detailed class view, right-click the enlarged class image and select **View detailed**. The pop-out shows the class image and its smoothed occupancy over the sequence of usable hits. The **Filter size** controls the smoothing window; it does not change the inferred occupancies. This is useful for spotting changes through an acquisition, but the hit axis is the retained sequence, not necessarily the original experiment's shot IDs or timestamps.

The zoom settings are retained as you browse between classes, which is useful for examining a smaller region such as the central one, more closely.

.. figure:: /_static/images/autoplot_2d_detail.png
   :alt: Detailed view of class 4 and its marginal occupancy percentage across the hit sequence, smoothed with a filter size of 1000
   :width: 100%

   Class 4 and its marginal occupancy over the hit sequence, smoothed with a filter size of 1,000. The trace is a percentage derived from class probabilities, not a count of hard assignments.

4. Select Classes and Export a Frame Blacklist
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After inspecting classes and their member frames:

#. Choose **Analysis > Mode selection > Toggle mode selection**.
#. Click thumbnails to mark the classes you want to retain. A translucent overlay marks selected classes; click again to deselect.
#. Choose **Analysis > Mode selection > Save blacklist file**.

The file contains one value per frame in concatenated input order: ``0`` for frames whose best joint assignment belongs to a selected class, and ``1`` for the others. It is a hard class-based selection, not an occupancy threshold.  Retain the original input-file order and any mapping to experimental shot IDs.

To analyze that subset in a separate run, add ``blacklist_file`` to a copy of the config and choose a new output folder and log file. A fresh classification can then relearn classes on the selected data without overwriting this run.

.. figure:: /_static/images/autoplot_2d_selection.png
   :alt: Progress viewer in mode-selection state, with translucent overlays on classes 0, 3, and 4 selected for retention
   :width: 100%

   Classes 0, 3, and 4 are marked for retention by translucent overlays. The exported blacklist uses ``0`` for frames whose best match belongs to a selected class and ``1`` for the others.

5. Optional Ways to Explore Class Images
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Display transformations.** Under **Analysis**, subtract the radial minimum to highlight anisotropic features, normalize high-q intensities to compare classes on a common outer-region scale (useful when there is significant unmodelled background), or align their principal axes. These operations change the images held in the viewer, not the saved EMC model or its assignments. **Reparse** reloads the original output. Reparse before quantitative comparisons or subsequent analysis if you want the unmodified class images; keep any transformations explicit in figure captions.

**Common-line PCA.** **Analysis > Open CLPCA** compares radial line-outs between class images and displays a PCA embedding of their common-line correlation matrix. Choose the angular sampling and inner mask radius, then **Calculate CC**. Clicking a point selects the corresponding class in autoplot. The calculation saves a ``CL_CC`` group in the selected output file when writable.  Common-line relationships have a particular interpretation for projections of a common object; an embedding of heterogeneous patterns is not proof that they arise from a single 3D structure.

.. note::

   Please cite Y. Zhuang et al, "Unsupervised learning approaches to characterizing heterogeneous samples using X-ray single-particle imaging", IUCrJ 9(2), 204-214 (2022) if you end up using the CLPCA method to help classify your data.

.. figure:: /_static/images/autoplot_2d_clpca.png
   :alt: Common-line PCA embedding of eight class images, with class 4 marked by a green cross and controls set to 180 angular bins and inner radius 20
   :width: 75%

   Common-line PCA of the eight class images, using 180 angular bins and an inner mask radius of 20. The green cross marks class 4; selecting a point also selects that class in autoplot.

**2D class phasing.** **Analysis > 2D Class Phaser**, also available from the enlarged image's right-click menu, opens a separate phase-retrieval tool. Select a class, set the preprocessing radii and smoothing width, and **Preprocess** before phasing. Set the support size, positivity assumption, algorithm sequence (for example ``50 ER 100 DM 100 ER``), and number of runs, then **Phase**.  Inspect the recovered density, support, and calculated intensity rather than treating a density image as confirmation of a successful classification.  Phase retrieval is optional and requires additional physical assumptions.

**Save** adds the preprocessed intensity, density, and support under ``phasing/`` in the selected iteration file. Work on a copy if you want to keep the original reconstruction output unchanged. This also applies to CLPCA's saved correlation matrix.

.. figure:: /_static/images/autoplot_2d_phasing.png
   :alt: Class 4 in the 2D phaser after phase retrieval, with a recovered-density inset, a 400-pixel support, positivity enabled, and five runs
   :width: 85%

   Optional phase retrieval for class 4, using a 400-pixel support, positivity, and five runs of ``50 ER 100 DM 100 ER``. The inset shows the recovered density, a downstream phase-retrieval result rather than the EMC intensity model itself.

Save Figures and Record Your Choices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use **Image > Save Slices Image** for the class display and **Save Log Plot** for iteration metrics. Pop-out plots have Matplotlib toolbar save controls.  Movie export requires ``ffmpeg`` and an available H.264 encoder. Record the config, iteration, selected classes, display transformations, and mapping to original frame IDs alongside exported figures and blacklists.

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
    In 3D mode, three orthogonal slices through the reconstructed volume at the current layer number. The slice planes can be rotated by Ctrl+dragging or Alt+dragging on a slice panel, and the normal vectors can be set explicitly via right-click.

    In 2D mode, a grid of all class averages is shown on the left with an enlarged view of the selected class on the right. Click a class thumbnail to select it.

**Metrics plots** (bottom left)
    Plots of reconstruction metrics parsed from the log file:

    - **RMS change** of the model between iterations
    - **Mutual information** between frames and orientations
    - **Average log-likelihood**
    - **Orientation convergence** showing the most likely orientation for each frame as a function of iteration (sorted and colored by the last iteration)

    White dashed lines mark changes in ``beta`` (annealing schedule). Orange dashed lines mark changes in the number of rotation samples.

    Use **Analysis > Frame Metrics Plot** for frame-wise scatter plots, 2D histograms, or hexbins from the current HDF5 output file. Available metrics include frame index, ``likelihood``, ``mutual_info``, ``scale`` when present, ``orientation``, ``mode`` for multi-mode runs, ``blacklist``, and per-mode ``occupancy_*`` columns when occupancies are saved.

**Options panel** (right)
    Controls for log file path, volume file path, color scale range, exponent (gamma), iteration/layer/mode sliders, and action buttons.

Controls
--------

**Iteration slider**
    Select which iteration's output to display.

**Layer slider** (3D mode)
    Select the layer number for the orthogonal slices.

**Mode slider** (multi-mode reconstructions)
    Select which mode to display. In 2D mode, clicking a class thumbnail also selects the mode.

**Check / Keep checking**
    ``Check`` reads the log file to detect new completed iterations and updates all plots. ``Keep checking`` polls the log file every 5 seconds for automatic updates during a running reconstruction.

**Plot / Reparse**
    ``Plot`` re-renders with current settings. ``Reparse`` forces re-reading of the volume file from disk.

**VRange / Exp**
    Set the minimum and maximum of the color scale, and the exponent (gamma) for the power-law normalization. Enter ``log`` for symmetric logarithmic normalization.

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
- **Color Map** -- Choose from ``coolwarm``, ``cubehelix``, ``CMRmap``, ``gray``, ``gray_r``, ``jet``.

**Analysis**

- **Open Frameviewer** -- Launch the :doc:`frameviewer` in masked comparison mode, showing members of the current class for multi-mode runs.
- **Subtract radial minimum** -- Remove the radial minimum from the current volume to highlight anisotropic features.
- **Normalize high q** (2D only) -- Normalize the outer region of all classes to 1.
- **Align models** (2D only) -- Align principal axes of all class averages.
- **Mode selection** -- Toggle mode selection for generating a blacklist file.  Click class thumbnails to mark good classes, then save a blacklist that flags frames not belonging to the selected classes.
- **Open CLPCA** -- Launch the Common-Line PCA analysis window.
- **Frame Metrics Plot** -- Launch a pop-out window for frame-wise scatter plots, 2D histograms, or hexbins from the selected HDF5 output file.
- **2D Class Phaser** (2D only) -- Launch 2D phasing of class averages.

Right-click the enlarged 2D class image for **View detailed** (image and occupancy trace) and **2D Phaser**. For 3D results, right-click a slice to set its normal vector, or use Ctrl/Alt-dragging to rotate the slice plane.
