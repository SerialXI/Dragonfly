The EMC Algorithm
=================

Overview
--------

Dragonfly uses the Expand-Maximize-Compress (EMC) algorithm to estimate models from noisy diffraction patterns with unknown assignments. Two common uses are:

* **3D reconstruction:** estimate a reciprocal-space intensity volume from patterns expected to arise from a common object in unknown orientations.  Model predictions sample that volume on the detector's Ewald surface, including its curvature.
* **2D classification:** learn multiple class intensity images and estimate each frame's class and in-plane orientation. This is useful for exploring serial diffractive imaging data even when a single reproducible 3D object is not an appropriate description or a 3D reconstruction is not the goal.

Algorithm Principle
-------------------

The inference problem has two coupled unknowns:

#. the intensity model or models
#. the orientation and, for multiple models, class of every measured frame

EMC addresses this by alternating between estimating how well each frame matches candidate model/orientation combinations and updating models from those weighted matches. No orientation or class labels are required in advance. The photon-counting path uses a Poisson likelihood; optional per-frame scale factors account for intensity variation, and sparse photon data can include an explicit background model.

The Three Steps
---------------

Each EMC iteration consists of:

**1. Expand**

Sample each current model on the orientation grid. For 3D, generate detector views from the volume; for 2D, generate in-plane rotated views of each class.  Compare each frame's observed counts with the candidate predictions.

**2. Maximize**

Convert comparison scores into joint model/orientation probabilities and accumulate weighted contributions from the data. Frames can contribute to multiple candidates rather than being assigned irrevocably to one class or orientation. Summing those probabilities over orientations gives a frame's marginal class occupancies.

**3. Compress**

Merge the orientation-specific estimates back onto the regular model grid: a volume for 3D, or a set of images for 2D. The next iteration starts from this updated representation.

What Dragonfly Reconstructs
---------------------------

EMC estimates intensities, not a real-space electron density. Its output is a 3D intensity volume or a stack of 2D class intensities, plus per-frame assignments and diagnostics. Phase retrieval is a separate downstream step; Dragonfly's viewer includes an optional 2D class phasing tool.

2D classes describe recurring patterns, not necessarily different structures.  They can capture different views, sample variation, backgrounds, or artifacts.  Use them to inspect, average, and select data, and validate their physical interpretation against measured frames and experimental metadata. A useful classification does not establish that the whole dataset supports a 3D model.

Orientation Sampling
--------------------

3D orientations are sampled on a quaternion grid controlled primarily by ``num_div``. A larger subdivision count means finer sampling and more orientations to test per iteration.

For 3D reconstructions, the number of sampled orientations is:

.. math::

    N_{\mathrm{orient}} = 10\,(5n^3 + n),

where :math:`n` is ``num_div``, before any point-group symmetry reduction.

This count grows quickly, so increasing ``num_div`` improves angular resolution at a significant computational cost.

Coarse-to-fine refinement avoids a full global search at the fine sampling by restricting each frame to neighborhoods of its likely coarse orientations. Its two-value notation is ``num_div = fine coarse``; see :ref:`orientation-refinement` for the workflow, the role of saved probabilities, and its current implementation status in the Python package.

For 2D, ``num_rot`` sets the number of in-plane angles and ``num_modes`` sets the number of classes. Angles span 360 degrees by default, or 180 degrees with ``friedel_sym = 1``. Centrosymmetry is not automatic for an arbitrary measured pattern, especially with appreciable Ewald curvature; impose it only when justified. See :doc:`../user-guides/config` for example configurations.

Deterministic Annealing
-----------------------

Dragonfly also supports deterministic annealing through the ``beta`` and ``beta_schedule`` settings. Early iterations can use broader orientation probabilities; later iterations gradually sharpen them. This often helps the reconstruction or classification avoid poor local optima, especially for cleaner data. ``beta = auto`` estimates starting values per frame; ``beta_factor`` adjusts them, and effective beta is capped at one.

Convergence and Practical Use
-----------------------------

In practice, users monitor quantities such as log-likelihood, RMS model change, and the visual stability of the volume or class images. EMC is iterative and its behavior depends strongly on data quality, detector geometry, masking, angular sampling, and annealing settings.

Inspect member frames and assignment uncertainty, and compare independent starts or data subsets. A stable model can still describe a detector artifact or an incorrect physical assumption. The :doc:`../user-guides/autoplot` walkthrough shows these diagnostics for 2D classification.

The algorithm is designed for large sets of sparse diffraction patterns, which is why Dragonfly parallelizes the work with MPI and OpenMP.

References
----------

* Loh, D. M. & Elser, V. (2009). Phys. Rev. E 80, 026705.
* Ayyer et al. (2016). J. Appl. Cryst. 49, 1320-1335.
