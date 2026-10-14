Welcome to Dragonfly's documentation!
=====================================

Dragonfly is an analysis toolkit for **serial diffractive imaging** data: large collections of diffraction patterns measured one exposure at a time. It combines tools for inspecting individual patterns, learning recurring pattern types, forming class averages, selecting frames, and reconstructing reciprocal-space intensity with the EMC (Expand-Maximize-Compress) algorithm.

You do not need to target a 3D reconstruction to use Dragonfly. **2D classification** can help explore a heterogeneous dataset, improve the signal-to-noise ratio of recurring patterns through rotationally aligned averaging, and identify subsets for further analysis. Class averages are statistical descriptions of the data, not automatically distinct structures.

For samples that can be described by a common 3D intensity distribution, Dragonfly can also recover that distribution from patterns with unknown orientations. Reconstruction and classification use MPI and OpenMP for CPU-based parallel computation. A simulator generates photon-counting data from PDB structures for learning the workflow and testing analysis choices.

Getting Started
---------------

Start with :doc:`user-guides/installation`, then choose the tools that match your question:

* **Inspect your data:** use the :doc:`user-guides/frameviewer` to browse frames and view their powder sum. The :doc:`user-guides/emc-format` and :doc:`user-guides/detector-format` describe the required data and geometry.
* **Explore and classify patterns:** the :doc:`user-guides/config` includes a starting configuration for 2D classification. The :doc:`user-guides/autoplot` walks through class inspection, measured/model comparisons, per-frame metrics, and class-based frame selection.
* **Reconstruct in 3D:** follow :doc:`user-guides/experimental-quickstart` for a worked experimental SPI example. This requires stronger assumptions about the sample than exploratory 2D classification.
* **Try simulated data:** follow :doc:`user-guides/simulation-quickstart`. Simulations are useful for learning the tools and exploring feasibility limits.

The :doc:`background/emc-algorithm` explains how the models and frame assignments are estimated. The API reference below provides Python interfaces for integrating Dragonfly into your own analysis workflows.

.. toctree::
   :maxdepth: 1
   :caption: User Guides

   user-guides/installation
   user-guides/simulation-quickstart
   user-guides/simulator
   user-guides/config
   user-guides/experimental-quickstart
   user-guides/emc-format
   user-guides/detector-format
   user-guides/autoplot
   user-guides/frameviewer
   user-guides/faq

.. toctree::
   :maxdepth: 1
   :caption: Background

   background/emc-algorithm
   background/oversampling

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
