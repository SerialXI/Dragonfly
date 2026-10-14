Quaternion Module
=================

The quaternion module handles orientation sampling for SO(3).

Generating Orientations
-----------------------

For 3D reconstruction, generate uniformly distributed quaternions:

.. code-block:: python

   from dragonfly.quaternion import Quaternion

   quat = Quaternion(num_div=10)

For 2D reconstruction, generate in-plane rotations:

.. code-block:: python

   quat = Quaternion(num_rot=360, point_group='1')  # 360 angles, no symmetry, default
   quat = Quaternion(num_rot=180, point_group='2')  # 180 angles, Friedel symmetry

The sampling in both cases above is 1 degree. The former should be used when the 2D class averages are expected to not be centrosymmetric due to Ewald sphere curvature.

Point Groups
------------

* ``'1'`` - No symmetry
* ``'2'`` - Friedel symmetry (2-fold), for 2D
* ``'A5'`` - Icosahedral symmetry
* ``'S4'`` - Octahedral symmetry

Configuration
-------------

Quaternions can also be generated from a configuration file:

.. code-block:: python

   quat = Quaternion()
   quat.from_config('config.ini', section_name='emc')

.. automodule:: dragonfly.quaternion
   :members:
   :undoc-members:
   :show-inheritance:
