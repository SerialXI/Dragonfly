Parameters Module
=================

The parameters module stores reconstruction parameters.

Configuration
-------------

Load parameters from a configuration file, which should be the most common way to set the various parameters:

.. code-block:: python

   from dragonfly.params import EMCParams

   param = EMCParams()
   param.from_config('config.ini', section_name='emc')

.. automodule:: dragonfly.params
   :members:
   :undoc-members:
   :show-inheritance:

Reconstruction Types
--------------------

* ``RECON3D`` - 3D volume reconstruction
* ``RECON2D`` - 2D slice reconstruction
* ``RECONRZ`` - Cylindrical (R-Z) reconstruction
