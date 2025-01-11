The HSD cases are configured to be run on NOAA Tier-1 platforms, and the configuration files for each platform are located at:

.. code-block:: console

   ${UFS_WM}/tests-dev/machine_config/machine_<PLATFORM>.config

where ``<PLATFORM>`` corresponds to the name of the platform. These configuration files load the necessary Python and Rocoto modules for each platform. Users generally do not need to make any changes to these files. 