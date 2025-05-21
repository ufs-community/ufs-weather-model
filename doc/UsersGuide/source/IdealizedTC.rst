.. role:: raw-html(raw)
    :format: html

.. _idealized-tc:

**************************************
Idealized Tropical Cyclone Test Case
**************************************

The idealized, regional tropical cyclone case is derived from the I-HAFS configuration (:cite:t:`Wang2024`) and is designed to support controlled studies of tropical cyclone dynamics and forecast development. This configuration removes real-world data assimilation and ocean coupling, focusing solely on atmospheric forecasts using idealized inputs.

Initial and lateral boundary conditions (ICs/LBCs) are derived from a large-scale, idealized global FV3-based atmospheric forecast. The initial vortex is constructed using the Reed and Jablonowski (2011) method, introducing a weak, balanced storm into an environment favorable for rapid intensification. The preprocessing system generates fixed distributions of geography-related variables and constructs the ICs/LBCs from ``tcvitals`` and GRIB input files.

The configuration mirrors the operational HAFS structure but simplifies terrain and surface properties. It includes:

- Preprocessing to set up the forecast and nest domains  
- Optional vortex initialization  
- FV3-based forecast integration  
- Postprocessing to generate GRIB2 and ATCF output files

A utility called ``cal_vortex`` is available to recalculate wind, temperature, and humidity fields based on user-defined vortex specifications. In a recent experiment, altering damping settings resulted in a stronger, more compact vortex and a rightward track shift after 48 hours of forecast time.

This test case provides a simplified environment to study TC dynamics and forecast behavior. Future development plans include incorporating idealized ocean and wave modules and expanding vortex customization options.

============================
Obtaining Data for HSD Cases
============================

.. include:: ./doc-snippets/hsd_data.rst

.. _run-TC:

=================================================
Running the Idealized Tropical Cyclone Test Case
=================================================

This section explains how to run the baroclinic wave case described above using the ``ufs_test.sh`` script.

Clone the Repository
--------------------

.. include:: ./doc-snippets/clone_hsd.rst

Machine Configuration
----------------------

.. include:: ./doc-snippets/hsd_machine_config.rst

.. _idealized-config:

Test Configuration
-------------------

The idealized TC case can be run as-is without adjusting the configuration. 

Running tests
-------------

.. include:: ./doc-snippets/hsd_run_tests.rst

Example:
^^^^^^^^^

Users with access to the ``epic`` account can run the ``tropical_cyclone`` test case with the ``intel`` compiler on :term:`RDHPCS` where they have access using the following command:

.. code-block:: console

   ./ufs_test.sh -a epic -s -c -k -r -n "tropical_cyclone intel"

Checking Results
-------------------

.. include:: ./doc-snippets/hsd_check_results.rst

For example, to monitor progress or check results for the ``baroclinic_wave`` case, run:

.. code-block:: console

   tail -f ${UFS_WM}/tests-dev/run_dir/tropical_cyclone_intel/err
   tail -f ${UFS_WM}/tests-dev/run_dir/tropical_cyclone_intel/out

.. include:: ./doc-snippets/hsd_notes.rst
