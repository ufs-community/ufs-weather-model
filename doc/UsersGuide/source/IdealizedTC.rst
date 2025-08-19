.. role:: raw-html(raw)
    :format: html

.. _idealized-tc:

**************************************
Idealized, Regional Tropical Cyclone Case
**************************************

The idealized, regional tropical cyclone case is derived from the I-HAFS configuration (:cite:t:`Wang2024`) and is designed to support controlled studies of tropical cyclone dynamics and forecast development. It uses the ``FV3_HAFS_v1_thompson_nonsst`` physics suite, which is documented in the `DTC UFS HAFS v1 Scientific Documentation <https://dtcenter.ucar.edu/GMTB/UFS_SRW_HSD_TC/scidoc/_h_a_f_sv1_page.html>`_. The configuration used in this case removes real-world data assimilation and ocean coupling, focusing solely on atmospheric forecasts using idealized inputs.

The case is configured to run at 4-km resolution, with 81 vertical levels. The forecast is initialized on 24 August 2019, and initial/lateral boundary conditions are provided for up to a five-day forecast duration. Initial and lateral boundary conditions (ICs/LBCs) are derived from a large-scale, idealized global FV3-based atmospheric forecast. The initial vortex is constructed using the Reed and Jablonowski (2011) method, introducing a weak, balanced storm into an environment favorable for rapid intensification. 

This idealized test case uses components derived from the I-HAFS configuration, but within the UFS HSD framework the workflow is simplified. The following capabilities are available to users running this case:

- FV3-based forecast integration 
- Adjustable physics suites, namelist settings, and computational parameters
- Optional vortex initialization is present in I-HAFS, but not invoked in this test case
- Postprocessing and preprocessing steps (e.g., IC/LBC generation, terrain setup) are handled outside of this test case and are not included in the ``ufs-weather-model`` workflow

.. note::

   While the I-HAFS system includes preprocessing to set up the forecast and nest domains and postprocessing to generate GRIB2 and ATCF output files, this UFS HSD test case **does not** perform those steps. It relies on pre-generated ICs/LBCs, which are provided as part of the test data.

The script below generates 10-m wind plots from the model's GRIB output and can create an animated GIF to visualize the tropical cyclone's evolution.

.. code-block:: console

   ufs-weather-model/tests-dev/test_cases/utils/plot_tc.sh 

This test case provides a simplified environment to study TC dynamics and forecast behavior. Future development plans include incorporating idealized ocean and wave modules and expanding vortex customization options.

============================
Obtaining Data for HSD Cases
============================

.. include:: ./doc-snippets/hsd_data.rst

.. _run-TC:

=================================================
Running the Idealized Tropical Cyclone Test Case
=================================================

This section explains how to run the Idealized Tropical Cyclone case described above using the ``ufs_test.sh`` script.

Clone the Repository
--------------------

.. include:: ./doc-snippets/clone_hsd.rst

Machine Configuration
----------------------

.. include:: ./doc-snippets/hsd_machine_config.rst

.. _idealized-config:

Test Configuration
-------------------

By default, the forecast length and runtime settings for this idealized tropical cyclone test case are conservative and may need adjustment to simulate a complete tropical cyclone lifecycle.

In the file:

``ufs-weather-model/tests-dev/test_cases/tests/tropical_cyclone``

the following variables can be modified:

.. code-block:: console

   FH_MAX=3

Change to:

.. code-block:: console

   FH_MAX=120

This sets the forecast length to 120 hours (approximately 5 days), which matches the length supported by the provided IC/LBC data.

Also, the wallclock time limit is set as:

.. code-block:: console

   WLCLK=00:30

Change to something like:

.. code-block:: console

   WLCLK=08:00

This allows enough time (6–8 hours recommended) for the full 120-hour simulation to run, depending on system performance.

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

For example, to monitor progress or check results for the ``tropical_cyclone`` case, run:

.. code-block:: console

   tail -f ${UFS_WM}/tests-dev/run_dir/tropical_cyclone_intel/err
   tail -f ${UFS_WM}/tests-dev/run_dir/tropical_cyclone_intel/out

.. include:: ./doc-snippets/hsd_notes.rst
