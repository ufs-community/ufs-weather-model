.. role:: raw-html(raw)
    :format: html

.. _baroclinic-wave:

*****************************
Baroclinic Instability Case
*****************************

The UFS WM baroclinic wave case adapts the test outlined in :cite:t:`Jablonowski&Williamson2006` (2006). This test is designed to evaluate the accuracy of various atmospheric models in simulating a baroclinic wave, which commonly forms in the Northern Hemisphere and influences weather patterns. This test aims to assess how well "dry dynamical cores," the foundational components of weather and climate models that handle air movement and temperature changes, perform in idealized conditions. 

The simulation sets the model's atmosphere to an initial steady state, designed to be a simple, realistic representation of atmospheric conditions using the adiabatic (no heat exchange) and inviscid (no friction) primitive equations. The test first checks whether each model can maintain this steady, zonal (west-to-east) state without developing any unintended changes. After verifying this, the next step is to introduce a small disturbance, or perturbation, which triggers the growth of a baroclinic wave. The wave then evolves over several simulated days, allowing the researchers to observe how accurately different models handle the wave's development and movement.

This test provides a standard way to assess how different atmospheric models handle the development of baroclinic waves. The results help identify which models are more accurate and can serve as benchmarks for model improvement, ultimately contributing to better simulations of atmospheric behavior in weather and climate predictions.

In the UFS WM, the idealized baroclinic wave test case is an atmosphere-only, :term:`dycore`-only forecast run at C192 resolution with 127 vertical levels. It uses default values from the WM's ``export_fv3`` function, along with default values for a tiled grid namelist (from ``export_tiled``). These initial values are set based on values from `default_vars.sh <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/default_vars.sh>`_. 

The test is set to run a dynamics-only 24-hour forecast from 2019-12-03 at 0z. However, it is recommended that users modify the case to run it as a 5-10 day forecast by setting the forecast length (``FHMAX``) to 120-240 hours in the test file (see :numref:`Section %s <bw-config>` for instructions). Users will also need to update ``OUTPUT_FH`` accordingly. This ensures that users will be able to see the disturbance (wave propagation) develop, starting about 120-150 hours into the forecast. 

=============================
Obtaining Data for HSD Cases
=============================

.. include:: ./doc-snippets/hsd_data.rst

.. _run-bw:

==================================
Running the Baroclinic Wave Case
==================================

This section explains how to run the baroclinic wave case described above using the ``ufs_test.sh`` script.

Clone the Repository
--------------------

.. include:: ./doc-snippets/clone_hsd.rst

Machine Configuration
-----------------------

.. include:: ./doc-snippets/hsd_machine_config.rst

.. _bw-config:

Test Configuration
----------------------

It is recommended that users adjust certain values in the baroclinic wave case. Currently, the forecast length (``FHMAX``) is set to 24 hours to maintain the regression testing baselines with reduced resource consumption. However, it is recommended that users run the case for 5-10 days (120-240 hours) to see the disturbance (wave propagation) develop about 120-150 hours into the forecast. To modify the case, open ``${UFS_WM}/tests-dev/test_cases/tests/baroclinic_wave`` using ``vi``/``vim`` or a code editor. Then, add ``FHMAX`` and update ``OUTPUT_FH`` to extend by increments of 6 to the new ``FHMAX``. 

.. code-block:: console

   export FHMAX=120      # (or 240) 
   export OUTPUT_FH='0 6 12 18 24 30 36 42 48 54 60 66 72 78 84 90 96 102 108 114 120'

In general, it is preferable to make ``FHMAX`` a multiple of 24. 

.. note:: 

   On Jet, users will also need to adjust ``${UFS_WM}/tests/fv3_conf/compile_slurm.IN_jet`` in order to manage memory requirements for longer runs of the ``baroclinic_wave`` test. Users will need to change the number of tasks per node from 8 to 6 and add ``#SBATCH --mem=0``. 

   The file should include: 

   .. code-block:: console 
      
      #SBATCH --ntasks-per-node=6
      #SBATCH --mem=0

Baseline Configuration
----------------------

.. include:: ./doc-snippets/hsd_baseline_config.rst

Running Tests
-------------

.. include:: ./doc-snippets/hsd_run_tests.rst

Example:
^^^^^^^^^

Users with access to the ``epic`` account can run the ``baroclinic_wave`` test case with the ``intel`` compiler on :term:`RDHPCS` where they have access using the following command:

.. code-block:: console

   ./ufs_test.sh -a epic -s -c -k -r -n "baroclinic_wave intel"

Running Multiple Cases
^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ./doc-snippets/hsd_run_multiple.rst

Checking Results
-----------------

.. include:: ./doc-snippets/hsd_check_results.rst

For example, to monitor progress or check results for the ``baroclinic_wave`` case, run:

.. code-block:: console

   tail -f ${UFS_WM}/tests-dev/run_dir/baroclinic_wave_intel/err
   tail -f ${UFS_WM}/tests-dev/run_dir/baroclinic_wave_intel/out

.. include:: ./doc-snippets/hsd_notes.rst
