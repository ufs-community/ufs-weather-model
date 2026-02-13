.. role:: raw-html(raw)
    :format: html

.. _aquaplanet:

*********************
AquaPlanet Test Case
*********************

The AquaPlanet test case is an idealized atmosphere-only forecast configuration designed to study atmospheric dynamics in a simplified Earth-like setting where all land is replaced with ocean. This configuration removes the complexity of land-surface interactions, topography, and regional variations, allowing researchers to focus on fundamental atmospheric processes such as tropical convection, Hadley circulation, jet stream dynamics, and the global energy budget.

The test case runs at C48 resolution with the ``FV3_GFS_v17_p8_ugwpv1`` physics suite. Initial conditions are created by modifying standard GFS data to represent an aquaplanet configuration: sea surface temperatures (SST) follow an idealized latitudinal profile, all land is converted to ocean, topography is set to zero, and sea ice is removed. The atmospheric initial state is configured with idealized vertical profiles of temperature, humidity, and winds appropriate for an aquaplanet simulation.

A key feature of this configuration is the 90-day spin-up period required to allow the model to adjust from its initial state to a balanced aquaplanet climate. After spin-up, the simulation can be run for extended periods to study seasonal variations, climate statistics, and atmospheric phenomena in the absence of land-surface influences.

============================
Obtaining Data for HSD Cases
============================

.. include:: ./doc-snippets/hsd_data.rst

.. _run-aquaplanet:

==============================
Running the AquaPlanet Case
==============================

This section explains how to run the AquaPlanet case using the ``ufs_test.sh`` script with pre-staged initial conditions. This is the recommended way to run the case for most users.

Clone the Repository
--------------------

.. include:: ./doc-snippets/clone_hsd.rst

Machine Configuration
---------------------

.. include:: ./doc-snippets/hsd_machine_config.rst

Running Tests
-------------

.. include:: ./doc-snippets/hsd_run_tests.rst

Example
^^^^^^^

Users with access to the ``epic`` account can run the ``aquaplanet`` test case with the ``intel`` compiler on :term:`RDHPCS` where they have access using the following command:

.. code-block:: console

   ./ufs_test.sh -a epic -s -c -k -r -n "aquaplanet intel"

Where:

- ``-s``: use tests-dev, symlink sharable test scripts
- ``-c``: prevents scripts from comparing with previous results
- ``-k``: keep run directory
- ``-r``: use rocoto scheduler (``-e`` will use ecFlow)
- ``-n "aquaplanet intel"``: use the aquaplanet test case with intel compiler

Running Extended Experiments
-----------------------------

After running the base test case, users can extend the simulation by running multiple restart segments. For example, to run a 1-year experiment broken into four 3-month segments:

Navigate to the run directory (linked as ``./tests-dev/run_dir``).

Use the ``input.nml`` file configured for restart runs (with ``warm_start = .true.``).

For each 3-month segment, update ``fhrot`` and ``nhours_fcst`` in ``model_configure``:

**First segment (months 4-6):**

.. code-block:: console

   fhrot:                   2160
   nhours_fcst:             4320

**Second segment (months 7-9):**

.. code-block:: console

   fhrot:                   4320
   nhours_fcst:             6480

**Third segment (months 10-12):**

.. code-block:: console

   fhrot:                   6480
   nhours_fcst:             8640

After each segment completes, rename and move restart files to the ``INPUT`` directory:

.. code-block:: console

   cd RESTART/
   for file in 20*.060000.*.nc; do mv "$file" "${file#20*.060000.}"; done
   mv 20*.060000.coupler.res coupler.res
   mv * ../INPUT/.

Checking Results
----------------

.. include:: ./doc-snippets/hsd_check_results.rst

For example, to monitor progress or check results for the ``aquaplanet`` case, run:

.. code-block:: console

   tail -f ${UFS_WM}/tests-dev/run_dir/aquaplanet_intel/err
   tail -f ${UFS_WM}/tests-dev/run_dir/aquaplanet_intel/out

   .. _plotting-aquaplanet:

=======================
Plotting Script
=======================

A plotting script is available to generate seasonal mean plots for key atmospheric variables. The script is located at:

.. code-block:: console

   ./tests-dev/test_cases/utils/plot_aq.sh

.. note::

   This plotting script is currently configured to run on Hera and Ursa only.

By default, this script creates seasonal means for three variables:

- Jet stream characteristics
- Precipitation patterns
- Temperature at 500 mb

The script uses staged data from a 1-year control simulation. If you want to use this dataset on other machines (which have access to :term:`HPSS`), you can retrieve it from HPSS:

.. code-block:: console

   /5year/NCEPDEV/emc-meso/Ratko.Vasic/AQUAPLANET/1yr-results.tar

Customizing the Plotting Script
--------------------------------

To use the plotting script with user-generated data:

1. Copy the script to your run directory.
2. Edit the variable ``out_pth`` to point to your output location:

   .. code-block:: bash

      out_pth=./

3. Adjust the seasonal timing variables (``winter_start``, ``spring_start``, etc.), which are given in hours from the start of the 90-day spin-up run (hour 0).
4. Modify ``season_len`` to set the length of each season in days (typically 90 days, but can be set to 365 for annual means).