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

.. _setup-aquaplanet:

==========================================================
Advanced: Setting Up the AquaPlanet Experiment from Scratch
==========================================================

.. note::

   This section is **optional**. Most users can run the AquaPlanet case using the ``ufs_test.sh`` method described above, which uses pre-staged initial conditions. The steps below are for advanced users who wish to create their own initial conditions from scratch.

The from-scratch setup involves editing orography, SST, ice, and sea-land mask files to create an aquaplanet configuration, then running a 90-day spin-up to allow the model to reach a balanced state. This process produces the same initial conditions that are provided pre-staged for the standard test case.

At a high level, the steps are:

1. Build the UFS Weather Model and run a baseline ``control_c48`` test to generate a run directory.
2. Build UFS_UTILS and use ``gdas_init`` to generate raw atmospheric initial conditions.
3. Use the `Aquaplanet tools <https://github.com/RatkoVasic-NOAA/Aquaplanet>`_ to modify SST, ice, sea-land mask, orography, and atmospheric/surface profiles to represent an aquaplanet.
4. Regenerate initial conditions with the modified files.
5. Apply minor source code changes and recompile the model.
6. Run a 90-day spin-up simulation (in three 30-day restart segments).

Clone Required Repositories
----------------------------

Clone the UFS Weather Model:

.. code-block:: console

   git clone https://github.com/ufs-community/ufs-weather-model.git
   cd ufs-weather-model/
   git submodule update --init --recursive

Clone UFS_UTILS:

.. code-block:: console

   git clone https://github.com/ufs-community/UFS_UTILS.git
   cd UFS_UTILS/
   git submodule update --init --recursive

Clone the aquaplanet tools:

.. code-block:: console

   git clone https://github.com/RatkoVasic-NOAA/Aquaplanet

Build UFS Weather Model
------------------------

Navigate to the UFS Weather Model test directory:

.. code-block:: console

   cd ufs-weather-model/tests/

Edit ``rt.conf`` to include only the following two lines:

.. code-block:: console

   COMPILE | atm_dyn32 | intel | -DAPP=ATM -DCCPP_SUITES=FV3_GFS_v17_p8_ugwpv1 -D32BIT=ON | | fv3 |
   RUN | control_c48 | | baseline |

Execute the regression test to compile the model and create a baseline run directory:

.. code-block:: console

   ./rt.sh -a <account> -k

Replace ``<account>`` with your project code (e.g., ``epic``). Save the location of the ``control_c48_intel`` run directory for later use.

Build UFS_UTILS
---------------

Compile the UFS_UTILS tools:

.. code-block:: console

   cd UFS_UTILS/
   module purge
   module use $PWD/modulefiles
   module load build.<platform>.intelllvm
   ./build_all.sh

Replace ``<platform>`` with your system name (e.g., ``ursa``).

Link the fix directories:

.. code-block:: console

   cd fix/
   ./link_fixdirs.sh emc <platform>

Since orography files need to be edited, remove the link and copy the files:

.. code-block:: console

   rm orog
   mkdir orog
   cp /scratch3/NCEPDEV/global/role.glopara/fix/orog/20240917/*.nc orog/.
   cp /scratch3/NCEPDEV/global/role.glopara/fix/orog/20240917/*.dat orog/.
   cp -r /scratch3/NCEPDEV/global/role.glopara/fix/orog/20240917/C48/ orog/.

Generate Initial Atmospheric Data
----------------------------------

Navigate to the ``gdas_init`` utility directory:

.. code-block:: console

   cd ../util/gdas_init/

Edit the ``config`` file to set the following variables:

.. code-block:: console

   EXTRACT_DIR=<path_to_UFS_UTILS>/input/
   EXTRACT_DATA=yes
   yy=2025
   mm=10
   dd=15
   CRES_HIRES=C48
   OUTDIR=<path_to_UFS_UTILS>/output/

Edit the driver script (e.g., ``driver.ursa.sh``) to set:

.. code-block:: console

   PROJECT_CODE=<your_account>

Comment out the line that removes the extract directory:

.. code-block:: console

   # rm -fr $EXTRACT_DIR

Run the driver script:

.. code-block:: console

   ./driver.<platform>.sh

Check that results are properly generated in ``$EXTRACT_DIR`` and ``$OUTDIR``.

Create Idealized SST Profile
-----------------------------

Navigate to the SST profile tool directory:

.. code-block:: console

   cd ~/Aquaplanet/sst-profile/

Copy the SST climatology file from your UFS Weather Model run directory:

.. code-block:: console

   cp <run_directory>/RTGSST.1982.2012.monthly.clim.grb .

Compile and run the tool:

.. code-block:: console

   ./compile.sh
   ./sst-profile.x

Copy the modified file back to your run directory:

.. code-block:: console

   cp new-RTGSST.1982.2012.monthly.clim.grb <run_directory>/RTGSST.1982.2012.monthly.clim.grb

Edit Global Ice Data
--------------------

Navigate to the glacier tool directory:

.. code-block:: console

   cd ~/Aquaplanet/glacier/

Copy the glacier file from your run directory:

.. code-block:: console

   cp <run_directory>/global_glacier.2x2.grb .

Compile and run the tool:

.. code-block:: console

   ./compile.sh
   ./glacier.x

Copy the modified file back:

.. code-block:: console

   cp new-global_glacier.2x2.grb <run_directory>/global_glacier.2x2.grb

Edit Monthly Ice Data
---------------------

Navigate to the ice monthly tool directory:

.. code-block:: console

   cd ~/Aquaplanet/ice-monthly/

Copy the ice climatology file from your run directory:

.. code-block:: console

   cp <run_directory>/IMS-NIC.blended.ice.monthly.clim.grb .

Compile and run the tool:

.. code-block:: console

   ./compile.sh
   ./ice-monthly.x

.. note::

   If the tool fails due to lack of memory, use the provided ``job_card`` file as a template and submit to a batch node.

Copy the modified file back:

.. code-block:: console

   cp new-IMS-NIC.blended.ice.monthly.clim.grb <run_directory>/IMS-NIC.blended.ice.monthly.clim.grb