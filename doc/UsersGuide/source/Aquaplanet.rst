.. role:: raw-html(raw)
    :format: html

.. _aquaplanet:

*********************
Aquaplanet Test Case
*********************

The Aquaplanet test case is an idealized atmosphere-only forecast configuration designed to study atmospheric dynamics in a simplified Earth-like setting where all land is replaced with ocean. This configuration removes the complexity of land-surface interactions, topography, and regional variations, allowing researchers to focus on fundamental atmospheric processes such as tropical convection, Hadley circulation, jet stream dynamics, and the global energy budget.

The test case runs at C48 resolution with the ``FV3_GFS_v17_p8_ugwpv1`` physics suite, which is documented in the `Developmental Testbed Center (DTC) UFS Aquaplanet scientific documentation <https://ncar.github.io/ccpp-physics/HSD/_g_f_s_v17_p8_ugwpv1_aquaplanet_page.html>`_. Initial conditions are created by modifying standard GFS data to represent an aquaplanet configuration: Sea surface temperatures (SST) follow an idealized latitudinal profile, all land is converted to ocean, topography is set to zero, and sea ice is removed. The atmospheric initial state is configured with idealized vertical profiles of temperature, humidity, and winds appropriate for an aquaplanet simulation.

A key feature of this configuration is the 90-day spin-up period required to allow the model to adjust from its initial state to a balanced aquaplanet climate. After spin-up, the simulation can be run for extended periods to study seasonal variations, climate statistics, and atmospheric phenomena in the absence of land-surface influences.

============================
Obtaining Data for HSD Cases
============================

.. include:: ./doc-snippets/hsd_data.rst

Note that the 90-day spin-up for the aquaplanet case has already been completed for this data. Advanced users may prefer to try the 90-day spin-up themselves using the instructions in :numref:`Section %s <setup-aquaplanet>`. However, it is recommended that all users try the case first with the prestaged data. 

.. _run-aquaplanet:

==============================
Running the Aquaplanet Case
==============================

This section explains how to run the Aquaplanet case using the ``ufs_test.sh`` script with prestaged initial conditions. This is the recommended way to run the case for most users.

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

where:

- ``-s``: uses tests-dev; symlinks sharable test scripts
- ``-c``: prevents scripts from comparing with previous results (by creating new baselines)
- ``-k``: keeps run directory
- ``-r``: uses rocoto scheduler (``-e`` will use ecFlow)
- ``-n "aquaplanet intel"``: uses the aquaplanet test case with intel compiler

Running Extended Experiments
-----------------------------

After running the base test case, users can extend the simulation by running multiple restart segments. For example, to run a 1-year experiment broken into four 3-month segments:

Navigate to the run directory: 

.. code-block:: console

   cd ${UFS_WM}/tests-dev/run_dir/aquaplanet_intel

Users may optionally save this location in an environment variable(e.g., ``${RUNDIR}``):

.. code-block:: console

   export RUNDIR=${UFS_WM}/tests-dev/run_dir/aquaplanet_intel

For each 3-month segment, update ``fhrot`` and ``nhours_fcst`` in ``model_configure``:

**First segment (months 4-6):**

.. code-block:: console

   nhours_fcst:             4320
   fhrot:                   2160

**Second segment (months 7-9):**

.. code-block:: console

   nhours_fcst:             6480
   fhrot:                   4320

**Third segment (months 10-12):**

.. code-block:: console

   nhours_fcst:             8640
   fhrot:                   6480

Then, update the walltime in the ``job_card`` to run for approximately 8 hours (480 minutes):  

.. code-block:: bash
   
   #SBATCH --time=480

Some systems may succeed in less time; others may need more time, but this is a reasonable starting point. 

Next, submit the job:

.. code-block:: console

   sbatch job_card

After each segment completes, rename and move restart files to the ``INPUT`` directory:

.. code-block:: console

   cd RESTART/
   for file in 20*.060000.*.nc; do mv "$file" "${file#20*.060000.}"; done
   mv 20*.060000.coupler.res coupler.res
   mv * ../INPUT/.
   cd ..

Then, repeat the steps above to adjust ``nhours_fcst`` and ``fhrot`` and submit the job card. 

Checking Results
----------------

To monitor progress or check results for the ``aquaplanet`` case, run:

.. code-block:: console

   tail -f ${UFS_WM}/tests-dev/run_dir/aquaplanet_intel/out
   tail -f ${UFS_WM}/tests-dev/run_dir/aquaplanet_intel/err

The ``out`` file contains messages output during the model run. The ``err`` file reports on errors and is useful if the case does not finish properly. Additionally, the model will produce ``atmf*.nc`` and ``sfcf*.nc`` files every 24 hours. 

When the test case finishes running, the ``out`` file should include a resource statistics summary at the bottom:

.. code-block:: console

   0: *****************RESOURCE STATISTICS*******************************
   0: The total amount of wall time                        = 22639.581802
   0: The total amount of time in user mode                = 22222.145691
   0: The total amount of time in sys mode                 = 189.322066
   0: The maximum resident set size (KB)                   = 1769600
   0: Number of page faults without I/O activity           = 89229617
   0: Number of page faults with I/O activity              = 99093
   0: Number of times filesystem performed INPUT           = 1368912
   0: Number of times filesystem performed OUTPUT          = 115520
   0: Number of Voluntary Context Switches                 = 347098
   0: Number of InVoluntary Context Switches               = 758822
   0: *****************END OF RESOURCE STATISTICS*************************
   0:
   Model ended: Thu Feb 19 10:46:07 UTC 2026

.. _plotting-aquaplanet:

=======================
Plotting Script
=======================

A plotting script is available to generate seasonal mean plots for key atmospheric variables. The script is located at:

.. code-block:: console

   ${UFS_WM}/tests-dev/test_cases/utils/plot_aq.sh

.. note::

   This plotting script is currently configured to run on Ursa only.

By default, this script creates seasonal means for three variables:

- Jet stream characteristics
- Precipitation patterns
- Temperature at 500 mb

The script uses staged data from a 1-year control simulation. If you want to use this dataset on other machines (which have access to :term:`HPSS`), you can retrieve it from HPSS at ``/5year/NCEPDEV/emc-meso/Ratko.Vasic/AQUAPLANET/1yr-results.tar``. Alternatively, users can adjust the ``data_path`` variable to point to their own model output from an extended run. 

Customizing the Plotting Script
--------------------------------

To use the plotting script with user-generated data:

1. Copy the script to your run directory.

   .. code-block:: bash

      cd ${UFS_WM}/tests-dev/run_dir/aquaplanet_intel
      cp ${UFS_WM}/tests-dev/test_cases/utils/plot_aq.sh .

2. Edit the variable ``data_path`` to point to either the prestaged model output data on Ursa or the data generated from your own model run. For example:

   .. code-block:: bash

      data_path=./

   This will use data from the run directory (rather than data from the prestaged model output) to produce the plots. 

3. Adjust the seasonal timing variables (``winter_start``, ``spring_start``, etc.), which are given in hours from the start of the 90-day spin-up run (hour 0).

   .. code-block:: bash

      winter_start=2160
      spring_start=4320
      summer_start=6480
      fall_start=8640

4. Modify ``season_len`` to set the length of each season in days (typically 90 days, but can be set to 365 for annual means).

   .. code-block:: bash

      season_len=90

.. _setup-aquaplanet:

============================================================
Advanced: Setting Up the Aquaplanet Experiment from Scratch
============================================================

.. note::

   This section is **optional**. Most users can run the Aquaplanet case using the ``ufs_test.sh`` method described above, which uses pre-staged initial conditions. The steps below are for advanced users who wish to create their own initial conditions from scratch. :term:`HPSS` access is required for the UFS_UTILS steps. 

The from-scratch setup involves editing orography, SST, ice, and sea-land mask files to create an aquaplanet configuration, then running a 90-day spin-up to allow the model to reach a balanced state. This process produces the same initial conditions that are provided pre-staged for the standard test case.

At a high level, the steps are:

1. Build the UFS Weather Model and run a baseline ``control_c48`` test to generate a run directory.
2. Build UFS_UTILS, and use ``gdas_init`` to generate raw atmospheric initial conditions.
3. Use the `Aquaplanet tools <https://github.com/NOAA-EPIC/Aquaplanet>`_ to modify SST, ice, sea-land mask, orography, and atmospheric/surface profiles to represent an aquaplanet.
4. Regenerate initial conditions with the modified files.
5. Apply minor source code changes and recompile the model.
6. Run a 90-day spin-up simulation (in three 30-day restart segments).

Create Working Directory (optional)
-------------------------------------

Users can create a new directory for their aquaplanet experiment or use an existing directory. To create a new directory, run: 

.. code-block:: console

   mkdir /path/to/ap-expt
   cd /path/to/ap-expt

where ``/path/to/ap-expt`` is the path to the directory where the user plans to perform the following steps (e.g., ``/Users/Joe.Schmoe/ap-expt``). 

Optionally, users can save this directory path in an environment variable (e.g., ``${AP_EXPT}``) to avoid typing out full path names later. 

.. code-block:: console

   export AP_EXPT=$PWD

In this documentation, ``${AP_EXPT}`` is used, but users are welcome to choose another name for this variable if they prefer. 

Clone Required Repositories
----------------------------

Clone the UFS Weather Model:

.. code-block:: console

   git clone --recursive https://github.com/ufs-community/ufs-weather-model.git

Clone UFS_UTILS:

.. code-block:: console

   git clone --recursive https://github.com/ufs-community/UFS_UTILS.git

Clone the Aquaplanet tools:

.. code-block:: console

   git clone https://github.com/NOAA-EPIC/Aquaplanet

Build UFS Weather Model
------------------------

Navigate to the UFS Weather Model test directory:

.. code-block:: console

   cd ${AP_EXPT}/ufs-weather-model/tests/

Edit ``rt.conf`` to include only the following two lines:

.. code-block:: console

   COMPILE | atm_dyn32 | intel | -DAPP=ATM -DCCPP_SUITES=FV3_GFS_v17_p8_ugwpv1 -D32BIT=ON | | fv3 |
   RUN | control_c48 | | baseline |

Execute the regression test to compile the model and create a baseline run directory:

.. code-block:: console

   ./rt.sh -a <account> -k -e

Replace ``<account>`` with your project code (e.g., ``epic``). Save the location of the ``control_c48_intel`` run directory for later use. The ``-e`` flag (or the ``-r``) flag are not required but can be useful. 

Build UFS_UTILS
---------------

Compile the UFS_UTILS tools:

.. code-block:: console

   cd ${AP_EXPT}/UFS_UTILS/
   module purge
   module use $PWD/modulefiles
   module load build.<platform>.intelllvm
   ./build_all.sh

Replace ``<platform>`` with your system name (e.g., ``ursa``).

Link the fix directories:

.. code-block:: console

   cd fix/
   ./link_fixdirs.sh emc <platform>

Since orography files need to be edited, remove the link and copy the files. The fix file locations for different platforms are listed in ``link_fixdirs.sh``. For example, on Ursa, users would run:

.. code-block:: console

   rm orog
   mkdir orog
   cp /scratch3/NCEPDEV/global/role.glopara/fix/orog/20240917/*.nc orog/.
   cp /scratch3/NCEPDEV/global/role.glopara/fix/orog/20240917/*.dat orog/.
   cp -r /scratch3/NCEPDEV/global/role.glopara/fix/orog/20240917/C48/ orog/.

.. _gen-initial-atmos-data:

Generate Initial Atmospheric Data With UFS_UTILS
--------------------------------------------------

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

where ``<path_to_UFS_UTILS>`` is replaced with the actual path to the UFS_UTILS repository. 

Edit the driver script (e.g., ``driver.ursa.sh``) to set:

.. code-block:: console

   PROJECT_CODE=<your_account>

Comment out the line that removes the extract directory:

.. code-block:: console

   # rm -fr $EXTRACT_DIR

Run the driver script:

.. code-block:: console

   ./driver.<platform>.sh

Check that results are properly generated in ``${EXTRACT_DIR}`` and ``${OUTDIR}``. This can take some time, so if there is only an empty ``${EXTRACT_DIR}``, run ``squeue -u $USER`` (on systems with Slurm) or ``qstat -u $USER`` (on systems with PBS Pro) to ensure that the job is still running. Eventually, both directories will contain a file named ``gfs.20251015``. 

Create Idealized SST Profile Using Aquaplanet Tools
-----------------------------------------------------

Navigate to the SST profile tool directory:

.. code-block:: console

   cd ${AP_EXPT}/Aquaplanet/sst-profile/

Copy the SST climatology file from your UFS Weather Model run directory, which will be referred to as ``${RUNDIR}`` in this documentation:

.. code-block:: console

   cp ${RUNDIR}/RTGSST.1982.2012.monthly.clim.grb .

Compile and run the tool:

.. code-block:: console

   ./compile.sh
   ./sst-profile.x

Copy the modified file back to your run directory:

.. code-block:: console

   cp new-RTGSST.1982.2012.monthly.clim.grb ${RUNDIR}/RTGSST.1982.2012.monthly.clim.grb

Edit Global Ice Data
--------------------

Navigate to the glacier tool directory:

.. code-block:: console

   cd ${AP_EXPT}/Aquaplanet/glacier/

Copy the glacier file from your run directory:

.. code-block:: console

   cp ${RUNDIR}/global_glacier.2x2.grb .

Compile and run the tool:

.. code-block:: console

   ./compile.sh
   ./glacier.x

Copy the modified file back:

.. code-block:: console

   cp new-global_glacier.2x2.grb ${RUNDIR}/global_glacier.2x2.grb

Edit Monthly Ice Data
---------------------

Navigate to the ice monthly tool directory:

.. code-block:: console

   cd ${AP_EXPT}/Aquaplanet/ice-monthly/

Copy the ice climatology file from your run directory:

.. code-block:: console

   cp ${RUNDIR}/IMS-NIC.blended.ice.monthly.clim.grb .

Compile and run the tool:

.. code-block:: console

   ./compile.sh
   ./ice-monthly.x

.. note::

   If the tool fails due to lack of memory, use the provided ``job_card`` file as a template and submit to a batch node.

Copy the modified file back:

.. code-block:: console

   cp new-IMS-NIC.blended.ice.monthly.clim.grb ${RUNDIR}/IMS-NIC.blended.ice.monthly.clim.grb

Edit Sea-Land Mask Data
-----------------------

Navigate to the sea-land mask tool directory:

.. code-block:: console

   cd ${AP_EXPT}/Aquaplanet/slmask/

Copy the sea-land mask file from your run directory:

.. code-block:: console

   cp ${RUNDIR}/global_slmask.t62.192.94.grb .

Compile and run the tool:

.. code-block:: console

   ./compile.sh
   ./slmask.x

Copy the modified file back:

.. code-block:: console

   cp new-global_slmask.t62.192.94.grb ${RUNDIR}/global_slmask.t62.192.94.grb

Setup Atmospheric Profiles
---------------------------

Navigate to the atmospheric profile tool directory:

.. code-block:: console

   cd ${AP_EXPT}/Aquaplanet/atmos-profile

Copy the atmospheric analysis file from ``${EXTRACT_DIR}`` (set in :numref:`Step %s <gen-initial-atmos-data>` to ``${AP_EXPT}/UFS_UTILS/input``):

.. code-block:: console

   cp ${EXTRACT_DIR}/gfs.20251015/06/atmos/gfs.t06z.atmanl.nc .

Compile and load required modules:

.. code-block:: console

   ./compile.sh
   module purge
   module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core
   module load stack-oneapi/2024.2.1
   module load stack-intel-oneapi-mpi/2021.13
   module load netcdf-fortran/4.6.1

Run the tool:

.. code-block:: console

   ./profile.x

Copy the modified atmospheric file back to the extract directory:

.. code-block:: console

   cp gfs.t06z.atmanl.nc ${EXTRACT_DIR}/gfs.20251015/06/atmos/gfs.t06z.atmanl.nc

Setup Surface Profiles
-----------------------

Navigate to the surface profile tool directory:

.. code-block:: console

   cd ${AP_EXPT}/Aquaplanet/sfc-profile/

Copy the surface analysis file from ``${EXTRACT_DIR}`` (set in :numref:`Step %s <gen-initial-atmos-data>` to ``${AP_EXPT}/UFS_UTILS/input``):

.. code-block:: console

   cp ${EXTRACT_DIR}/gfs.20251015/06/atmos/gfs.t06z.sfcanl.nc .

Compile and load required modules:

.. code-block:: console

   ./compile.sh
   module purge
   module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core
   module load stack-oneapi/2024.2.1
   module load stack-intel-oneapi-mpi/2021.13
   module load netcdf-fortran/4.6.1

Run the tool:

.. code-block:: console

   ./profile.x

Copy the modified surface file back to the extract directory:

.. code-block:: console

   cp gfs.t06z.sfcanl.nc ${EXTRACT_DIR}/gfs.20251015/06/atmos/gfs.t06z.sfcanl.nc

Prepare Orography Fields
-------------------------

Navigate to the orography profile tool directory:

.. code-block:: console

   cd ${AP_EXPT}/Aquaplanet/orog-profile/

Compile the tool:

.. code-block:: console

   ./compile.sh

Verify that ``list.txt`` contains the correct relative path to fix orography files, then run:

.. code-block:: console

   ./patch.sh

This will modify the orography files in place to set land mask to ocean and height to zero.

Regenerate Initial Conditions
------------------------------

Now that the correct fixed files are in place, return to UFS_UTILS and run ``chgres`` again.

Navigate back to the ``gdas_init`` utility directory:

.. code-block:: console

   cd ${AP_EXPT}/UFS_UTILS/util/gdas_init/

.. important::

   In the ``config`` file, set ``EXTRACT_DATA`` to ``no`` since data is already extracted.

Delete or rename the output directory. For example:

.. code-block:: console

   rm -rf $OUTDIR

or:

.. code-block:: console

   mv $OUTDIR $OUTDIR.backup

Run the driver script:

.. code-block:: console

   ./driver.<platform>.sh

Once the job has finished running, check the results in:

.. code-block:: console

   ../../output/gfs.20251015/06/model/atmos/input/

Add Noah-MP Variables to Surface Files
---------------------------------------

Since the GFS_v17 physics suite is used, empty fields (``sheleg``, ``snwdph``, and ``zorl``) must be added to the surface files:

.. code-block:: console

   cd ${AP_EXPT}/Aquaplanet/noah-MP-vars/

Copy the surface files from the output directory (set in :numref:`Step %s <gen-initial-atmos-data>` to ``${AP_EXPT}/UFS_UTILS/output``):

.. code-block:: console

   cp $OUTDIR/gfs.20251015/06/model/atmos/input/sfc_* .

Compile and load required modules:

.. code-block:: console

   ./compile.sh
   module purge
   module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core
   module load stack-oneapi/2024.2.1
   module load stack-intel-oneapi-mpi/2021.13
   module load netcdf-fortran/4.6.1

Run the tool:

.. code-block:: console

   ./add-vars.x

Copy the modified files back:

.. code-block:: console

   cp sfc_* $OUTDIR/gfs.20251015/06/model/atmos/input/.

Copy Initial Conditions to Run Directory
-----------------------------------------

Copy all generated initial condition files from ``${OUTDIR}`` (set in to ``${AP_EXPT}/UFS_UTILS/output`` in :numref:`Step %s <gen-initial-atmos-data>`) to your run directory:

.. code-block:: console

   cp $OUTDIR/gfs.20251015/06/model/atmos/input/* ${RUNDIR}/INPUT/.

Copy the modified orography files:

.. code-block:: console

   cd ${RUNDIR}/INPUT/
   cp ${AP_EXPT}/UFS_UTILS/fix/orog/C48/C48.mx500_oro_data.tile1.nc oro_data.tile1.nc
   cp ${AP_EXPT}/UFS_UTILS/fix/orog/C48/C48.mx500_oro_data.tile2.nc oro_data.tile2.nc
   cp ${AP_EXPT}/UFS_UTILS/fix/orog/C48/C48.mx500_oro_data.tile3.nc oro_data.tile3.nc
   cp ${AP_EXPT}/UFS_UTILS/fix/orog/C48/C48.mx500_oro_data.tile4.nc oro_data.tile4.nc
   cp ${AP_EXPT}/UFS_UTILS/fix/orog/C48/C48.mx500_oro_data.tile5.nc oro_data.tile5.nc
   cp ${AP_EXPT}/UFS_UTILS/fix/orog/C48/C48.mx500_oro_data.tile6.nc oro_data.tile6.nc

Modify Model Source Code
-------------------------

Two minor modifications are needed to the UFS Weather Model source code for the aquaplanet configuration.

**Fix Solar Zenith Angle**

Edit the file ``./ufs-weather-model/UFSATM/ccpp/physics/physics/Radiation/radiation_astronomy.f``. After line 600, add:

.. code-block:: console

   ! AQUAPLANET TEST CASE
         solcon=1370.
         dlt=0.
   ! AQUAPLANET TEST CASE

**Reduce Land-Surface Warnings**

Edit the file ``./ufs-weather-model/UFSATM/fv3/atmos_cubed_sphere/tools/fv_diagnostics.F90``. After lines 4292 and 4349, where it says ``if ( present(bad_range) ) then``, add:

.. code-block:: console

   ! AQUAPLANET TEST CASE
         return
   ! AQUAPLANET TEST CASE

Recompile the model after making these changes. To do this, run the following commands:

.. code-block:: console

   module use ${AP_EXPT}/ufs-weather-model/modulefiles
   module load ufs_<platform>_intel
   cd ${AP_EXPT}/ufs-weather-model
   mkdir build
   cd build
   cmake .. -DAPP=ATM -D32BIT=ON -DCCPP_SUITES=FV3_GFS_v17_p8_ugwpv1
   make -j4

Note that users can leave off the ``-j`` flag after the ``make`` command or change the integer number after it to run with fewer or more processors. 

Running these commands will produce an executable called ``ufs_model``. Rename this executable to ``fv3.exe``. Then, update the executable in your run directory by copying the newly compiled executable in your build directory to your run directory. 

.. code-block:: console

   cp ufs_model fv3.exe
   mv fv3.exe ${RUNDIR}

Running the 90-Day Spin-Up
---------------------------

The aquaplanet simulation requires a 90-day spin-up period to reach equilibrium. This is run in three 30-day (720-hour) segments with restarts.

**Initial 30-Day Run**

Edit the ``job_card`` file in your run directory:

.. code-block:: console

   #SBATCH --time=480

Edit the ``model_configure`` file:

.. code-block:: console

   start_year:              2025
   start_month:             10
   start_day:               15
   nhours_fcst:             720
   quilting_restart:        .false.
   output_fh:               24 -1

Submit the job:

.. code-block:: console

   sbatch job_card

In general, each 30-day run will take 5.5-6 hours, but this can vary. User can view progress with the ``tail`` command:

.. code-block:: console

   tail -f out

**First Restart (Days 31-60)**

After the first month completes, prepare for the restart:

.. code-block:: console

   cd RESTART/
   for file in 20*.060000.*.nc; do mv "$file" "${file#20*.060000.}"; done
   mv 20*.060000.coupler.res coupler.res
   mv * ../INPUT/.
   cd ..

Edit ``input.nml`` and set:

.. code-block:: console

   make_nh = .false.
   na_init = 0
   external_ic = .false.
   nggps_ic = .false.
   mountain = .true.
   warm_start = .true.

Edit ``model_configure`` and set:

.. code-block:: console

   nhours_fcst:             1440
   fhrot:                   720

Submit the job:

.. code-block:: console

   sbatch job_card

**Second Restart (Days 61-90)**

After the second month completes:

.. code-block:: console

   cd RESTART/
   for file in 20*.060000.*.nc; do mv "$file" "${file#20*.060000.}"; done
   mv 20*.060000.coupler.res coupler.res
   mv * ../INPUT/.
   cd ..

Keep ``input.nml`` settings from the first restart. Edit ``model_configure`` and set:

.. code-block:: console

   nhours_fcst:             2160
   fhrot:                   1440

Submit the job:

.. code-block:: console

   sbatch job_card

After this final segment completes, the files in the ``RESTART`` directory represent spun-up initial conditions. Rename and copy these files to the ``INPUT`` directory --- these are the new initial conditions for your experiment (spun up for 90 days). Save them for future use.

To continue running experiments with these self-generated initial conditions, follow the instructions in :numref:`Section %s <run-aquaplanet>`, using the ``input.nml`` from the restart runs and updating ``fhrot`` and ``nhours_fcst`` in ``model_configure`` as described in the extended experiments section.
