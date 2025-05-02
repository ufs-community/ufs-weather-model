.. _FAQ:

***
FAQ
***

==============================================================
How do I build and run a single test of the UFS Weather Model?
==============================================================

An efficient way to build and run the UFS Weather Model is to use the regression test (RT) script
(``rt.sh``). This script is widely used by model developers on :wm-wiki:`Tier 1 <Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>` and 2 platforms
and is described in :numref:`Section %s <run-wm>`. 

.. note::
   
   Users on Level 2-4 systems may need to perform additional steps prior to following the steps below. For example, they may need to :ref:`download data <GetData>` and :ref:`update files <other-systems>` with platform-specific information. 

For all systems, users will need to:

   #. Clone the source code and submodules as described in :numref:`Section %s <DownloadingWMCode>`; then
      navigate to the ``tests`` directory:

      .. code-block:: console

         git clone --recursive https://github.com/ufs-community/ufs-weather-model.git
         cd ufs-weather-model/tests

   #. Modify the ``rt.sh`` script to put the output in a run directory where they have write permissions. For example, on Hercules, users would update ``dprefix``:

      .. code-block:: console

         case ${MACHINE_ID} in
         ...
            hercules)
         ...
               dprefix="/work2/noaa/stmp/${USER}"
               DISKNM="/work/noaa/epic/hercules/UFS-WM_RT"
               STMP="${dprefix}/stmp"
               PTMP="${dprefix}/stmp"

   #. Run the ``rt.sh`` script: 
      
      * To run one specific test, such as ``control_c48``, use the ``-n`` flag to designate the name of the test and the type of compiler: 

         .. code-block:: console

            ./rt.sh -a <account_name> -k -n "control_c48 intel"

         where ``<account_name>`` is replaced with the name of an account where the user can charge computational resources. 
         The ``-k`` option will preserve the run directory after the forecast finishes. The :wm-repo:`rt.conf <blob/develop/tests/rt.conf>` file contains all of the currently maintained RTs. 

      * Users can run the entire RT suite using the ecFlow workflow manager:

         .. code-block:: console

            ./rt.sh -a <account_name> -e -k -l rt.conf
      
      * To run ``rt.sh`` using a custom configuration file and the Rocoto workflow manager, create a configuration file (e.g., ``my_tests.conf``) based on ``rt.conf``. For example, to run only a few S2S tests, create a file called ``s2s.conf``. 
         
         .. code-block:: console

            COMPILE | s2s | intel | -DAPP=S2S -DCCPP_SUITES=FV3_GFS_v17_coupled_p8,FV3_GFS_v17_coupled_p8_ugwpv1 |   | fv3 |
            RUN | cpld_control_c48         |                            | baseline |
            RUN | cpld_warmstart_c48       | - noaacloud                | baseline |
            RUN | cpld_restart_c48         | - noaacloud                |          | cpld_warmstart_c48

         Then run:

         .. code-block:: console

            ./rt.sh -a <account_name> -r -k -l s2s.conf

         adding additional arguments as desired. 

   #. Check ``${STMP}/FV3_RT/rt_PID/<test_name>`` for the model run, where ``PID`` is a process ID, and ``<test_name>`` refers to a specific test, such as ``control_c48_intel``.
      A successful test will produce an ``out`` file with an exit code at the bottom. ``exit code 0:0`` indicates a successful run. For example: 

      .. code-block:: console 

         Job 7430255 finished for user Joe.Schmoe in partition hera with exit code 0:0

      There is also a RESOURCE STATISTICS summary at the end of the test's ``out`` file. Errors will appear in the ``err`` file. Users can find log files with more detailed information in ``ufs-weather-model/tests/logs/log_<platform>`` (where platform is the name of the machine the user is running on, e.g., ``log_hercules``).
   
   #. When the build and run are complete, users can modify the namelist or ``model_configure`` files in the run directory (``${STMP}``) 
      and re-run their forecast/test with modifications by submitting the ``job_card`` file:

      .. code-block:: console

         qsub job_card
         # OR
         sbatch job_card

============================================
How do I change the length of the model run?
============================================
For individual RT tests, users can add the ``FHMAX`` variable to the test configuration file. For example, in the :wm-repo:`control_c48 <blob/develop/tests/tests/control_c48>` case, 
users can increase the forecast duration from the default (``DAYS*24`` --- or 24 hours, in this case) to 48 hours by adding the statement:

.. code-block:: console

   export FHMAX=48

Alternatively, users can set a different default value in :wm-repo:`default_vars.sh <blob/develop/tests/default_vars.sh>` by changing the ``FHMAX`` variable directly in ``default_vars.sh``
or they can modify the ``nhours_fcst`` variable in the ``model_configure*`` file that their experiment uses. 

To rerun a previously run test with a different forecast length, go to the run directory (usually ``${STMP}/FV3_RT/rt_PID/<test_name>``), and open the file named ``model_configure``.  
Change the variable ``nhours_fcst`` to the desired number of hours for the forecast and rerun by sumbitting the job card, as described above.

==============================================================
How do I set the output history interval?
==============================================================

The interval at which output (history) files are written is controlled via the ``model_configure*`` files. 
When using the RT framework, users can adjust values in the test file for the test they plan to run, and these will be fed into the appropriate ``model_configure`` file.
To adjust the default values for entire sets of tests, values can be modified in the ``tests/default_vars.sh`` script. 
:numref:`Table %s <OutputControl>` describes the relevant variables.  

.. _OutputControl:

.. list-table:: *Variables used to control the output file frequency*
   :widths: 15 10 10 30
   :header-rows: 1

   * - Namelist variable
     - Location
     - Default Value in ``export_fv3``
     - Description
   * - OUTPUT_FH
     - ``model_configure``
     - "12 -1"
     - Array listing the forecast output frequency; this can either be a list of times after initialization or an interval. 
   * - nhours_fcst
     - ``model_configure`` (uses ``FHMAX`` value set in the test file or ``default_vars.sh``)
     - 24
     - The maximal output time for the forecast.

=============================================================
How do I turn off IO for the components of the coupled model?
=============================================================

FV3atm restart and history files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To turn off FV3atm restart files, set the ``restart_interval`` in
``model_configure*.IN`` to a value greater than the forecast length.

To turn off history files, in ``model_configure`` there are two
options:

* Set ``quilting`` to .false., then in ``diag_table``, remove the history
  output file definitions ``fv3_history`` and ``fv3_history2d`` and the
  associated fields. This will turn off the write_grid component and the
  number of tasks used by FV3atm must also be adjusted to remove the
  tasks assigned to the write grid component.

* Set ``quilting`` to .true., then in ``model_configure`` set
  ``write_dopost`` to .false. and set ``output_fh`` to a value greater
  than the forecast length. This will turn off the writing of output but
  the write grid component tasks will still be necessary.

MOM6, CICE6 and CMEPS restart files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In ``ufs.configure*.IN``, set the ALLCOMP_attribute ``restart_n`` to a
value greater than the forecast length.

MOM6 history files
^^^^^^^^^^^^^^^^^^

In the ``diag_table*`` file, remove the ``ocn`` and ``SST`` history
output file definitions and fields.

MOM6 history output speed can also be increased by setting the
``IO_LAYOUT`` parameter in the relevant ``parm/MOM_input*.IN`` file.

::

   IO_LAYOUT = 4,2

CICE history files
^^^^^^^^^^^^^^^^^^

In the CICE namelist ``ice_in.IN``, set the ``histfreq`` to none with

::

   histfreq = 'x','x','x','x','x'

The initial condition file can be turned off using

::

   write_ic = .false.

GOCART history files
^^^^^^^^^^^^^^^^^^^^

In ``parm/gocart/AERO_HISTORY.rc.IN``, remove all the fields listed in ``COLLECTIONS``.

WW3 history and restart files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In ``ww3_shel.inp``, change the output interval for gridded frequency from
3600 to 0 on `line 68
<https://github.com/NOAA-EMC/WW3/blob/f1f14d582835ef77626aaf44146a53d60920a0f7/model/inp/ww3_shel.inp#L68>`_. To
turn off point output, change the output frequency from 900 to 0 on
`line 298
<https://github.com/NOAA-EMC/WW3/blob/f1f14d582835ef77626aaf44146a53d60920a0f7/model/inp/ww3_shel.inp#L298>`_. To
turn off restart files, change the frequency from 3600 to 0 on `line
323
<https://github.com/NOAA-EMC/WW3/blob/f1f14d582835ef77626aaf44146a53d60920a0f7/model/inp/ww3_shel.inp#L323>`_.

==============================================================
How do I set the total number of tasks for my job?
==============================================================

In the UFS WM, each component's MPI task information, including the
starting and ending tasks and the number of threads, are specified
using the component-specific ``*petlist_bounds`` and
``*omp_num_threads`` in ``ufs.configure``. In general, the total
number of MPI tasks required is the sum of all the sub-component
tasks, as long as those components do not overlap (i.e., share the
same PETs). An example of a global five-component coupled configuration
``ufs.configure`` appears at the end of this section.

FV3atm
^^^^^^

The FV3atm component consists of one or more forecast grid components
and write grid components.

The MPI tasks for the forecast grid components are specified in the
layout variable in one or more namelist files ``input*.nml``
(e.g., ``input.nml`` and ``input_nest02.nml``). The total number of MPI tasks
required is given by the product of the specified layout, summed over
all domains. For example, for a global domain with 6 tiles and
``layout = 6,8``, the total number required is ``6*6*8 = 288``. For
two regional domains using ``input.nml`` and ``input_nest02.nml``,
each with ``layout = 6,10``, the total required is the sum ``6*10 +
6*10 = 120``.

For the global configuration, an additional requirement is that the
layout specified must be a multiple of the ``blocksize`` parameter in
``input.nml``.  For example, using ``layout=8,8`` for C96 yields
subdomains of ``12 x 12``. The subdomain product is ``12*12 = 144``,
which is not divisible by a ``blocksize=32``. Therefore, the C96 does
not support an ``8,8`` layout for a blocksize of 32. If ``layout =
4,6``, the subdomain product is ``24*16 = 384``, which is divisible by
a ``blocksize=32``. A layout of ``4,6`` is supported for C96 with a
blocksize of 32.

The FV3atm will utilize the write grid component if ``quilting`` is
set to .true. In this case, the required MPI tasks for the
write grid component are the product of the ``write_groups`` and the
``write_tasks_per_group`` in the ``model_configure`` file.

::

   quilting:                .true.
   write_groups:            1
   write_tasks_per_group:   60


In the above case, the write grid component requires 60 tasks.

The total number of MPI ranks for FV3atm is the sum of the forecast tasks and any
write grid component tasks.

::

   total_tasks_atm = forecast tasks +  write grid component tasks

If ESMF-managed threading is used, the total number of PETs for the
atmosphere component is given by the product of the number of threads
requested and the total number of MPI ranks (both forecast and write
grid component). If ``num_threads_atm`` is the number of threads
specified for the FV3atm component, in ``ufs.configure`` the ATM PET
bounds are given by:

::

   ATM_petlist_bounds     0 total_tasks_atm*num_threads_atm-1
   ATM_omp_num_threads    num_threads_atm

Note that in UFS WM, the ATM component is normally listed first in
``ufs.configure`` so that the starting PET for the ATM is 0.

GOCART
^^^^^^

GOCART shares the same grid and forecast tasks as FV3atm, but it does
not have a separate write grid component in its NUOPC CAP. Also, while
GOCART does not have threading capability, it shares the same data
structure as FV3atm and so it has to use the same number of threads
used by FV3atm. Therefore, the total number of MPI ranks and threads
in GOCART is the same as the those for the FV3atm forecast component
(i.e., excluding any write grid component). Currently, GOCART only runs
on the global forecast grid component, for which only one namelist is
needed.

::

   total_tasks_chm = FV3atm forecast tasks

   CHM_petlist_bounds:             0 total_tasks_chm*num_threads_atm-1
   CHM_omp_num_threads:            num_threads_atm

CMEPS
^^^^^

The mediator MPI tasks can overlap with other components and in UFS
the tasks are normally shared on the FV3atm forecast tasks. However, a
large number of tasks for the mediator is generally not recommended
since it may cause slow performance. This means that the number of
MPI tasks for CMEPS is given by

::

   total_tasks_med = smaller of (300, FV3atm forecast tasks)

and in ``ufs.configure``

::

   MED_petlist_bounds:             0 total_tasks_med*num_threads_atm-1
   MED_omp_num_threads:            num_threads_atm

MOM6
^^^^

For MOM6 the only restriction currently on the number of MPI ranks
used by MOM6 is that it is divisible by 2. The starting PET in
``ufs.configure`` will be the last PET of the preceding component,
incremented by one. Threading in MOM6 is not recommended at this time.

::

   OCN_petlist_bounds:             starting_OCN_PET  total_tasks_ocn+starting_OCN_PET-1
   OCN_omp_num_threads:            1

CICE
^^^^

CICE requires setting the decomposition shape, the number of requested
processors and the calculated block sizes in the ``ice_in``
namelist. In UFS, the decomposition shape is always ``SlenderX2``,
except for the 5-degree configuration, which is ``SlenderX1``.

For ``SlenderX2`` decomposition, a given ``nprocs``, and global domain
``nx_global``, ``ny_global``, the block sizes are given by:

::

  block_size_y = ny_global/2
  block_size_x = nx_global/(nprocs/2)

Similarily, for ``SlenderX1``:

::

   block_size_y = ny_global
   block_size_x = nx_global/nprocs


For the 1-degree CICE domain for example, ``ice_in`` would be:

::

    nprocs            = 10
    nx_global         = 360
    ny_global         = 320
    block_size_x      = 72
    block_size_y      = 160
    max_blocks        = -1
    processor_shape   = 'slenderX2'


In the UFS, only a single thread is used for CICE so for ``nprocs`` set in
``ice_in``, the tasks in ``ufs.configure`` are set as:

::

   ICE_petlist_bounds:            starting_ICE_PET  nprocs+starting_ICE_PET-1
   ICE_omp_num_threads:           1

The starting ICE PET in ``ufs.configure`` will be the last PET of the
preceding component, incremented by one.

WW3
^^^

The WW3 component requires setting only the MPI ranks available
for WW3 and the number of threads to be used.

::

   WAV_petlist_bounds:         starting_WAV_PET  num_tasks_wav*num_threads_wav+starting_WAV_PET-1
   WAV_omp_num_threads:        num_threads_wav

The starting WAV PET in ``ufs.configure`` will be the last PET of the
preceding component, incremented by one.


Example: 5-component ufs.configure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A sample ``ufs.configure`` is shown below for the :wm-repo:`cpld_control_gefs <blob/develop/tests/tests/cpld_control_gefs>` test, which is a fully coupled S2SWA run. This test uses the :wm-repo:`ufs.configure.s2swa.IN <blob/develop/tests/parm/ufs.configure.s2swa.IN>` template to generate ``ufs.configure``.


.. code-block:: console

		#############################################
		####  UFS Run-Time Configuration File  #####
		#############################################

		# ESMF #
		logKindFlag:            ESMF_LOGKIND_MULTI
		globalResourceControl:  true

		# EARTH #
		EARTH_component_list: MED ATM CHM OCN ICE WAV
		EARTH_attributes::
		   Verbosity = 0
		::

		# MED #
		MED_model:                      cmeps
		MED_petlist_bounds:             0 599
		MED_omp_num_threads:            2
		::


		# ATM #
		ATM_model:                      fv3
		ATM_petlist_bounds:             0 959
		ATM_omp_num_threads:            2
		ATM_attributes::
         Verbosity = 0
         DumpFields = false
         ProfileMemory = false
         OverwriteSlice = true
		::

		# CHM #
		CHM_model:                      gocart
		CHM_petlist_bounds:             0 767
		CHM_omp_num_threads:            2
		CHM_attributes::
         Verbosity = 0
		::

		# OCN #
		OCN_model:                      mom6
		OCN_petlist_bounds:             960 1079
		OCN_omp_num_threads:            1
		OCN_attributes::
		   Verbosity = 0
		   DumpFields = false
		   ProfileMemory = false
		   OverwriteSlice = true
		   mesh_ocn = mesh.mx025.nc
         use_coldstart = false
         use_mommesh = true
		::

		# ICE #
		ICE_model:                      cice6
		ICE_petlist_bounds:             1080 1127
		ICE_omp_num_threads:            1
		ICE_attributes::
		   Verbosity = 0
		   DumpFields = false
		   ProfileMemory = false
		   OverwriteSlice = true
		   mesh_ice = mesh.mx025.nc
         eps_imesh = 1.0e-1
		   stop_n = 3
		   stop_option = nhours
		   stop_ymd = -999
		::

		# WAV #
		WAV_model:                      ww3
		WAV_petlist_bounds:             1128 1367
		WAV_omp_num_threads:            2
		WAV_attributes::
		   Verbosity = 0
		   OverwriteSlice = false
         mesh_wav = mesh.glo_025.nc
		   user_histname = true
         use_historync = true
         use_restartnc = true
         restart_from_binary = true
         pio_typename = pnetcdf
         pio_numiotasks = -99
         pio_stride = 4
         pio_rearranger = box
         pio_root = -99
		::

		 # CMEPS warm run sequence
		runSeq::
		@1800
         MED med_phases_prep_wav_avg
         MED med_phases_prep_ocn_avg
         MED -> WAV :remapMethod=redist
         MED -> OCN :remapMethod=redist
         WAV
         OCN
         @300
            MED med_phases_prep_atm
            MED med_phases_prep_ice
            MED -> ATM :remapMethod=redist
            MED -> ICE :remapMethod=redist
            ATM phase1
            ATM -> CHM
            CHM
            CHM -> ATM
            ATM phase2
            ICE
            ATM -> MED :remapMethod=redist
            MED med_phases_post_atm
            ICE -> MED :remapMethod=redist
            MED med_phases_post_ice
            MED med_phases_ocnalb_run
            MED med_phases_prep_ocn_accum
            MED med_phases_prep_wav_accum
         @
         OCN -> MED :remapMethod=redist
         WAV -> MED :remapMethod=redist
         MED med_phases_post_ocn
         MED med_phases_post_wav
         MED med_phases_restart_write
		@
		::

		# CMEPS variables

		DRIVER_attributes::
		::

		MED_attributes::
		   ATM_model = fv3
         ICE_model = cice6
         OCN_model = mom6
         WAV_model = ww3
         coupling_mode = ufs.frac
         pio_rearranger = box
         ocean_albedo_limit = 0.06
		::
		ALLCOMP_attributes::
		   ScalarFieldCount = 3
         ScalarFieldIdxGridNX = 1
         ScalarFieldIdxGridNY = 2
         ScalarFieldIdxGridNTile = 3
         ScalarFieldName = cpl_scalars
         start_type = continue
         restart_dir = ./RESTART/
         case_name = ufs.cpld
         restart_n = 3
         restart_option = nhours
         restart_ymd = -999
         write_restart_at_endofrun = .false.
         dbug_flag = 0
         stop_n = 12
         stop_option = nhours
         stop_ymd = -999
         orb_eccen = 1.e36
         orb_iyear = 2000
         orb_iyear_align = 2000
         orb_mode = fixed_year
         orb_mvelp = 1.e36
         orb_obliq = 1.e36
		::

========================================================
How can I get the UFS WM to output physics tendencies?
========================================================

Users will need to:

#. Update ``input.nml`` by setting ``ldiag3d`` and ``qdiag3d`` to ``.true.``. 
#. Update the ``diag_table`` according to the instructions in :numref:`Section %s <diag_tableFile>`.
 
Although it may seem counterintuitive, the physics tendencies will be output in ``sfc*.nc`` files once the ``diag_table`` changes have been made. Even 3D fields will appear there. 

Users may find the following GitHub Discussions on this topic informative: 

* :wm-repo:`WM Discussion #1867<discussions/1867>` 
* `SRW App Discussion #862 <https://github.com/ufs-community/ufs-srweather-app/discussions/862>`_ 
* :wm-repo:`WM Discussion #1862 <discussions/1862>`

===================================================================================================================
How can I output a particular variable (e.g., accumulated precipitation) from the UFS WM atmospheric model (FV3)?
===================================================================================================================

To output a particular variable from FV3, users must update the field section of the ``diag_table`` file, which specifies the fields to be output at run time. 
Only fields registered with ``register_diag_field()``, which is an API in the FMS ``diag_manager`` routine, can be used in the ``diag_table``. 
A line in the field section of the ``diag_table`` file contains eight variables with the following format:

.. code-block:: console

   "module_name", "field_name", "output_name", "file_name", "time_sampling", "reduction_method", "regional_section", packing

These variables are defined in :numref:`Table %s <diag-table-options>` of the UFS WM documentation on the ``diag_table`` file. 

For example, to output accumulated precipitation, the following line must appear in the ``diag_table`` file: 

.. code-block:: console

   "gfs_phys", "totprcp_ave", "prate_ave", "fv3_history2d", "all", .false., "none", 2

Users may refer to ``diag_table`` examples in the UFS WM repository. These files are used to configure groups of regression tests. 

View GitHub :wm-repo:`Discussion #2016 <discussions/2016/>` for the question that inspired this FAQ. 

===========================================================================================================
Where can I find up-to-date documentation for the ``diag_table`` variables used in the UFS Weather Model?
===========================================================================================================

Information on ``diag_table`` variables has been added to the :ref:`diag_table section <diag-table-options>` of the UFS Weather Model documentation. 
Currently, only variables coming from fv3atm and MOM6 are included, but ``diag_table`` variables from other components will be added as time permits. 

* :ref:`FV3ATM diag_table variables <fv3diagtable>`
* `MOM6 diag_table variables <https://ncar.github.io/MOM6/APIs/namespacemom__diagnostics.html>`_

See ufs-community `Discussion #33 <https://github.com/orgs/ufs-community/discussions/33>`_ for the question that inspired this FAQ.


