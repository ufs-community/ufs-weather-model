.. _FAQ:

***
FAQ
***

==============================================================
How do I build and run a single test of the UFS Weather Model?
==============================================================

An efficient way to build and run the UFS Weather Model is to use the regression test
(``rt.sh``).  This script is widely used by model developers on Tier 1 and 2 platforms
and is described in the UFS WM GitHub `wiki <https://github.com/ufs-community/ufs-weather-model/wiki/Making-code-changes-in-the-UFS-weather-model-and-its-subcomponents>`_.  The advantages to this approach are:

   * It does not require a workflow, pre- or post-processing steps.
   * The batch submission script is generated.
   * Any required input data is already available for machines used by the regression test.
   * Once the ``rt.sh`` test completes, you will have a working copy in your run directory where you can
     make modifications to the namelist and other files, and then re-run the executable.

The steps are:

   #. Clone the source code and all the submodules as described in :numref:`Section %s <DownloadingWMCode>`, then
      go into the ``tests`` directory:

      .. code-block:: console

         cd ufs-weather-model (or the top level where you checked out the code)
         cd tests

   #. Find a configure (``*.conf``) file that contains the machine and compiler you are using. For this
      example, the Intel compiler on Derecho is used.  To create a custom configure file, two lines are
      needed:  a ``COMPILE`` line and a ``RUN`` line.   The ``COMPILE`` line should contain the name
      of the machine and compiler ``derecho.intel`` and the desired ``SUITES`` for the build.  Choose a
      ``RUN`` line under this ``COMPILE`` command that uses the desired ``SUITE``.  For example:

      .. code-block:: console

         COMPILE | 32BIT=Y CCPP=Y STATIC=Y SUITES=FV3_GFS_v15p2,FV3_GFS_v16beta,FV3_GFS_v15p2_no_nsst,FV3_GFS_v16beta_no_nsst                     | standard    | derecho.intel | fv3
         RUN     | fv3_ccpp_gfs_v16beta                                                                                                           | standard    |                | fv3         |

      Put these two lines into a file called ``my_test.conf``.  The parameters used in this run can be
      found in the ``fv3_ccpp_gfs_v16beta`` file in the ``ufs-weather-model/tests/tests`` directory.

      .. note::  These two lines are long and may not appear in entirety in your browser. Scroll to the right to see
               the entire line.

   #. Modify the ``rt.sh`` script to put the output in a run directory where you have write permission:

      .. code-block:: console

         if [[ $MACHINE_ID = derecho.* ]]; then stanza:
         ...
         dprefix=/glade/scratch

      This works for Derecho, since ``$USER/FV3_RT`` will be appended.  Also check that ``RTPWD``
      points to a diretory that exists:

      .. code-block:: console

         if [[ $MACHINE_ID = derecho.* ]]; then
            RTPWD=${RTPWD:-$DISKNM/ufs-public-release-20200224/${COMPILER^^}}

   #. Run the ``rt.sh`` script from the ``tests`` directory:

      .. code-block:: console

         ./rt.sh -k -l my_test.conf >& my_test.out &

      Check ``my_test.out`` for build and run status, plus other standard output. Check
      ``/glade/scratch/$USER/FV3_RT/rt_PID`` for the model run, where ``PID`` is a process ID.
      The build will take about 10-15 minutes and the run will be fast, depending on how long
      it waits in the queue.  A message ``"REGRESSION TEST WAS SUCCESSFUL"`` will be written to this
      file, along with other entertainment: ``'Elapsed time: 00h:14m:12s. Have a nice day!'``.

   #. When the build and run are complete, modify the namelist or ``model_configure`` files
      and re-run by submitting the ``job_card`` file:

      .. code-block:: console

         qsub job_card

============================================
How do I change the length of the model run?
============================================
In your run directory, there is a file named ``model_configure``.  Change the
variable ``nhours_fcst`` to the desired number of hours.

==============================================================
How do I set the output history interval?
==============================================================

The interval at which output (history) files are written is controlled in two
places, and depends on whether you are using the write component to generate your output files.
:numref:`Table %s <OutputControl>` describes the relevant variables.  If the write_component is used, then the variables listed as ``model_configure`` are required.  It is however, also required that the settings in ``input.nml`` match those same settings in ``model_configure``.  If these settings are inconsistent, then unpredictable output files and intervals may occur!

.. _OutputControl:

.. list-table:: *Namelist variables used to control the output file frequency.*
   :widths: 15 10 10 30
   :header-rows: 1

   * - Namelist variable
     - Location
     - Default Value
     - Description
   * - fdiag
     - input.nml
     - 0
     - Array with dimension ``maxhr`` = 4096 listing the diagnostic output times (in hours) for the GFS physics.
       This can either be a list of times after initialization, or an interval if only the first entry is
       nonzero. The default setting of 0 will result in no outputs.
   * - fhmax
     - input.nml
     - 384
     - The maximal forecast time for output.
   * - fhmaxhf
     - input.nml
     - 120
     - The maximal forecast hour for high frequency output.
   * - fhout
     - input.nml
     - 3
     - Output frequency during forecast time from 0 to ``fhmax``, or from ``fhmaxhf`` to ``fhmax`` if ``fhmaxf>0``.
   * - fhouthf
     - input.nml
     - 1
     - The high frequency output frequency during the forecast time from 0 to ``fhmaxhf`` hour.
   * - nfhmax_hf
     - model_configure
     - 0
     - forecast length of high history file
   * - nfhout_hf
     - model_configure
     - 1
     - high history file output frequency
   * - nfhout
     - model_configure
     - 3
     - history file output frequency

=============================================================
How do I turn off IO for the components of the coupled model?
=============================================================

FV3atm restart and history files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To turn off FV3atm restart files, set the ``restart_interval`` in
``model_configure`` to a value greater than the forecast length.

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

In ``ufs.configure``, set the ALLCOMP_attribute ``restart_n`` to a
value greater than the forecast length.

MOM6 history files
^^^^^^^^^^^^^^^^^^

In the ``diag_table`` file, remove the ``ocn`` and ``SST`` history
output file definitions and fields.

MOM6 history output speed can also be increased by setting the
``IO_LAYOUT`` parameter in ``INPUT/MOM_input``.

::

   IO_LAYOUT = 4,2

CICE history files
^^^^^^^^^^^^^^^^^^

In the CICE namelist ``ice_in``, set the ``histfreq`` to none with

::

   histfreq = 'x','x','x','x','x'

The initial condition file can be turned off using

::

   write_ic = .false.

GOCART history files
^^^^^^^^^^^^^^^^^^^^

In AERO_HISTORY.rc, remove all the fields listed in ``COLLECTIONS``

::

   COLLECTIONS:
   ::

WW3 history and restart files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In ``ww3_shel.inp``, change the output interval for gridded frequency from
3600 to 0 on `line 68
<https://github.com/NOAA-EMC/WW3/blob/5ebed915755da0b21cf4d20e21726411fb2948c4/model/inp/ww3_shel.inp#L68>`_. To
turn off point output, change the output frequency from 900 to 0 on
`line 296
<https://github.com/NOAA-EMC/WW3/blob/5ebed915755da0b21cf4d20e21726411fb2948c4/model/inp/ww3_shel.inp#L296>`_. To
turn off restart files, change the frequency from 3600 to 0 on `line
321
<https://github.com/NOAA-EMC/WW3/blob/5ebed915755da0b21cf4d20e21726411fb2948c4/model/inp/ww3_shel.inp#L321>`_.



==============================================================
How do I set the total number of tasks for my job?
==============================================================

In the UFS WM, each component's MPI task information, including the
starting and ending tasks and the number of threads, are specified
using the component-specific ``petlist_bounds`` and
``omp_num_threads`` in ``ufs.configure``. In general, the total
number of MPI tasks required is the sum of all the sub-component
tasks, as long as those components do not overlap (i.e., share the
same PETs). An example of a global 5 component coupled configuration
ufs.configure at the end of this section.

FV3atm
^^^^^^

The FV3atm component consists of one or more forecast grid components
and write grid components.

The MPI tasks for the forecast grid components are specified in the
layout variable in one or more namelist files ``input*.nml``
(e.g. input.nml and input_nest02.nml). The total number of mpi tasks
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
set to .true. In this case, the required mpi tasks for the
write grid components is the product of the ``write_groups`` and the
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
bounds are given by

::

   ATM_petlist_bounds     0 total_tasks_atm*num_threads_atm-1
   ATM_omp_num_threads    num_threads_atm

Note that in UWM, the ATM component is normally listed first in
``ufs.configure`` so that the starting PET for the ATM is 0.

GOCART
^^^^^^

GOCART shares the same grid and forecast tasks as FV3atm but it does
not have a separate write grid component in its NUOPC CAP. Also, while
GOCART does not have threading capability, it shares the same data
structure as FV3atm and so it has to use the same number of threads
used by FV3atm. Therefore, the total number of MPI ranks and threads
in GOCART is the same as the those for the FV3atm forecast component
(i.e., excluding any write grid component). Currently GOCART only runs
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
except for the 5 deg configuration, which is ``SlenderX1``.

For ``SlenderX2`` decomposition, a given ``nprocs``, and global domain
``nx_global``, ``ny_global``, the block sizes are given by

::

  block_size_y = ny_global/2
  block_size_x = nx_global/(nprocs/2)

Similarily, for ``SlenderX1``

::

   block_size_y = ny_global
   block_size_x = nx_global/nprocs


For the 1-deg CICE domain for example, ``ice_in`` would be

::

    nprocs            = 10
    nx_global         = 360
    ny_global         = 320
    block_size_x      = 72
    block_size_y      = 160
    max_blocks        = -1
    processor_shape   = 'slenderX2'


In UFS, only a single thread is used for CICE so for ``nprocs`` set in
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

For the fully coupled S2SWA application, a sample ``ufs.configure`` is shown below :


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
		MED_petlist_bounds:             0 767
		MED_omp_num_threads:            2
		::


		# ATM #
		ATM_model:                      fv3
		ATM_petlist_bounds:             0 863
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
		 OCN_petlist_bounds:             864 983
		 OCN_omp_num_threads:            1
		 OCN_attributes::
		   Verbosity = 0
		   DumpFields = false
		   ProfileMemory = false
		   OverwriteSlice = true
		   mesh_ocn = mesh.mx025.nc
		 ::

		 # ICE #
		 ICE_model:                      cice6
		 ICE_petlist_bounds:             984 1031
		 ICE_omp_num_threads:            1
		 ICE_attributes::
		   Verbosity = 0
		   DumpFields = false
		   ProfileMemory = false
		   OverwriteSlice = true
		   mesh_ice = mesh.mx025.nc
		   stop_n = 3
		   stop_option = nhours
		   stop_ymd = -999
		 ::

		 # WAV #
		 WAV_model:                      ww3
		 WAV_petlist_bounds:             1032 1191
		 WAV_omp_num_threads:            2
		 WAV_attributes::
		   Verbosity = 0
		   OverwriteSlice = false
		   diro = "."
		   logfile = wav.log
		   mesh_wav = mesh.gwes_30m.nc
		   multigrid = false
		 ::

		 CMEPS warm run sequence
		 runSeq::
		 @1800
		 MED med_phases_prep_ocn_avg
		 MED -> OCN :remapMethod=redist
		 OCN
		 @300
		   MED med_phases_prep_atm
		   MED med_phases_prep_ice
		   MED med_phases_prep_wav_accum
		   MED med_phases_prep_wav_avg
		   MED -> ATM :remapMethod=redist
		   MED -> ICE :remapMethod=redist
		   MED -> WAV :remapMethod=redist
		   ATM phase1
		   ATM -> CHM
		   CHM
		   CHM -> ATM
		   ATM phase2
		   ICE
		   WAV
		   ATM -> MED :remapMethod=redist
		   MED med_phases_post_atm
		   ICE -> MED :remapMethod=redist
		   MED med_phases_post_ice
		   WAV -> MED :remapMethod=redist
		   MED med_phases_post_wav
		   MED med_phases_prep_ocn_accum
		 @
		 OCN -> MED :remapMethod=redist
		 MED med_phases_post_ocn
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
		  history_n = 1
		  history_option = nhours
		  history_ymd = -999
		  coupling_mode = nems_frac
		  history_tile_atm = 384
		::
		ALLCOMP_attributes::
		  ScalarFieldCount = 2
		  ScalarFieldIdxGridNX = 1
		  ScalarFieldIdxGridNY = 2
		  ScalarFieldName = cpl_scalars
		  start_type = startup
		  restart_dir = RESTART/
		  case_name = ufs.cpld
		  restart_n = 3
		  restart_option = nhours
		  restart_ymd = -999
		  dbug_flag = 0
		  use_coldstart = false
		  use_mommesh = true
		  eps_imesh = 1.0e-1
		  stop_n = 6
		  stop_option = nhours
		  stop_ymd = -999
		::
