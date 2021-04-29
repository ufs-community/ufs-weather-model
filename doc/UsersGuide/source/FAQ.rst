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

- It does not require a workflow, pre- or post-processing steps.
- The batch submission script is generated.
- Any required input data is already available for machines used by the regression test.
- Once the ``rt.sh`` test completes, you will have a working copy in your run directory
  where you can make modifications to the namelist and other files, and then re-run the
  executable.

The steps are:

1. Clone the source code and all the submodules as described in :numref:`Section %s <DownloadingWMCode>`, then
   go into the ``tests`` directory:

   .. code-block:: console

       cd ufs-weather-model (or the top level where you checked out the code)
       cd tests

2. Find a configure (``*.conf``) file that contains the machine and compiler you are using. For this
   example, the Intel compiler on Cheyenne is used.  To create a custom configure file, two lines are
   needed:  a ``COMPILE`` line and a ``RUN`` line.   The ``COMPILE`` line should contain the name
   of the machine and compiler ``cheyenne.intel`` and the desired ``SUITES`` for the build.  Choose a
   ``RUN`` line under this ``COMPILE`` command that uses the desired ``SUITE``.  For example:

   .. code-block:: console

       COMPILE | 32BIT=Y CCPP=Y STATIC=Y SUITES=FV3_GFS_v15p2,FV3_GFS_v16beta,FV3_GFS_v15p2_no_nsst,FV3_GFS_v16beta_no_nsst                     | standard    | cheyenne.intel | fv3
       RUN     | fv3_ccpp_gfs_v16beta                                                                                                           | standard    |                | fv3         |

   Put these two lines into a file called ``my_test.conf``.  The parameters used in this run can be
   found in the ``fv3_ccpp_gfs_v16beta`` file in the ``ufs-weather-model/tests/tests`` directory.

   .. note::  These two lines are long and may not appear in entirety in your browser. Scroll to the right to see
              the entire line.

3. Modify the ``rt.sh`` script to put the output in a run directory where you have write permission:

   .. code-block:: console

       if [[ $MACHINE_ID = cheyenne.* ]]; then stanza:
       ...
       dprefix=/glade/scratch

   This works for Cheyenne, since ``$USER/FV3_RT`` will be appended.  Also check that ``RTPWD``
   points to a diretory that exists:

   .. code-block:: console

       if [[ $MACHINE_ID = cheyenne.* ]]; then
         RTPWD=${RTPWD:-$DISKNM/ufs-public-release-20200224/${COMPILER^^}}

4. Run the ``rt.sh`` script from the ``tests`` directory:

   .. code-block:: console

       ./rt.sh -k -l my_test.conf >& my_test.out &

   Check ``my_test.out`` for build and run status, plus other standard output. Check
   ``/glade/scratch/$USER/FV3_RT/rt_PID`` for the model run, where ``PID`` is a process ID.
   The build will take about 10-15 minutes and the run will be fast, depending on how long
   it waits in the queue.  A message ``"REGRESSION TEST WAS SUCCESSFUL"`` will be written to this
   file, along with other entertainment: ``'Elapsed time: 00h:14m:12s. Have a nice day!'``.

5. When the build and run are complete, modify the namelist or ``model_configure`` files
   and re-run by submitting the ``job_card`` file:

   .. code-block:: console

       qsub job_card

============================================
How do I change the length of the model run?
============================================
In your run directory, there is a file named ``model_configure``.  Change the
variable ``nhours_fcst`` to the desired number of hours.

========================================================================
How do I select the file format for the model output (netCDF or NEMSIO)?
========================================================================
In your run directory, there is a file named ``model_configure``.  Change the
variable ``output_file`` to ``'netcdf'`` or ``'nemsio'``. The variable ``output_file``
is only valid when the write component is activated by setting ``quilting`` to .true.
in the ``model_configure`` file.

==============================================================
How do I set the output history interval?
==============================================================
The interval at which output (history) files are written is controlled in two
places, and depends on whether you are using the write component to generate your output files.
:numref:`Table %s <OutputControl>` describes the relevant variables.  If the write_component is used, then the variables listed as *model_configure* are required.  It is however, also required that the settings in *input.nml* match those same settings in *model_configure*.  If these settings are inconsistent, then unpredictable output files and intervals may occur!

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

==============================================================
How do I set the total number of tasks for my job?
==============================================================
The total number of MPI tasks used by the UFS Weather Model is a combination of compute and quilt tasks, and can be calculated using the following relationship:

- total tasks = compute tasks + quilt tasks
- compute tasks = x layout * y layout * number of tiles
- quilt tasks = write_groups * write_tasks_per_group if quilting==.true.

The layout and tiles settings are in ``input.nml``, and the quilt task settings are in ``model_configure``
