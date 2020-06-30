.. _FAQ:

***
FAQ
***

==============================================================
How do I build and run a single test of the UFS Weather Model?
==============================================================

An efficient way to build and run the UFS Weather Model is to use the regression test
(``rt.sh``).  This script is widely used by model developers on Tier 1 and 2 platforms
and is described in :numref:`Section %s <ConductingRegTests>`.  The advantages to this approach are:

- You can bypass the workflow, pre- and post-processing steps.
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

=====================================================
How do I change the model output to netcdf or nemsio?
=====================================================
In your run directory, there is a file named ``model_configure``.  Change the
variable ``output_file`` to ``'netcdf'`` or ``'nemsio'``.

