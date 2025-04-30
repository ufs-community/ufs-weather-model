.. role:: raw-html(raw)
    :format: html

.. _hsd-container:

**********************************************
Running UFS WM Idealized Cases in a Container
**********************************************

This chapter provides instructions for running the Unified Forecast System (:term:`UFS`) Weather Model (WM) Hierarchical System Develop (HSD) cases using a `Singularity/Apptainer <https://apptainer.org/docs/user/latest/>`_ container. Normally, the details of building and running Earth system models will vary based on the computing platform because there are many possible combinations of operating systems, compilers, :term:`MPIs <MPI>`, and package versions available. Installation via Singularity/Apptainer container reduces this variability and allows for a smoother experience building and running the UFS WM. This approach is recommended for users not running the UFS WM on a supported :wm-wiki:`Level 1 <Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>` system (e.g., Hera, Orion). 

This chapter provides instructions for building and running the Unified Forecast System UFS WM HSD cases using a container. Currently, users can select from the following cases: 

   * The :ref:`July 2020 CAPE Case <cape-2020>`
   * The :ref:`Baroclinic Instability Case <baroclinic-wave>`

.. attention::

   This chapter of the User's Guide should **only** be used for container builds. For non-container builds, see the chapters above for instructions on how to run individual cases. These chapters describe the steps for configuring and running the UFS WM HSD cases on a :wm-wiki:`Level 1 System <Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>` **without** a container.

.. _Prereqs:

Prerequisites 
*****************

The containerized version of the UFS WM requires: 

   * `Installation of Apptainer <https://apptainer.org/docs/admin/latest/installation.html>`_ (or its predecessor, Singularity)
   * At least 240 CPU cores with individual nodes capable of running more than 40 tasks
   * An **Intel** compiler and :term:`MPI` (available for `free here <https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html>`_) 
   * The `Slurm <https://slurm.schedmd.com/quickstart.html>`_ job scheduler


Install Apptainer
==================

.. note::

   As of November 2021, the Linux-supported version of Singularity has been `renamed <https://apptainer.org/news/community-announcement-20211130/>`_ to *Apptainer*. Apptainer has maintained compatibility with Singularity, so ``singularity`` commands should work with either Singularity or Apptainer (see `compatibility details here <https://apptainer.org/docs/user/1.2/introduction.html>`_.)

To configure and run the HSD cases using a Apptainer container, first install the software according to the `Apptainer Installation Guide <https://apptainer.org/docs/admin/1.2/installation.html>`_. This will include the installation of all dependencies. 

.. attention:: 
   Docker containers can only be run with root privileges, and users generally do not have root privileges on :term:`HPCs <HPC>`. However, an Apptainer image may be built directly from a Docker image for use on the system.

.. _DownloadContainer:

Build the Container
**********************

.. _CloudHPC:

Set Environment Variables
=============================

For users working on systems with limited disk space in their ``/home`` directory, it is important to set the ``SINGULARITY_CACHEDIR`` and ``SINGULARITY_TMPDIR`` environment variables to point to a location with adequate disk space. For example:

.. code-block:: 

   export SINGULARITY_CACHEDIR=/absolute/path/to/writable/directory/cache
   export SINGULARITY_TMPDIR=/absolute/path/to/writable/directory/tmp

where ``/absolute/path/to/writable/directory/`` refers to a writable directory (usually a project or user directory within ``/lustre``, ``/work``, ``/scratch``, or ``/glade`` on NOAA :term:`RDHPCS` systems). If the ``cache`` and ``tmp`` directories do not already exist, they must be created with the ``mkdir`` command. 

On NOAA Cloud systems, the ``sudo su`` command may also be required. For example, users would run:
   
.. code-block:: 

   mkdir /lustre/cache
   mkdir /lustre/tmp
   sudo su
   export SINGULARITY_CACHEDIR=/lustre/cache
   export SINGULARITY_TMPDIR=/lustre/tmp
   exit

.. note:: 
   ``/lustre`` is a fast but non-persistent file system used on NOAA Cloud systems. To retain work completed in this directory, `tar the files <https://www.howtogeek.com/248780/how-to-compress-and-extract-files-using-the-tar-command-on-linux/>`_ and move them to the ``/contrib`` directory, which is much slower but persistent.

.. COMMENT:

.. _ContainerBuild:

Build the Container
======================

Set a top-level directory location for UFS WM work, and navigate to it. For example:

.. code-block:: console 

   mkdir /path/to/hsd
   cd /path/to/hsd
   export HSD=`pwd`

where ``/path/to/hsd`` is the path to this top-level directory (e.g., ``/Users/Joe.Schmoe/hsd``). 

.. hint::
   If a ``singularity: command not found`` error message appears in any of the following steps, try running: ``module load singularity`` or ``module load apptainer``.

NOAA RDHPCS Systems
----------------------

On many NOAA :term:`RDHPCS`, a container named ``ubuntu22.04-intel-wm-dev-hsd-test.img`` has already been built, and users may access the container at the locations in :numref:`Table %s <PreBuiltContainers>`.

.. _PreBuiltContainers:

.. table:: Locations of Pre-Built Containers

   +--------------------+--------------------------------------------------------+
   | Machine            | File location                                          |
   +====================+========================================================+
   | Gaea               | /gpfs/f5/epic/world-shared/containers                  |
   +--------------------+--------------------------------------------------------+
   | Hera               | /scratch1/NCEPDEV/nems/role.epic/containers            |
   +--------------------+--------------------------------------------------------+
   | Jet                | /mnt/lfs5/HFIP/hfv3gfs/role.epic/containers            |
   +--------------------+--------------------------------------------------------+
   | NOAA Cloud [#fn]_  | /contrib/EPIC/containers                               |
   +--------------------+--------------------------------------------------------+
   | Orion/Hercules     | /work/noaa/epic/role-epic/contrib/containers           |
   +--------------------+--------------------------------------------------------+

.. [#fn] The CAPE case can run on the NOAA Cloud ParallelWorks (PW) platforms, but the baroclinic wave case cannot.

Users can simply set an environment variable to point to the container: 

.. code-block:: console

   export img=path/to/ubuntu22.04-intel-wm-dev-hsd-test.img

If users prefer, they may copy the container to their local working directory. For example, on Jet:

.. code-block:: console

   cp /mnt/lfs5/HFIP/hfv3gfs/role.epic/containers/ubuntu22.04-intel-wm-dev-hsd-test.img .

Other Systems
----------------

On other systems, users can build the Singularity container from a public Docker :term:`container` image or download the ``ubuntu22.04-intel-wm-dev-hsd-test.img`` container from the `UFS Hierarchical Testing Framework (HTF) Data Bucket <https://registry.opendata.aws/noaa-ufs-htf-pds/>`_. Downloading may be faster depending on the download speed on the user's system. Note that the container in the data bucket is from the November 20, 2024 ``develop`` branch.

To download from the data bucket, users can run:

.. code-block:: console

   wget https://noaa-ufs-htf-pds.s3.amazonaws.com/develop-20241115/ubuntu22.04-intel-wm-dev-hsd-test.img

To build the container from a Docker image, users can run:

.. code-block:: console

   singularity build --force ubuntu22.04-intel-wm-dev-hsd-test.img docker://noaaepic/ubuntu22.04-intel21.10-wm:ue160-fms202401-dev

This process may take several hours depending on the system. 

.. note:: 

   Some users may need to issue the ``singularity build`` command with ``sudo`` (i.e., ``sudo singularity build...``). Whether ``sudo`` is required is system-dependent. If ``sudo`` is required (or desired) for building the container, users should set the ``SINGULARITY_CACHEDIR`` and ``SINGULARITY_TMPDIR`` environment variables with ``sudo su``, as in the NOAA Cloud example from :numref:`Section %s <CloudHPC>` above.

.. _GetDataC:

Get Data
***********

In order to run the UFS WM HSD cases, users will need both fix files and model input data. These files are already present on :wm-wiki:`Level 1 <Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>` systems. 

Users on any system may download and untar the data from the `UFS Hierarchical Testing Framework (HTF) Data Bucket <https://registry.opendata.aws/noaa-ufs-htf-pds/>`_ into their ``$HSD`` directory. 

.. code-block:: console

   cd $HSD
   wget https://noaa-ufs-htf-pds.s3.amazonaws.com/develop-20241115/HSD_fix_files_and_case_data.tar.gz
   tar xvfz HSD_fix_files_and_case_data.tar.gz

.. _RunContainer:

Run the Container
********************

To run the container, users must:

   #. :ref:`Set up the container <SetUpContainer>`
   #. :ref:`Configure the experiment <ConfigureExptC>`
   #. :ref:`Run the experiment <RunExptC>`

.. _SetUpContainer:

Set Up the Container
=======================

Save the location of the container in an environment variable.

.. code-block:: console

   export img=/path/to/ubuntu22.04-intel-wm-dev-hsd-test.img

Users may convert a container ``.img`` file to a writable sandbox. This step is optional and unnecessary on most systems (it can take several hours):

.. code-block:: console

   singularity build --sandbox ubuntu22.04-intel-wm-dev-hsd-test $img

When making a writable sandbox on NOAA :term:`RDHPCS`, the following warnings commonly appear and can be ignored:

.. code-block:: console

   INFO:    Starting build...
   INFO:    Verifying bootstrap image ubuntu22.04-intel-wm-dev-hsd-test.img
   WARNING: integrity: signature not found for object group 1
   WARNING: Bootstrap image could not be verified, but build will continue.

From within the ``$HSD`` directory, copy the ``stage-rt.sh`` script out of the container. 

.. code-block:: console

   singularity exec -H $PWD $img cp /opt/stage-rt.sh .

The ``stage-rt.sh`` script should now be in the ``$HSD`` directory. If for some reason, the previous command was unsuccessful, users may try a version of the following command instead: 

.. code-block:: console

   singularity exec -B /<local_base_dir>:/<container_dir> $img cp /opt/stage-rt.sh .

where ``<local_base_dir>`` and ``<container_dir>`` are replaced with a top-level directory on the local system and in the container, respectively. Additional directories can be bound by adding another ``-B /<local_base_dir>:/<container_dir>`` argument before the container location (``$img``). Note that if previous steps included a ``sudo`` command, ``sudo`` may be required in front of this command. 

.. note::

   Sometimes binding directories with different names can cause problems. In general, it is recommended that the local base directory and the container directory have the same name. For example, if the host system's top-level directory is ``/user1234``, the user may want to convert the ``.img`` file to a writable sandbox and create a ``user1234`` directory in the sandbox to bind to. 

Run the ``stage-rt.sh`` script with the proper arguments. 

.. code-block:: console

   ./stage-rt.sh -c=<compiler> -m=<mpi_implementation> [-p=<platform>] -i=$img

where:

   * ``-c`` is the compiler on the user's local machine (e.g., ``intel/2022.1.2``)
   * ``-m`` is the :term:`MPI` on the user's local machine (e.g., ``impi/2022.1.2``)
   * ``-p`` refers to the local machine/platform (e.g., ``hera``, ``jet``, ``gaea``, ``noaacloud``). Required for Gaea and Jet only. 
   * ``-i`` is the full path to the container image (e.g., ``$img`` or ``$HSD/ubuntu22.04-intel-wm-dev-hsd-test.img``).

.. note::

   When using a Singularity container, Intel compilers and Intel :term:`MPI` (preferably 2020 versions or newer) need to be available on the host system to properly launch MPI jobs. Generally, this is accomplished by loading a module with a recent Intel compiler and then loading the corresponding Intel MPI. 

When this command runs, ``stage-rt.sh`` will print the following message to the console: 

.. code-block:: console

   Copying out ufs-weather-model repo from the container
   Set run_test.sh to use exe in the container
   Updating compiler and mpi in fv3_slurm.IN_singularity
   Creating ufs_singularity.intel.lua
   Tricking ufs_test.sh file
   Updating various files with host paths
   Done

Additionally, the user should see the ``ufs-weather-model`` directory in the ``$HSD`` directory (``ls``). 

.. note::

   Gaea and Jet:
      * Gaea uses a different compiler and MPI to run with the container: ``-c=intel-classic/2023.2.0 -m=cray-mpich/8.1.28``
      * On Jet, ``cd`` to ``/mnt`` first before navigating to individual user workspaces to use the container.

.. _ConfigureExptC:

Configure the Experiment
===========================

To configure the experiment, users may need to update the ``default_vars.sh`` script and/or the ``machine_singularity.config`` files. 

Module Modification
--------------------

The machine configuration file is located at ``ufs-weather-model/tests-dev/machine_config/machine_singularity.config``. It assumes that Rocoto can be loaded via ``module load`` command from the host machine's initial state. If an additional path or module needs to be loaded, modify the ``machine_singularity.config`` to reflect those additions. For example, if the Rocoto package is found within the ``contrib`` module, add ``module load contrib`` before the ``module load rocoto`` statement in the machine configuration file.

Host Machine Modifications
---------------------------

Default variables for regression tests and HSD tests are set in the ``default_vars.sh`` script in the ``ufs-weather-model/tests`` directory copied *from the container*. The individual test scripts (e.g., ``baroclinic_wave``, ``2020_CAPE``) override these variables where necessary. However, when running the HSD cases in a container, the tasks-per-node (TPN) variables in the singularity section need to be modified to reflect the user's host machine TPN configuration. 

Test Configuration
--------------------

Additional configuration may be needed for the specific test the user plans to run. For information on test-specific configuration, view the information for specific tests: 

   * The :ref:`July 2020 CAPE Test Configuration <cape-config>`
   * The :ref:`Baroclinic Instability Test Configuration <bw-config>`

.. _RunExptC:

Run the Experiment
=====================

To start the experiment, run: 

.. code-block:: console
   
   cd $HSD/ufs-weather-model/tests-dev
   ./ufs_test.sh -a <ACCOUNT> -s -c -k -r -n "<CASE_NAME> <COMPILER>"

where:

* ``<ACCOUNT>``: Account/project number for batch jobs.
* ``<CASE_NAME>``: Name of the test case (e.g., ``2020_CAPE`` or ``baroclinic_wave``).
* ``<COMPILER>``: Compiler used for the tests (``intel`` or ``gnu``).

The script will loop until it runs both tasks or crashes. ``rococtostat`` can be used to track its progress; see the :ref:`Track Progress <TrackProgress>` section for details.

.. _TrackProgress:

Track Progress
----------------

To check on the job status, users on a system with a Slurm job scheduler may run (usually in a separate terminal window): 

.. code-block:: console

   squeue -u $USER

To view the experiment status, make sure that rocoto is loaded and run:

.. code-block:: console

   rocotostat -w rocoto_workflow.xml -d rocoto_workflow.db -v 10

It will print a status table:

.. code-block:: console

          CYCLE                     TASK   JOBID     STATE   EXIT STATUS  TRIES      DURATION
   ===========================================================================================
   197001010000  compile_atm_dyn32_intel       1   RUNNING             -      0           0.0
   197001010000          2020_CAPE_intel       -         -             -      -             -

If the job hangs or otherwise fails, stop the workflow in the active terminal using ``(Ctrl+C)``. To resubmit the experiment, remove the ``rocoto_workflow*`` files and lock directory before rerunning the ``ufs_test.sh`` script again:

.. code-block:: console

	rm -rf rocoto_workflow* lock

.. _CheckExptOutput:

Check Experiment Output
-------------------------

If the experiment completes successfully, the loop will exit with output similar to the following:

.. code-block:: console
   :emphasize-lines: 1, 2, 11

   Rocoto workflow has completed.
   + return 0
   + [[ true == true ]]
   + [[ '' != '' ]]
   ++ date '+%Y%m%d %T'
   + TEST_END_TIME='20241115 16:43:41'
   + export TEST_END_TIME
   + python -c 'import create_log; create_log.finish_log()'
   running: /usr/bin/singularity exec --env-file /scratch1/NCEPDEV/stmp4/User.Name/hsd-test/new-cont/ufs-weather-model/container-scripts/ufswm.env -B /scratch1:/scratch1 /scratch1/NCEPDEV/stmp4/User.Name/hsd-test/new-cont/ubuntu22.04-intel-wm-dev-hsd-test.img python tmp_arg_file.py
   Performing Cleanup...
   REGRESSION TEST RESULT: SUCCESS
   + echo 'ufs_test.sh finished'
   ufs_test.sh finished
   + cleanup
   ++ awk '{print $2}'
   + PID_LOCK=2947803
   + [[ 2947803 == \2\9\4\7\8\0\3 ]]
   + rm -rf /scratch1/NCEPDEV/stmp4/User.Name/hsd-test/new-cont/ufs-weather-model/tests-dev/lock
   + [[ false == true ]]
   + trap 0
   + exit

The experiment output can be found under the run directory (``${PTMP}/${USER}/FV3_RT/rt_${pid}``), which will contain two subdirectories, two log files, and two environment variable files (one for the compile task and one for the experiment task). For example: 

.. code-block:: console

   $ ls run_dir/
   baroclinic_wave_intel      compile_atm_dyn32_intel       compile_atm_dyn32_intel.log
   baroclinic_wave_intel.log  compile_atm_dyn32_intel.env   run_test_baroclinic_wave_intel.env
