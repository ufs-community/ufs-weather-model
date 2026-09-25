.. _container-rt-tests:

*************************************************************************
Container and Community Platform Workflows for the UFS Weather Model
*************************************************************************

This chapter describes two build-and-run workflows for the UFS Weather Model (:term:`WM`)
driven by the driver script ``tests/community.sh``:

* **Container workflow** (default): The model is compiled and run inside a
  Singularity/Apptainer software container that bundles the prerequisite compilers,
  MPI, and all required third-party libraries. The UFS WM source code is checked out
  from standard GitHub repositories and built inside the container environment,
  eliminating the need to install :term:`spack-stack` or host-specific modules.

* **Community platform workflow** (``-p`` flag): The container is bypassed entirely.
  The model is compiled and run natively using a software stack already installed on
  the host system and exposed through a user-provided Lmod modulefile. This mode
  requires natively installed software stack to already exist on a community platform.

.. attention::

   This chapter covers the container-based and community platform workflows driven by
   ``community.sh``. For the standard, Tier-1-oriented RT framework driven by ``rt.sh``,
   see :numref:`Section %s <UsingRegressionTest>`.

.. _container-rt-vs-rt:

==========================================================
Relationship to the Regression Test (RT) Workflow
==========================================================

``community.sh`` was intentionally designed to follow a logic and structure similar to
the Regression Test (:term:`RT`) framework (``rt.sh``/``rt.conf``). The goal of this
design choice is to facilitate testing and adoption of this new workflow. UFS WM
developers already familiar with running RTs would be able to use ``community.sh``
without learning an entirely new set of conventions. To that end, ``community.sh``
reuses several pieces of the existing RT infrastructure:

* Test cases are selected from the same test definition files under ``tests/tests``
  that ``rt.conf``/``rt.sh`` draw from, and CMake build options follow the same
  conventions as ``rt.conf``.
* ``rt_utils.sh`` is reused to configure some of the shell environment variables needed
  to run the tests.
* Job card templates in ``tests/fv3_conf/`` follow the same naming pattern used by the
  RT framework.

That said, ``community.sh`` is a separate driver script with its own configuration
file, ``community.conf``, and its own goal.

``community.sh`` is intended for portability testing on a contributor's own platform.
It gives the UFS-WMcommunity members and users outside the core UFS WM development team 
a simple starting point for building and running the model. This platform may be a Tier 1
system, another HPC center, a cloud instance, or a laptop/workstation, with or without
container software. Running the tests successfully with ``community.sh`` confirms that
the model has been **ported, built, and run successfully**, rather than reproducing the full
baseline-comparison capability maintained on Tier 1.

By contrast, ``rt.sh`` targets the officially supported NOAA Tier 1 RDHPC platforms
exclusively (such as Ursa, Gaea, Orion, Hercules, Derecho, and NOAA Cloud at the moment of
writing). It confirms that code changes preserve bit-for-bit baseline results, using
Rocoto/ECFlow workflow management.

Keeping ``community.sh`` and ``community.conf`` separate from ``rt.sh`` also means that
container support and community-platform complexity do not burden the production
``rt.sh`` CI/CD pipeline, and vice versa. ``community.sh`` runs sequentially — compiling each
configuration and then running its tests in turn, with no Rocoto or ECFlow workflow
manager involved — which keeps it simple for interactive debugging on whatever platform
the user has available.

Users are expected to modify the ``community.conf`` configuration to suit their computing
platform standards and job scheduler (if any). The locations of staged input data, the 
container image, the runtime directory, as well as host-system runtime modules need to be
adjusted to fit local data paths and user's environment. Users are encouraged to further
tailor the workflow to fit their own modeling needs beyond running predefined test cases.

.. _container-rt-prereqs:

=============
Prerequisites
=============

.. _container-rt-apptainer:

Singularity/Apptainer
-----------------------

Users running the workflow in container mode (the default, without the ``-p`` community
platform flag) must have **Singularity** or **Apptainer** software installed on their
compute platform. `Singularity/Apptainer <https://en.wikipedia.org/wiki/Apptainer#History>`_ container software is widely used in HPC environments to provide portable and reproducible software environments. It provides OS-level virtualization by packaging an application, its dependencies, and selected runtime environment components into a container image.  For MPI workflows, the host HPC system typically coordinates process launch and task initialization through its scheduler and runtime services, allowing the containerized application to integrate with compute nodes, interconnects, and parallel file systems.


For further information of container software, see:
*SingularityCE* `https://sylabs.io/singularity/ <https://sylabs.io/singularity/>`_ and
*Apptainer* `https://apptainer.org/ <https://apptainer.org/>`_


On many HPC systems, Singularity/Apptainer is available as a loadable module:


.. code-block:: console

   module load singularity
   # or
   module load apptainer

When not available system-wide, Apptainer can be installed on a Linux-based system by following the `Apptainer Installation Guide <https://apptainer.org/docs/admin/latest/installation.html>`__.

The following table lists the container software and the module load command on NOAA RDHPC Tier 1 platforms:

.. list-table:: Container software on NOAA RDHPC Tier 1 platforms
   :widths: 25 25 30
   :header-rows: 1

   * - Machine
     - Container command
     - Module to load
   * - Ursa
     - ``apptainer``
     - none required
   * - Gaea-C6
     - ``apptainer``
     - none required
   * - Hercules/Orion
     - ``singularity``
     - ``module load singularity``
   * -
     - ``apptainer`` (*)
     - ``module load spack-managed-x86-64_v3/v1.0 apptainer``
   * - Derecho
     - ``apptainer``
     - ``module load apptainer``
   * - NOAA Cloud (AWS/Azure)
     - ``singularity``
     - none required

(*) - The ``apptainer`` module on Hercules/Orion is a Spack-managed
install that loads a separate environment, which may not combine well
with other system modules. The ``apptainer`` enables certain container
build features
that are otherwise limited in ``singularity`` module by security
constraints. The ``singularity`` module could further be used for compile
and runtime environments.

.. note::

   Apptainer is fully compatible with Singularity, and commands shown with ``singularity`` may be replaced with ``apptainer`` as appropriate.
   When using Apptainer, prefer the ``APPTAINER_`` environment-variable prefix
   instead of the legacy ``SINGULARITY_`` prefix.

Further information on Singularity/Apptainer is available at:

- `Apptainer documentation <https://apptainer.org/docs/>`__
- `SingularityCE documentation <https://docs.sylabs.io/guides/latest/user-guide/>`__
- `NOAA RDHPCS container documentation <https://docs.rdhpcs.noaa.gov/software/containers>`__

.. _container-rt-image:

Container Image
---------------

The container RT workflow requires a Singularity/Apptainer image (``*.sif``). Both GNU-based and Intel-based images are supported.

.. note::

   The method for accessing pre-staged container images on NOAA RDHPC Tier 1 platforms and
   for building containers from Docker Hub on other community platforms is identical to the
   procedure described in the `UFS Short-Range Weather (SRW) App Container Quick Start Guide
   <https://ufs-srweather-app.readthedocs.io/en/latest/BuildingRunningTesting/ContainerQuickstart.html>`__.
   Users already familiar with SRW App containers may refer to that guide for additional
   context, troubleshooting tips, and platform-specific notes.

.. _container-rt-image-tier1:

On NOAA RDHPC Tier 1 Platforms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pre-built images for both GNU and Intel toolchains are staged at system-specific shared directories:

.. list-table:: Pre-built container image locations on Tier 1 platforms
   :widths: 20 50
   :header-rows: 1

   * - Machine
     - Directory
   * - Ursa
     - ``/scratch3/NCEPDEV/nems/role.epic/containers``
   * - Gaea-C6
     - ``/gpfs/f6/bil-fire8/world-shared/containers``
   * - Hercules / Orion
     - ``/work/noaa/epic/role-epic/contrib/containers``
   * - Derecho
     - ``/glade/work/epicufsrt/contrib/containers``
   * - NOAA Cloud
     - ``/contrib/EPIC/containers``

The container image file names are:

.. list-table:: Container image file names
   :widths: 15 25 40
   :header-rows: 1

   * - Toolchain
     - Platform
     - Image file name
   * - GNU (GCC 13.3.1 / OpenMPI 4.1.6)
     - Ursa, Gaea-C6, Hercules, Orion, NOAA Cloud
     - ``rocky9-gcc13-ss192-ompi416.sif``
   * - GNU (GCC 13.3.1 / OpenMPI 5.0.7)
     - Derecho
     - ``rocky9-gcc13-ss192-ompi507.sif``
   * - Intel (oneAPI 2024.2 / Intel MPI 2021.13)
     - Ursa, Gaea-C6, Hercules, Orion, NOAA Cloud
     - ``rocky9-oneapi2024.2-ss192.sif``

.. note::

   Derecho uses a GNU container image built with **OpenMPI 5.0.7** (``ompi507``) rather
   than 4.1.6 (``ompi416``) used on other platforms, due to MPI compatibility requirements
   on that system. The Intel container image has **not** been tested on Derecho to date;
   only the GNU container is currently supported there.

For example, on Hercules or Orion the Intel image is at:

.. code-block:: console

   /work/noaa/epic/role-epic/contrib/containers/rocky9-oneapi2024.2-ss192.sif

.. _container-rt-image-build:

Building a Container Image on Other Systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On systems where a pre-staged image is not available, build a Singularity/Apptainer image from Docker Hub.

.. note::

   Building a container image or sandbox requires temporary disk space in ``/tmp`` or ``~/.singularity/cache``.
   On systems with limited home-directory space, redirect the cache and temp directories before building:

   .. code-block:: console

      export SINGULARITY_CACHEDIR=/path/to/large/filesystem/cache
      export SINGULARITY_TMPDIR=/path/to/large/filesystem/tmp

   When using Apptainer, use ``APPTAINER_CACHEDIR`` and ``APPTAINER_TMPDIR`` instead.

**Option 1: Build a GNU-based container from Docker Hub**

On most platforms where the container does not already exist, run:

.. code-block:: console

   singularity build rocky9-gcc13-ss192-ompi416.sif \
       docker://noaaepic/rocky9-gcc13.3.1-spack-stack:v1.9.2-ufs-env-ompi416

On **Derecho**, use the OpenMPI 5.0.7 variant instead:

.. code-block:: console

   singularity build rocky9-gcc13-ss192-ompi507.sif \
       docker://noaaepic/rocky9-gcc13.3.1-spack-stack:v1.9.2-ufs-env-ompi507

**Option 2: Build an Intel-capable container from Docker Hub**

The Intel oneAPI software cannot be distributed inside Docker Hub images due to Intel's End User License Agreement (EULA). An Intel-capable software-stack image is available on Docker Hub, but the Intel oneAPI compilers and MPI must be reinstalled locally into a writable sandbox. The steps below produce a fully functional Intel image.

.. note::

   Site-specific SingularityCE installations may restrict sandbox builds more than Apptainer.
   If you encounter errors with SingularityCE, use Apptainer for the build steps and
   SingularityCE for runtime afterwards. On Hercules/Orion, a build-capable Apptainer
   is available via:

   .. code-block:: console

      module load spack-managed-x86-64_v3/v1.0 apptainer/1.3.3

#. Create a writable sandbox from the Docker Hub Intel-capable image. Bind all host
   top-level filesystems that contain your working directories (replace ``</top_dir>``
   and optional ``</bind_add>`` with paths appropriate for your system — see
   :numref:`Section %s <container-rt-binddirs>`):

   .. code-block:: console

      singularity build --sandbox --fix-perms  rocky9-oneapi2024.2-ss192 \
          docker://noaaepic/rocky9-oneapi2024.2-spack-stack:v1.9.2-ufs-wm-env

#. Copy the helper scripts out of the sandbox:

   .. code-block:: console

      singularity exec rocky9-oneapi2024.2-ss192 cp /opt/intel-sandbox.sh .
      singularity exec rocky9-oneapi2024.2-ss192 cp /opt/compilers_cp.sh .

   The scripts ``intel-sandbox.sh`` and ``compilers_cp.sh`` retrieve and reinstall the
   Intel compiler and MPI components.

#. Create the Intel oneAPI source sandbox:

   .. code-block:: console

      ./intel-sandbox.sh

   This produces an additional ``intel-sandbox`` directory containing the Intel oneAPI
   compilers and MPI.

#. Copy the Intel compilers and MPI into the software-stack sandbox. Provide only the
   sandbox names (not full paths):

   .. code-block:: console

      ./compilers_cp.sh intel-sandbox rocky9-oneapi2024.2-ss192

   After this step, the software-stack sandbox contains the full Intel toolchain.
   The ``intel-sandbox`` directory can then be removed.

#. Convert the sandbox into a compressed SIF image:

   .. code-block:: console

      singularity build --fix-perms rocky9-oneapi2024.2-ss192.sif rocky9-oneapi2024.2-ss192

.. _container-rt-binddirs:

Bind Directories for Tier 1 Platforms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following table lists the typical bind directories for NOAA RDHPC Tier 1 platforms.
These paths should be provided as a comma-separated list in the ``BIND_DIRS`` field of
``community.conf`` header line 1 (see :numref:`Section %s <container-rt-conf>`):

.. list-table:: Typical bind directories on NOAA RDHPC Tier 1 platforms
   :widths: 25 35 40
   :header-rows: 1

   * - Machine
     - Main bind directory
     - Additional bind directory
   * - Derecho
     - ``/glade``
     - none
   * - Ursa
     - ``/scratch3``
     - ``/scratch4``
   * - Gaea-C6
     - ``/gpfs``
     - ``/ncrc/home2``
   * - Hercules / Orion
     - ``/work``
     - ``/work2``, ``/local``
   * - NOAA Cloud (AWS/Azure)
     - ``/contrib``
     - ``/lustre``

.. _container-rt-data:

Input and Baseline Data
-----------------------

The container and community workflow uses the same input datasets as the standard RT framework.
On Level 1 and Level 2 systems these are pre-staged; see :numref:`Section %s <DataLocations>` for
the ``DISKNM`` and ``INPUTDATA_ROOT`` paths for each platform. These paths are set in
``community.conf`` (see :numref:`Section %s <container-rt-conf>`).

For Level 3–4 systems, input data is publicly available in the `UFS WM Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`__.

.. _container-rt-setup:

==================
Repository Setup
==================

Clone the UFS Weather Model repository and its submodules:

.. code-block:: console

   git clone --recursive https://github.com/ufs-community/ufs-weather-model.git
   cd ufs-weather-model

All further steps in this section assume the working directory is the root of the cloned repository.

.. _container-rt-modulefiles:

Required Modulefiles
--------------------

The modulefiles required depend on the workflow selected. The container option (default)
requires a user-adapted modulefile to load any host system modules during the runtime.
The community platform option (``-p`` flag) requires a modulefile to load all the required
software stack libraries on the platform. All modulefiles are placed in the
``modulefiles/`` directory at the root of the repository.

.. _container-rt-modulefiles-container:

Container Option (Default)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The container workflow uses two modulefiles. Only ``ufs_container.runtime.lua``
**must be adapted** by the user. The ``ufs_container.<compiler>.lua`` build module
depends on the software stack inside the container image and does not
require user changes.

.. _container-rt-runtime-mod:

``ufs_container.runtime.lua`` — Host-Side Runtime Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This modulefile is loaded **on the host** by ``community.sh`` and by the compile and run
job cards. **Users must create and adapt this file** for their platform. Its content
depends on the MPI launch method:

* **Slurm** (``srun``): ``srun`` coordinates MPI rank launch across compute nodes via the
  host Process Management Interface. The GNU-based image with OpenMPI 4.1.6 supports PMI2
  (``--mpi=pmi2``); the image with OpenMPI 5.0.7 supports PMIx (``--mpi=pmix``). Run
  ``srun --mpi=list`` to confirm availability. In this case no host MPI libraries are
  required and the modulefile only needs to load the Singularity/Apptainer module.

* **PBS** (``mpirun``/``mpiexec``): the host MPI launcher requires ABI-compatible MPI
  libraries on the host. Load the Singularity/Apptainer module together with compiler and
  MPI modules that match the container’s toolchain.

.. warning::

   Mismatched MPI implementations, incompatible MPI versions, or incompatible PMI/PMIx
   support may lead to runtime failures, hangs, or incorrect behavior.

A minimal example for Hercules or Orion, where only the Singularity module needs to be loaded:

.. code-block:: lua

   -- modulefiles/ufs_container.runtime.lua
   -- Host-side runtime environment for container-based UFS-WM RTs (Hercules/Orion)
   whatis("Host runtime module: loads singularity for container RT jobs")

   load("singularity")

On Derecho, which uses a PBS Pro job scheduler with ``mpirun``/``mpiexec`` as the MPI
launcher, GNU and OpenMPI host modules that are ABI-compatible with the container must be
loaded alongside the Apptainer module:

.. code-block:: lua

   -- modulefiles/ufs_container.runtime.lua
   -- Host-side runtime environment for container-based UFS-WM RTs (Derecho)
   whatis("Host runtime module: loads apptainer, GNU compilers, and host OpenMPI for container RT jobs")

   load("apptainer")
   load("gcc/14.3.0")
   load("openmpi/5.0.9")

Adapt the module names as needed. On systems where Singularity/Apptainer is already in
``PATH`` (e.g., Gaea-C6, NOAA Cloud), this file may be left empty or only load
supplementary host libraries needed by the MPI launcher.

.. _container-rt-build-mod:

``ufs_container.<compiler>.lua`` — Inside-Container Build Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This modulefile is loaded **inside the container** during both the compile and run stages.
It sets up the compiler toolchain, MPI library, and all required software libraries
available within the container image.

.. note::

   This file is pre-configured to match the software stack inside the container image and
   **does not normally require user changes**. A working example is provided in
   ``modulefiles/ufs_container.<compiler>.lua`` within the repository. Users should only
   modify it if they need to customize the inside-container environment (e.g., to override
   specific library versions).


A minimal example for an Intel-based container image:

.. code-block:: lua

   -- modulefiles/ufs_container.intel.lua
   -- Inside-container build/run environment for Intel-based UFS-WM container
   whatis("Inside-container software environment: Intel compilers, MPI, and spack-stack libraries")

   -- The container image ships with a self-contained module system.
   -- Adjust paths to match the software stack inside the image.
   prepend_path("MODULEPATH", "/opt/spack-stack/envs/ufs-wm/install/modulefiles/Core")

   load("stack-intel/2024.2.1")
   load("stack-intel-oneapi-mpi/2021.13")
   load("ufs-weather-model-env")

The exact module names depend on the container image being used. To discover available
modules, open an interactive shell inside the container (example for the Intel-based
container on Hercules/Orion):

.. code-block:: console

   singularity shell -B /work -B /work2 -B /local <container-image>
   # inside the container:
   source /opt/spack-stack/spack-stack-1.9.2/.bashenv   # or equivalent init script
   module avail
   module load stack-oneapi/2024.2.1
   module load stack-intel-oneapi-mpi/2021.13
   module avail

where ``<container-image>`` is the path to the container image file (``*.sif``) on the host.

.. _container-rt-modulefiles-community:

Community Platform Option (``-p`` flag)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When running with the ``-p`` flag, the model is built and run natively on the host — no
container is involved. Natively installed software stack is expected to be present. 

.. _container-rt-community-mod:

``ufs_<MACHINE_ID>.<compiler>.lua`` — Platform Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Users must create and adapt this file** for their platform. It is loaded during both
the compile and run stages and sets up the compiler toolchain, MPI library, and all
required software libraries available on the host system.

The file must be named ``modulefiles/ufs_<MACHINE_ID>.<compiler>.lua``, where
``<MACHINE_ID>`` matches the value in header line 1 of the conf file and ``<compiler>``
is ``intel`` or ``gnu``.

An example for a GNU-based native stack:

.. code-block:: lua

   -- modulefiles/ufs_myplatform.gnu.lua
   whatis("Native build/run environment for UFS-WM (GNU)")

   prepend_path("MODULEPATH", "/path/to/spack-stack/envs/ufs-wm/install/modulefiles/Core")
   load("stack-gcc/13.3.0")
   load("stack-openmpi/4.1.6")
   load("ufs-weather-model-env")

.. _container-rt-conf:

====================================
Configuring ``community.conf``
====================================

The file ``tests/community.conf`` controls what gets compiled and tested.
It begins with four mandatory header lines followed by one or more compile configuration blocks,
each with a list of test cases.

The file uses ``|`` as a field separator and ``#`` for comments. Blank lines between configuration blocks are ignored.

**Container mode example:**

.. code-block:: text

   # tests/community.conf — container mode example

   # Header line 1: MACHINE_ID | RT_COMPILER | CONTAINER_IMG | BIND_DIRS
   container | intel | /work/noaa/epic/role-epic/contrib/containers/rocky9-oneapi2024.2-ss192.sif | /work,/work2,/local

   # Header line 2: TPN | SCHEDULER | ACCNR | PARTITION | QUEUE | MPI_LAUNCH
   80 | slurm | epic | hercules | batch |

   # Header line 3: RUNDIR_ROOT
   /work2/noaa/epic/nperlin/hercules/UFS-WM/ufs-weather-model/tests/run_container

   # Header line 4: INPUTDATA_ROOT | INPUTDATA_ROOT_WW3 | INPUTDATA_LM4 | INPUTDATA_GFSv17opn
   /work2/noaa/epic/hercules/UFS-WM_RT/NEMSfv3gfs/input-data-20251015 | | | /work2/noaa/epic/hercules/UFS-WM_RT/NEMSfv3gfs/GFSv17opn_20251014

   # Compile configuration block: compile_id | MAKE_OPT
   atm | -DAPP=ATM -DCCPP_SUITES=FV3_GFS_v16,FV3_GFS_v17_p8
   control_c48
   control_p8

**Community platform mode example** (``-p`` flag, no container image needed):

.. code-block:: text

   # tests/community.conf — community platform mode example

   # Header line 1: MACHINE_ID | RT_COMPILER | CONTAINER_IMG | BIND_DIRS
   myplatform | gnu |  |

   # Header line 2: TPN | SCHEDULER | ACCNR | PARTITION | QUEUE | MPI_LAUNCH
   96 |  |  |  |  | mpirun

   # Header line 3: RUNDIR_ROOT
   /scratch/nperlin/ufs-weather-model/tests/run_myplatform

   # Header line 4: INPUTDATA_ROOT | INPUTDATA_ROOT_WW3 | INPUTDATA_LM4 | INPUTDATA_GFSv17opn
   /data/UFS-WM_INPUT/input-data-20251015 | | | /data/UFS-WM_INPUT/GFSv17opn_20251014

   atm | -DAPP=ATM -DCCPP_SUITES=FV3_GFS_v16,FV3_GFS_v17_p8
   control_c48

Header Line Fields
------------------

**Header line 1:**

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Field
     - Description
   * - ``MACHINE_ID``
     - Platform identifier. Use ``container`` for the container workflow, or a custom
       name matching the modulefile name for community platform (``-p``) runs.
       See :numref:`Section %s <container-rt-community-mod>`.
   * - ``RT_COMPILER``
     - Compiler toolchain: ``intel`` or ``gnu``.
   * - ``CONTAINER_IMG``
     - Absolute path to the Singularity/Apptainer image file (``*.sif``) on the host.
       Leave blank when using the ``-p`` flag.
   * - ``BIND_DIRS``
     - Comma-separated list of host directories to bind/mount to the container.
       Include all filesystems containing the source tree, input data, and run directory.
       See :numref:`Section %s <container-rt-binddirs>` for typical values on Tier 1 platforms.
       Leave blank when using the ``-p`` flag.

**Header line 2:**

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Field
     - Description
   * - ``TPN``
     - MPI tasks per node (default: 40).
   * - ``SCHEDULER``
     - Job scheduler: ``slurm``, ``pbs``, or leave blank for interactive/no-scheduler runs.
   * - ``ACCNR``
     - Scheduler account or project name (leave blank if not required).
   * - ``PARTITION``
     - Slurm partition name (leave blank for PBS or interactive runs).
   * - ``QUEUE``
     - Slurm QOS / PBS queue name (leave blank for interactive runs).
   * - ``MPI_LAUNCH``
     - MPI launch command used when no scheduler is set: ``mpirun`` or ``mpiexec``.
       Defaults to ``mpirun`` if omitted.

**Header line 3:**

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Field
     - Description
   * - ``RUNDIR_ROOT``
     - Top-level directory where compile and test run directories will be created.
       This should be a user-writable path (preferably on a scratch or work filesystem).
       If ``RUNDIR_ROOT`` differs from ``${PATHRT}/run_dir``, the driver automatically
       creates a convenience symlink ``tests/run_dir`` pointing to ``RUNDIR_ROOT``.

**Header line 4:**

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Field
     - Description
   * - ``INPUTDATA_ROOT``
     - Input data directory (required).
   * - ``INPUTDATA_ROOT_WW3``
     - Input data directory for WaveWatch III data (optional; leave blank if not needed).
   * - ``INPUTDATA_LM4``
     - Input data directory for LM4 land model data (optional; leave blank if not needed).
   * - ``INPUTDATA_GFSv17opn``
     - Input data directory for GFS v17 operational data (optional; leave blank if not needed).

Compile Configuration Blocks
-----------------------------

After the four header lines, each compile configuration block consists of:

1. A **compile line** containing the configuration name and CMake options, separated by ``|``:

   .. code-block:: text

      <compile_id> | <MAKE_OPT>

2. One or more **test case lines**, each naming a single test (no ``|`` character):

   .. code-block:: text

      <test_case_name>

A blank line between blocks closes the current configuration. Configuration names must be unique. CMake options follow the same conventions as the standard ``rt.conf`` file.

.. note::

   The workflow is **sequential**: the driver compiles configuration 1, runs all its tests in order,
   then moves to configuration 2, and so on. A startup summary listing all configurations and their
   test cases is printed before any work begins.

.. _container-rt-run:

====================
Running the Tests
====================

All tests are launched by running ``community.sh`` from the ``tests/`` directory:

.. code-block:: console

   cd ${WM_HOME}/tests
   ./community.sh [options] <conf_file>

The configuration file is a required positional argument. Provide the path to the
``community.conf``-style file as the last argument.

Command-Line Options
--------------------

.. list-table::
   :widths: 10 60
   :header-rows: 1

   * - Option
     - Description
   * - ``-p``
     - Community platform mode: build and run natively without a container.
       Requires a ``ufs_<MACHINE_ID>.<compiler>.lua`` modulefile in ``modulefiles/``;
       see :numref:`Section %s <container-rt-community-mod>`.
   * - ``-d``
     - Delete each test run directory after the test completes.
   * - ``-n <name>``
     - Run only the single test named ``<name>`` (the compile step that owns it still runs).
       ``<name>`` must match a test case name listed in ``<conf_file>``.
   * - ``-o``
     - Compile only; skip all test cases.
   * - ``-v``
     - Verbose output: enables shell tracing (``set -x``) in the driver and all sub-scripts,
       prints full configuration detail on startup, and prompts before starting work.
   * - ``-h``
     - Print help and exit.

Job Script Templates
--------------------

When ``SCHEDULER`` is set to ``slurm`` or ``pbs``, the driver uses job script templates
from ``tests/fv3_conf/``. The template name is based on the scheduler type and the
``MACHINE_ID`` value from the conf file:

- ``fv3_slurm.IN_<MACHINE_ID>`` — Slurm run job card
- ``fv3_qsub.IN_<MACHINE_ID>`` — PBS run job card
- ``compile_slurm.IN_<MACHINE_ID>`` — Slurm compile job card
- ``compile_qsub.IN_<MACHINE_ID>`` — PBS compile job card

For the container workflow (``MACHINE_ID=container``), templates are provided in the
repository. For community platform runs with a scheduler, users must create
platform-specific templates following the pattern of the container or Tier 1 templates
in ``tests/fv3_conf/``. When ``SCHEDULER`` is blank, no job template is used and
compile and run steps execute directly on the current host.

These templates contain scheduler directives and environment setup that may require
platform-specific adjustments. Scheduler adjutments may include account/project names,
partition or queue names, wall-clock limits, node counts, and any
platform-specific environment variables required before launching the model.

Running with a Job Scheduler (Slurm or PBS)
--------------------------------------------

When ``SCHEDULER`` is set to ``slurm`` or ``pbs`` in ``community.conf``, set ``ACCNR``, ``PARTITION``, and ``QUEUE`` appropriately and run the driver from a login node:

.. code-block:: console

   ./community.sh [-p] community.conf

The driver submits each compile and test job to the scheduler and blocks until the job finishes before submitting the next one. Progress is reported on the terminal; full output is captured in ``${RUNDIR_ROOT}/logs/``.

Running Interactively (No Scheduler)
--------------------------------------

When ``SCHEDULER`` is blank in ``community.conf``, jobs run directly on the current host —
suitable for an allocated compute node or single-workstation development.

If not allowed using a login node for runtime tests, request an interactive compute node
allocation before running
the driver. On **Slurm** systems the command may look similar to:

.. code-block:: console

   salloc -N 1 -n <cores> -A <account> -t <time> -q <qos> --partition=<partition>

On **PBS** systems the command may look similar to:

.. code-block:: console

   qsub -I -l walltime=<time> -A <account> -q <queue> -l select=1:ncpus=<cores>:mpiprocs=<cores>

After the allocation is granted (and connecting via ``ssh`` to the compute node if required), run the driver.
A command set as ``MPI_LAUNCH`` in header line 2 of the ``community.conf`` to start MPI tasks will be used.

.. code-block:: console

   cd tests
   ./community.sh [-p] community.conf


.. note::

    For a **container workflow**, the driver starts the software container first, and then runs MPI tasks
    entirely inside it, which is the correct approach for single-node interactive runs.

   The ``--mpi=pmi2`` flag is a Slurm ``srun``-specific option and should **not** be used
   with ``mpirun`` or ``mpiexec``. When no scheduler is set, ``srun`` is not used and
   no container is launched for community platform (``-p``) runs.

.. _container-rt-output:

===================
Run Directory
===================

After the driver starts, it creates the following structure under ``RUNDIR_ROOT``.

**Container mode** (``MACHINE_ID=container``):

.. code-block:: text

   ${RUNDIR_ROOT}/
   ├── logs/                         # per-job log files and timestamps
   ├── compile_<compile_id>/         # compile working directory
   │   ├── job_card                  # generated compile job script (scheduler) or absent
   │   ├── container_compile.sh      # script executed inside the container
   │   ├── modulefiles/              # modulefile staged for the build
   │   └── out / err                 # job stdout and stderr files
   └── <test_id>_<compiler>/         # test working directory
       ├── job_card                  # generated test job script (scheduler) or absent
       ├── fv3_container_run.sh      # wrapper executed inside the container
       ├── modulefiles/              # modulefile staged for the run
       └── out / err                 # job stdout and stderr files

**Community platform mode** (``-p`` flag):

.. code-block:: text

   ${RUNDIR_ROOT}/
   ├── logs/                         # per-job log files and timestamps
   ├── compile_<compile_id>/         # compile working directory
   │   ├── job_card                  # generated compile job script (if scheduler is used)
   │   ├── modulefiles/              # modulefile staged for the build
   │   └── out / err                 # job stdout and stderr files
   └── <test_id>_<compiler>/         # test working directory
       ├── job_card                  # generated test job script (if scheduler is used)
       ├── fv3_run.sh                # native run wrapper (when no scheduler)
       ├── modulefiles/              # modulefile staged for the run
       └── out / err                 # job stdout and stderr files

A symlink at ``tests/run_dir`` is created pointing to ``RUNDIR_ROOT``, making it easy to navigate to the run directory without knowing the full path. This symlink is not created if the user has already set ``RUNDIR_ROOT`` to ``${PATHRT}/run_dir``.

A PASS/FAIL summary is printed to the terminal when all tests have finished. The driver exits with status 1 if any compile or test failed, and 0 if all succeeded.
