.. _container-rt-tests:

******************************************************************************
Container and Community Platform Options to Test and Run the UFS Weather Model
******************************************************************************

This chapter describes two ways to build and run the UFS Weather Model (:term:`WM`)
outside of the standard Tier 1 software stacks. Both use the regression test script
``tests/rt.sh`` with the ``-P <platform.def>`` option, where ``<platform.def>`` is a
*platform definition file* that describes the platform (see :numref:`Section %s <container-rt-conf>`):

* **Container option**: The model is compiled and run inside a
  Singularity/Apptainer software container that bundles the prerequisite compilers,
  MPI, and all required third-party libraries. The UFS WM source code is checked out
  from standard GitHub repositories and built inside the container environment,
  eliminating the need to install :term:`spack-stack` or host-specific modules.
  This option is selected when the platform definition file gives a container image.

* **Community platform option** (native software stack): No container is used.
  The model is compiled and run natively using a software stack already installed on
  the host system and exposed through a user-provided Lmod modulefile. This option
  is selected when the container image field in the platform definition file is left blank.

The container option has been tested on the community platform Stampede3 (TACC, University of
Texas) with both GNU-based and Intel-based software containers. An example platform definition
file for Stampede3, ``tests/stampede.def``, is included (see :numref:`Section %s <container-rt-conf>`).
The container option has also been tested on NOAA RDHPC Tier 1 platforms with the same method
(``-P <platform.def>``), by adapting the platform definition file for each platform and container
(GNU- or Intel-based).

.. attention::

   This chapter covers ``rt.sh`` runs with the ``-P <platform.def>`` option. For standard Tier 1
   regression testing, see :numref:`Section %s <UsingRegressionTest>`. A summary of
   the ``-P`` option is also given in :numref:`Section %s <rt-container>`.

.. _container-rt-vs-rt:

==========================================================
Relationship to the Regression Test (RT) Framework
==========================================================

Container and community platform runs use the same Regression Test (:term:`RT`) framework
as Tier 1 platforms:

* Tests are selected from ``tests/rt.conf`` (or from another file given with ``-l``), and the
  test definitions are the same files under ``tests/tests``. CMake build options follow the
  same conventions as ``rt.conf``.
* A ``COMPILE`` or ``RUN`` line in ``rt.conf`` runs on the platform only when its
  **Machines** column contains ``+<PLATFORM_NAME>``, where ``<PLATFORM_NAME>`` is the name
  given in line 1/ field 1 of the platform definition file. ``-<PLATFORM_NAME>`` excludes the line
  even if ``+<PLATFORM_NAME>`` is also present. The tests currently enabled for the container
  option are tagged ``+container`` in ``rt.conf``.
* The same ``rt.sh`` options are used, for example ``-a``, ``-l``, ``-n``, ``-o``, ``-c``,
  ``-m``, ``-r``, and ``-e`` (see :numref:`Section %s <container-rt-run>`).

The goal of a ``-P`` run differs from a Tier 1 run. By default, a ``-P`` run is a
**portability check**: it confirms that the model has been **ported, built, and run
successfully** on the platform, with the same input data, model configuration(s), and
test case definitions as the Tier 1 platforms. The ``-P`` run is not intended to be a
The platform may be a Tier 1 system, another community/HPC center system, a cloud instance,
or a laptop/workstation, with or without container software. Baselines can still be created
(``-c``) and compared against (``-m``) on the same platform. These baselines are kept
under the platform user's runtime directory, and are separate from the Tier 1 baselines.

By default, a ``-P`` run is **sequential**: ``rt.sh`` compiles each configuration and then
runs its tests in turn, which keeps it simple for interactive debugging. The Rocoto (``-r``)
or ecFlow (``-e``) workflow managers can be used instead when the platform definition file
provides the optional line 5 (see :numref:`Section %s <container-rt-conf>`).

Users are expected to adapt the platform definition file to their computing platform and job
scheduler (if any): the container image, the run directory, the location of the staged input
data, and the host-system runtime modules all need to match the local environment. Users are
encouraged to further tailor these runs to fit their own modeling needs beyond running
predefined test cases.

.. _container-rt-prereqs:

=============
Prerequisites
=============

.. _container-rt-apptainer:

Singularity/Apptainer
-----------------------

Users running with the container option (a platform definition file with a container image)
must have **Singularity** or **Apptainer** software installed on their
compute platform. `Singularity/Apptainer <https://en.wikipedia.org/wiki/Apptainer#History>`_
container software is widely used in HPC environments to provide portable and reproducible
software environments. Multi-node MPI runs require the MPI library inside the container to be
binary (ABI) compatible with the host MPI (for ``mpirun``/``mpiexec``) or with the host's PMI/PMIx
(for ``srun``), so that the MPI ranks can communicate across nodes; a mismatch can cause failures
or hangs. See `Hybrid MPI Model (Host + Container MPI)
<https://docs.rdhpcs.noaa.gov/software/containers/index.html#hybrid-mpi-model-host-container-mpi>`__
on the `Containers <https://docs.rdhpcs.noaa.gov/software/containers/index.html>`__
page of the NOAA RDHPCS documentation.


For further information of container software, see:
*SingularityCE* `https://sylabs.io/singularity/ <https://sylabs.io/singularity/>`_ and
*Apptainer* `https://apptainer.org/ <https://apptainer.org/>`_


On many HPC systems, Singularity/Apptainer is available as a loadable module:


.. code-block:: console

   module load singularity
   # or
   module load apptainer

When not available system-wide, Apptainer can be installed on a Linux-based system by following the `Apptainer Installation Guide <https://apptainer.org/docs/admin/latest/installation.html>`__.

The following table lists the container software and the module load command on several platforms:

.. list-table:: Examples of container software used on various platforms
   :widths: 25 25 30
   :header-rows: 1

   * - Machine
     - Container command
     - Module to load
   * - Stampede3
     - ``apptainer``
     - ``module load tacc-apptainer``
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

The container software module is loaded on the host system by the
``modulefiles/ufs_container.runtime.lua`` module before the container is started (see
:numref:`Section %s <container-rt-runtime-mod>`). Examples of the module loading lines are
given in that file: uncomment the lines needed for the platform, or add your own. For example,
to run a container on Stampede3, add:

.. code-block:: lua

   load("tacc-apptainer")

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

The container option requires a Singularity/Apptainer image (``*.sif``). Both GNU-based and
Intel-based images are supported. Users can either build their own image from the Docker Hub
images (see :numref:`Section %s <container-rt-image-build>`), or use an existing ``*.sif`` image
if one is already staged on the system (see :numref:`Section %s <container-rt-image-tier1>`).
Pre-built images on Stampede3 a staged in a user's directory providing an example. 
Container images are also staged on NOAA RDHPC Tier 1 platforms in common locations maintained by EPIC.

.. note::

   The method for accessing pre-staged container images on NOAA RDHPC Tier 1 platforms and
   for building containers from Docker Hub on other community platforms is identical to the
   procedure described in the `UFS Short-Range Weather (SRW) App Container Quick Start Guide
   <https://ufs-srweather-app.readthedocs.io/en/latest/BuildingRunningTesting/ContainerQuickstart.html>`__.
   Users already familiar with SRW App containers may refer to that guide for additional
   context, troubleshooting tips, and platform-specific notes.

.. _container-rt-image-build:

Building a Container Image on a Community Platform
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On a community platform, or on any system where a pre-staged image is not available, build a
Singularity/Apptainer image from the Docker Hub images maintained by EPIC.

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

.. _container-rt-image-tier1:

Pre-built Container Images
~~~~~~~~~~~~~~~~~~~~~~~~~~

If a ``*.sif`` image already exists on the system, it can be used directly. Pre-built images for
both GNU and Intel toolchains are staged in the following directories:

.. list-table:: Pre-built container image locations on various platforms
   :widths: 20 50
   :header-rows: 1

   * - Machine
     - Directory
   * - Stampede3 (*)
     - ``/work2/10000/nperlin/stampede3``
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

(*) - The Stampede3 images are staged in a user's directory, not in a common location
maintained by EPIC. The Tier 1 locations are common locations maintained by EPIC.

The container image file names are:

.. list-table:: Container image file names
   :widths: 15 25 40
   :header-rows: 1

   * - Toolchain
     - Platform
     - Image file name
   * - GNU (GCC 13.3.1 / OpenMPI 4.1.6)
     - Stampede3, Ursa, Gaea-C6, Hercules, Orion, NOAA Cloud
     - ``rocky9-gcc13-ss192-ompi416.sif``
   * - GNU (GCC 13.3.1 / OpenMPI 5.0.7)
     - Derecho
     - ``rocky9-gcc13-ss192-ompi507.sif``
   * - Intel (oneAPI 2024.2 / Intel MPI 2021.13)
     - Stampede3, Ursa, Gaea-C6, Hercules, Orion, NOAA Cloud
     - ``rocky9-oneapi2024.2-ss192.sif``

.. note::

   Derecho uses a GNU container image built with **OpenMPI 5.0.7** (``ompi507``) rather
   than 4.1.6 (``ompi416``) used on other platforms, due to MPI compatibility requirements
   on that system. The Intel container image has **not** been tested on Derecho to date;
   only the GNU container is currently supported there.

For example, on Hercules or Orion the Intel image is at:

.. code-block:: console

   /work/noaa/epic/role-epic/contrib/containers/rocky9-oneapi2024.2-ss192.sif

.. _container-rt-binddirs:

Bind Directories for Tier 1 Platforms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following table lists the typical bind directories for NOAA RDHPC Tier 1 platforms.
These paths should be provided as a comma-separated list in the ``BIND_DIRS`` field of
line 1 of the platform definition file (see :numref:`Section %s <container-rt-conf>`):

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

Input Data
----------

Container and community platform runs use the same input datasets as the standard RT framework.
On Level 1 and Level 2 systems these are pre-staged; see :numref:`Section %s <DataLocations>` for
the ``DISKNM`` and ``INPUTDATA_ROOT`` paths for each platform. These paths are set in
line 4 of the platform definition file (see :numref:`Section %s <container-rt-conf>`).

For Level 3–4 systems, input data is publicly available in the `UFS WM Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`__.
The current input data sets are ``input-data-20260617`` (``INPUTDATA_ROOT``), with WaveWatch III data in
``input-data-20260617/WW3_input_data_20260811`` (``INPUTDATA_ROOT_WW3``) and LM4 data in
``input-data-20260617/LM4_input_data`` (``INPUTDATA_LM4``).

The complete ``input-data-20260617`` data set is large. To run only the tests tagged ``+container`` in
``tests/rt.conf``, a subset of about 121 GiB is enough. The table below lists, for each ``+container``
test, the subdirectories of ``INPUTDATA_ROOT`` (``input-data-20260617`` in the data bucket) that the test
copies input data from. Only those subdirectories need to be downloaded. "(files only)" means that only
the files directly in that directory are used, not its subdirectories. The compilers in parentheses are
those with which the test is tagged ``+container``. The input data is the same for both compilers,
except for ``control_c48``.

.. list-table:: Input data subdirectories used by each ``+container`` test
   :widths: 18 30 52
   :header-rows: 1

   * - Configuration (compilers)
     - Test
     - Input data subdirectories
   * - ``s2swa_32bit`` (intel, gnu)
     - ``cpld_control_p8``
     - | ``CICE_FIX/100``
       | ``CICE_IC/100``
       | ``CPL_FIX/aC96o100``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT``
       | ``FV3_input_data/INPUT_L127_gfsv17``
       | ``FV3_input_data/INPUT_L127_mx100``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/100``
       | ``MOM6_IC`` (files only)
       | ``WW3_input_data_20260811`` (files only)
   * - ``s2sw_32bit_pdlib`` (intel, gnu)
     - ``cpld_control_gfsv17``
     - | ``CICE_FIX/100``
       | ``CICE_IC/100``
       | ``CPL_FIX/aC96o100``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT``
       | ``FV3_input_data/INPUT_L127_gfsv17``
       | ``FV3_input_data/INPUT_L127_mx100``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/100``
       | ``MOM6_IC`` (files only)
       | ``WW3_input_data_20260811`` (files only)
   * - ``s2s_32bit_sfs`` (intel)
     - ``cpld_control_sfs``
     - | ``CICE_FIX/025``
       | ``CPL_FIX/aC192o025``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C192mx025``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data192`` (files only)
       | ``FV3_input_data192/INPUT``
       | ``FV3_input_data192/INPUT_L127_gfsv17``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/025``
       | ``SFS/1994050100``
   * - ``s2s_32bit_sfs_debug`` (intel, gnu)
     - ``cpld_debug_sfs``
     - | ``CICE_FIX/025``
       | ``CPL_FIX/aC192o025``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C192mx025``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data192`` (files only)
       | ``FV3_input_data192/INPUT``
       | ``FV3_input_data192/INPUT_L127_gfsv17``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/025``
       | ``SFS/1994050100``
   * - ``s2swl`` (intel)
     - ``cpld_control_p8_lnd``
     - | ``CICE_FIX/100``
       | ``CICE_IC/100``
       | ``CPL_FIX/aC96o100``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT``
       | ``FV3_input_data/INPUT_L127_gfsv17``
       | ``FV3_input_data/INPUT_L127_mx100``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/100``
       | ``MOM6_IC`` (files only)
       | ``WW3_input_data_20260811`` (files only)
   * - ``s2s_aoflux`` (intel)
     - | ``cpld_control_noaero_p8_``
       | ``agrid``
     - | ``CICE_FIX/100``
       | ``CICE_IC/100``
       | ``CPL_FIX/aC96o100``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT``
       | ``FV3_input_data/INPUT_L127_gfsv17``
       | ``FV3_input_data/INPUT_L127_mx100``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/100``
       | ``MOM6_IC`` (files only)
   * - ``s2sw_pdlib`` (intel)
     - ``cpld_control_c48_5deg``
     - | ``CICE_FIX/500``
       | ``CICE_IC/C48mx500/2021032206``
       | ``CPL_FIX/aC48o500``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C48mx500``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data48`` (files only)
       | ``FV3_input_data48/INPUT``
       | ``FV3_input_data48/INPUT_L127_gfsv17``
       | ``FV3_input_data48/INPUT_L127_mx500/2021032206``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/500``
       | ``MOM6_IC/C48mx500/2021032206``
       | ``WW3_input_data_20260811`` (files only)
   * - ``s2sw_pdlib`` (intel)
     - ``cpld_warmstart_c48_5deg``
     - | ``CICE_FIX/500``
       | ``CICE_IC/C48mx500/2021032306``
       | ``CMEPS_IC/C48mx500/2021032306``
       | ``CPL_FIX/aC48o500``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C48mx500``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data48`` (files only)
       | ``FV3_input_data48/INPUT``
       | ``FV3_input_data48/INPUT_L127_gfsv17``
       | ``FV3_input_data48/INPUT_L127_mx500/2021032306``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/500``
       | ``MOM6_IC/C48mx500/2021032306``
       | ``WW3_IC/C48mx500/2021032306``
       | ``WW3_input_data_20260811`` (files only)
   * - ``s2sw_pdlib`` (intel)
     - ``cpld_control_c24_5deg``
     - | ``CICE_FIX/500``
       | ``CICE_IC/C24mx500/2021032206``
       | ``CPL_FIX/aC24o500``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C24mx500``
       | ``FV3_input_data24`` (files only)
       | ``FV3_input_data24/INPUT``
       | ``FV3_input_data24/INPUT_L41_mx500/2021032206``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/500``
       | ``MOM6_IC/C24mx500/2021032206``
       | ``WW3_input_data_20260811`` (files only)
   * - ``s2sw_pdlib`` (intel)
     - ``cpld_warmstart_c24_5deg``
     - | ``CICE_FIX/500``
       | ``CICE_IC/C24mx500/2021032306``
       | ``CMEPS_IC/C24mx500/2021032306``
       | ``CPL_FIX/aC24o500``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C24mx500``
       | ``FV3_input_data24`` (files only)
       | ``FV3_input_data24/INPUT``
       | ``FV3_input_data24/INPUT_L41_mx500/2021032306``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/500``
       | ``MOM6_IC/C24mx500/2021032306``
       | ``WW3_IC/C24mx500/2021032306``
       | ``WW3_input_data_20260811`` (files only)
   * - ``s2sw_pdlib`` (intel)
     - ``cpld_control_c24_9deg``
     - | ``CICE_FIX/900``
       | ``CICE_IC/C24mx900/2021032206``
       | ``CPL_FIX/aC24o900``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C24mx900``
       | ``FV3_input_data24`` (files only)
       | ``FV3_input_data24/INPUT``
       | ``FV3_input_data24/INPUT_L41_mx900/2021032206``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/900``
       | ``MOM6_IC/C24mx900/2021032206``
       | ``WW3_input_data_20260811`` (files only)
   * - ``s2sw_pdlib`` (intel)
     - ``cpld_control_c12_9deg``
     - | ``CICE_FIX/900``
       | ``CICE_IC/C12mx900/2021032206``
       | ``CPL_FIX/aC12o900``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C12mx900``
       | ``FV3_input_data12`` (files only)
       | ``FV3_input_data12/INPUT``
       | ``FV3_input_data12/INPUT_L41_mx900/2021032206``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/900``
       | ``MOM6_IC/C12mx900/2021032206``
       | ``WW3_input_data_20260811`` (files only)
   * - ``atm_dyn32`` (intel, gnu)
     - ``control_CubedSphereGrid``
     - | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127_gfsv17``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
   * - ``atm_dyn32`` (intel)
     - ``control_c48``
     - | ``FV3_fix``
       | ``FV3_fix_tiled/C48mx500``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data48`` (files only)
       | ``FV3_input_data48/INPUT_L127_gfsv17``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
   * - ``atm_dyn32`` (gnu)
     - ``control_c48``
     - | ``FV3_fix``
       | ``FV3_fix_tiled/C48mx500``
       | ``FV3_input_data48`` (files only)
       | ``FV3_input_data48/INPUT_L127_gfsv17``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
   * - ``atm`` (gnu)
     - ``control_c48``
     - | ``FV3_fix``
       | ``FV3_fix_tiled/C48mx500``
       | ``FV3_input_data48`` (files only)
       | ``FV3_input_data48/INPUT_L127_gfsv17``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
   * - ``atm_dyn32`` (intel, gnu)
     - ``control_c192``
     - | ``FV3_fix``
       | ``FV3_fix_tiled/C192mx050``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data192`` (files only)
       | ``FV3_input_data192/INPUT_L127_gfsv17``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
   * - ``atm_dyn32`` (intel, gnu)
     - ``control_p8``
     - | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127_gfsv17``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
   * - ``atm`` (gnu)
     - ``control_p8``
     - | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127_gfsv17``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
   * - ``atm_dyn32`` (intel, gnu)
     - ``control_p8.v2.sfc``
     - | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127_v2_sfc``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
   * - ``atm_dyn32_rad32`` (intel)
     - ``control_p8_rrtmgp_rad32``
     - | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127_gfsv17``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``FV3_input_data_RRTMGP``
   * - ``rrfs`` (intel, gnu)
     - ``rap_control``
     - | ``FV3_fix``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127_gfsv17``
   * - ``rrfs`` (intel, gnu)
     - ``hrrr_control``
     - | ``FV3_fix``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127``
       | ``FV3_input_data_gsd`` (files only)
       | ``lake_p8_water_fraction2020``
   * - ``rrfs`` (intel, gnu)
     - ``rrfs_v1beta``
     - | ``FV3_fix``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127_gfsv17``
   * - ``wam`` (intel)
     - ``control_wam``
     - | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``FV3_input_data_L149_wam`` (files only)
       | ``FV3_input_data_L149_wam/INPUT_C``
   * - ``wam_debug`` (intel, gnu)
     - ``control_wam_debug``
     - | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``FV3_input_data_L149_wam`` (files only)
       | ``FV3_input_data_L149_wam/INPUT_C``
   * - ``rrfs_dyn32_phy32`` (intel, gnu)
     - ``rap_control_dyn32_phy32``
     - | ``FV3_fix``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127_gfsv17``
   * - ``rrfs_dyn32_phy32`` (intel, gnu)
     - ``hrrr_control_dyn32_phy32``
     - | ``FV3_fix``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127``
       | ``FV3_input_data_gsd`` (files only)
       | ``lake_p8_water_fraction2020``
   * - ``rrfs_dyn32_phy32`` (intel, gnu)
     - | ``hrrr_control_2threads_``
       | ``dyn32_phy32``
     - | ``FV3_fix``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127``
       | ``FV3_input_data_gsd`` (files only)
       | ``lake_p8_water_fraction2020``
   * - ``rrfs_dyn32_phy32`` (intel, gnu)
     - ``conus13km_control``
     - | ``FV3_aeroclim``
       | ``FV3_fix``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data_conus13km/INPUT``
   * - ``rrfs_dyn64_phy32`` (intel)
     - ``rap_control_dyn64_phy32``
     - | ``FV3_fix``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127_gfsv17``
   * - ``hafsw`` (intel)
     - ``hafs_regional_atm``
     - | ``FV3_fix``
       | ``FV3_hafs_input_data`` (files only)
       | ``FV3_hafs_input_data/INPUT_hafs_regional_atm``
   * - ``hafsw`` (intel)
     - ``hafs_regional_atm_wav``
     - | ``FV3_fix``
       | ``FV3_hafs_input_data`` (files only)
       | ``FV3_hafs_input_data/INPUT_hafs_regional_atm``
       | ``WW3_input_data_20260811`` (files only)
   * - ``hafsw`` (intel)
     - ``hafs_global_1nest_atm``
     - | ``FV3_fix``
       | ``FV3_hafs_input_data`` (files only)
       | ``FV3_hafs_input_data/INPUT_hafs_global_1nest_atm``
   * - ``hafsw`` (intel)
     - | ``hafs_global_multiple_``
       | ``4nests_atm``
     - | ``FV3_fix``
       | ``FV3_hafs_input_data`` (files only)
       | ``FV3_hafs_input_data/INPUT_hafs_global_multiple_4nests_atm``
   * - ``hafs_mom6w`` (intel)
     - | ``hafs_regional_storm_``
       | ``following_1nest_atm_ocn_``
       | ``wav_mom6``
     - | ``FV3_fix``
       | ``FV3_hafs_input_data`` (files only)
       | ``FV3_hafs_input_data/CDEPS_input_data``
       | ``FV3_hafs_input_data/INPUT_hafs_regional_storm_following_1nest_atm``
       | ``FV3_hafs_input_data/MOM6_regional_input_data``
       | ``FV3_hafs_input_data/WW3_hafs_regional_input_data``
       | ``WW3_input_data_20260811`` (files only)
   * - ``hafs_all`` (intel, gnu)
     - ``hafs_regional_docn``
     - | ``DOCN_MOM6_input_data``
       | ``FV3_fix``
       | ``FV3_hafs_input_data`` (files only)
       | ``FV3_hafs_input_data/INPUT_hafs_regional_atm``
   * - ``hafs_all`` (intel, gnu)
     - ``hafs_regional_docn_oisst``
     - | ``DOCN_OISST_input_data``
       | ``FV3_fix``
       | ``FV3_hafs_input_data`` (files only)
       | ``FV3_hafs_input_data/INPUT_hafs_regional_atm``
   * - ``datm_cdeps`` (intel, gnu)
     - ``datm_cdeps_control_cfsr``
     - | ``CICE_FIX/100``
       | ``CICE_IC/100``
       | ``DATM_CDEPS`` (files only)
       | ``DATM_CDEPS/CFSR/201110``
       | ``MOM6_FIX/100``
       | ``MOM6_FIX_DATM/100``
       | ``MOM6_IC/100/2011100100``
   * - ``datm_cdeps`` (intel)
     - ``datm_cdeps_control_gefs``
     - | ``CICE_FIX/100``
       | ``CICE_IC/100``
       | ``DATM_CDEPS`` (files only)
       | ``DATM_CDEPS/GEFS_NEW/201110``
       | ``MOM6_FIX/100``
       | ``MOM6_FIX_DATM/100``
       | ``MOM6_IC/100/2011100100``
   * - ``datm_cdeps`` (intel)
     - ``datm_cdeps_ciceC_cfsr``
     - | ``CICE_FIX/100``
       | ``CICE_IC/100``
       | ``DATM_CDEPS`` (files only)
       | ``DATM_CDEPS/CFSR/201110``
       | ``MOM6_FIX/100``
       | ``MOM6_FIX_DATM/100``
       | ``MOM6_IC/100/2011100100``
   * - ``datm_cdeps`` (intel)
     - ``datm_cdeps_gfs``
     - | ``CICE_FIX/100``
       | ``CICE_IC/100``
       | ``DATM_CDEPS`` (files only)
       | ``DATM_CDEPS/GFS/202103``
       | ``MOM6_FIX/100``
       | ``MOM6_FIX_DATM/100``
       | ``MOM6_IC`` (files only)
   * - ``datm_cdeps_land`` (intel)
     - ``datm_cdeps_lnd_gswp3``
     - | ``CPL_FIX/aC96o100``
       | ``DATM_GSWP3_input_data``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data/INPUT``
       | ``FV3_input_data/INPUT_L127_gfsv17``
       | ``NOAHMP_IC/CLMNCEP``
   * - ``datm_cdeps_land`` (intel)
     - ``datm_cdeps_lnd_era5``
     - | ``CPL_FIX/aC96o100``
       | ``DATM_ERA5_input_data_v2``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data/INPUT``
       | ``FV3_input_data/INPUT_L127_gfsv17``
       | ``NOAHMP_IC/ERA5``
   * - ``atm_ds2s_docn_pcice`` (intel, gnu)
     - ``atm_ds2s_docn_pcice``
     - | ``CICE_FIX/100``
       | ``CICE_IC/100``
       | ``CPL_FIX/aC96o100``
       | ``DOCN_DICE_ERA5``
       | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT``
       | ``FV3_input_data/INPUT_L127_gfsv17``
       | ``FV3_input_data/INPUT_L127_mx100``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``MOM6_FIX/100``
       | ``MOM6_IC`` (files only)
   * - ``atmaero`` (intel)
     - ``atmaero_control_p8``
     - | ``FV3_fix``
       | ``FV3_fix_tiled/C96mx100``
       | ``FV3_input_data`` (files only)
       | ``FV3_input_data/INPUT_L127_gfsv17``
       | ``FV3_input_data_INCCN_aeroclim/MERRA2_y14_24``
       | ``FV3_input_data_INCCN_aeroclim/aer_data/LUTS``
       | ``GOCART/p8``

The ``+container`` tests do not use ``GEFS``, ``FV3_input_data384``, ``LM4_input_data``, or ``FV3_regional``.
They also do not use the ``GFSv17opn_20251014`` and ``BM_IC-20220207`` data sets.

The data can be downloaded with the AWS CLI; no AWS account is needed when ``--no-sign-request`` is used.
For example:

.. code-block:: console

   export BUCKET=s3://noaa-ufs-regtests-pds/input-data-20260617
   export INPUTDATA_ROOT=/path/to/input-data-20260617

   aws s3 sync ${BUCKET}/FV3_fix/ ${INPUTDATA_ROOT}/FV3_fix/ --no-sign-request
   for combo in C96mx100 C192mx050 C192mx025 C48mx500 C24mx500 C24mx900 C12mx900; do
     aws s3 sync ${BUCKET}/FV3_fix_tiled/${combo}/ ${INPUTDATA_ROOT}/FV3_fix_tiled/${combo}/ --no-sign-request
   done
   # Skip subdirectories that the +container tests do not use
   aws s3 sync ${BUCKET}/MOM6_IC/ ${INPUTDATA_ROOT}/MOM6_IC/ --no-sign-request --exclude "025/*"
   aws s3 sync ${BUCKET}/CICE_IC/ ${INPUTDATA_ROOT}/CICE_IC/ --no-sign-request --exclude "025/*" --exclude "050/*"
   aws s3 sync ${BUCKET}/WW3_input_data_20260811/ ${INPUTDATA_ROOT}/WW3_input_data_20260811/ --no-sign-request

Repeat the ``aws s3 sync`` command for each remaining subdirectory in the table above.
Set ``INPUTDATA_ROOT`` (and, if needed, ``INPUTDATA_ROOT_WW3``) in line 4 of the platform definition file
(see :numref:`Section %s <container-rt-conf>`) to the download location.

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

The modulefiles required depend on the option selected. The container option
requires a user-adapted modulefile to load any host system modules needed at runtime.
The community platform option (native software stack) requires a modulefile that loads
all the required software stack libraries on the platform. All modulefiles are placed in the
``modulefiles/`` directory at the root of the repository.

.. _container-rt-modulefiles-container:

Container Option
~~~~~~~~~~~~~~~~

The container option uses two modulefiles. Only ``ufs_container.runtime.lua``
**must be adapted** by the user. The ``ufs_container.<compiler>.lua`` build module
depends on the software stack inside the container image and does not
require user changes. ``rt.sh`` checks that both files exist before it compiles or runs
a test in a container.

.. _container-rt-runtime-mod:

``ufs_container.runtime.lua`` — Host-Side Runtime Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This modulefile is loaded **on the host** by the compile and run job cards (or by
``rt.sh`` itself when no scheduler is used). **Users must create and adapt this file**
for their platform. Its content depends on the MPI launch method:

* **Slurm** (``srun``): ``srun`` coordinates MPI rank launch across compute nodes via the
  host Process Management Interface. The GNU-based image with OpenMPI 4.1.6 supports PMI2
  (``--mpi=pmi2``); the image with OpenMPI 5.0.7 supports PMIx (``--mpi=pmix``). Run
  ``srun --mpi=list`` to confirm availability. In this case no host MPI libraries are
  required and the modulefile only needs to load the Singularity/Apptainer module.

* **PBS** (``mpirun``/``mpiexec``): the host MPI launcher requires ABI-compatible MPI
  libraries on the host. Load the Singularity/Apptainer module together with compiler and
  MPI modules that match the container's toolchain.

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

Community Platform Option (Native Software Stack)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the container image field of the platform definition file is left blank, the model is
built and run natively on the host — no container is involved. A natively installed software
stack is expected to be present.

.. _container-rt-community-mod:

``ufs_<PLATFORM_NAME>.<compiler>.lua`` — Platform Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Users must create and adapt this file** for their platform. It is loaded during both
the compile and run stages and sets up the compiler toolchain, MPI library, and all
required software libraries available on the host system.

The file must be named ``modulefiles/ufs_<PLATFORM_NAME>.<compiler>.lua``, where
``<PLATFORM_NAME>`` matches the platform name in line 1 of the platform definition file and
``<compiler>`` is ``intel`` or ``gnu``. ``compile.sh`` loads this module for the build and
saves a copy next to the executable, which is then loaded for each test run.

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
The Platform Definition File
====================================

The platform definition file passed to ``rt.sh -P <platform.def>`` describes the platform: its name,
compiler, container image (if any), scheduler, run directory, and input data locations.
It does **not** list tests; tests are always selected from ``rt.conf`` (or the file given
with ``-l``), as described in :numref:`Section %s <container-rt-vs-rt>`.

The following platform definition files are provided in the ``tests/`` directory:

.. list-table::
   :widths: 25 60
   :header-rows: 1

   * - File
     - Description
   * - ``platform.def``
     - Template with a description of every field. Copy it and edit the copy for a new platform.
   * - ``stampede.def``
     - Example for GNU or Intel container runs on Stampede3 (TACC, University of Texas).

To set up a new platform:

#. Copy ``platform.def``, keeping the original as a reference, and edit the data lines of the
   copy for the new platform.
#. Tag the ``COMPILE`` and ``RUN`` lines to run on the platform with ``+<PLATFORM_NAME>`` in
   the **Machines** column of ``rt.conf`` (or of the file given with ``-l``), using the
   name given in line 1 of the platform definition file. ``-<PLATFORM_NAME>``
   excludes a line from the platform even if ``+<PLATFORM_NAME>`` is also present. With
   ``PLATFORM_NAME`` set to ``container``, the lines already tagged ``+container`` are used.
#. Run ``rt.sh`` with the new file:

   .. code-block:: console

      ./rt.sh -a <account> -P myplatform.def -l rt.conf

The file uses ``|`` as a field separator. Blank lines and lines starting with ``#`` are
skipped. ``rt.sh`` reads four required data lines, in the order described in
:numref:`Section %s <container-rt-conf-fields>`, and an optional fifth line that is
needed only for the Rocoto (``-r``) and ecFlow (``-e``) workflow managers. Tests are not
listed in this file; they are selected from ``rt.conf``.

By default, a run with a platform definition file is sequential and does not compare results
against baselines: it is a portability check that the model builds and runs on the platform,
not a regression test. The ``-c`` and ``-m`` options create and compare against baselines,
and ``-r`` and ``-e`` use Rocoto or ecFlow (these need line 5).

**Template** ``tests/platform.def``. The comments describe every entry, followed by the data
lines to edit:

.. literalinclude:: ../../../tests/platform.def
   :language: text

**Container option example** (Stampede3, GNU container). These are the data lines of
``tests/stampede.def``; the full file also contains the same entry descriptions as the
template:

.. code-block:: text

   #container | intel | /work2/10000/nperlin/stampede3/rocky9-oneapi2024.2-ss192.sif | /work,/work2,/scratch
   container | gnu | /work2/10000/nperlin/stampede3/rocky9-gcc13-ss192-ompi416.sif  | /work,/work2,/scratch
   48 | slurm | skx-dev | myqueue | mpirun
   /scratch/10000/nperlin/UFS-WM/RUNDIR_RT
   /work2/10000/nperlin/stampede3/ufs-wm_input/input-data-20260617 |  |  |
   #slurm | module load rocoto/1.3.7

In this example:

* Line 1 selects the GNU container image and binds ``/work``, ``/work2``, and
  ``/scratch`` into the container. To use the Intel container instead, comment out this
  line and uncomment the ``intel`` line above it.
* Line 2 gives 48 MPI tasks per node (Stampede3 SKX nodes) and the Slurm partition
  ``skx-dev``.
* Line 3 is the run directory, and line 4 gives only ``INPUTDATA_ROOT``; the
  other input data directories use their defaults.
* Line 5 is commented out, so the tests run sequentially. Uncomment it to use ``-r``.

The container software must be available before ``rt.sh`` is started (on Stampede3:
``module load tacc-apptainer``; see :numref:`Section %s <container-rt-apptainer>`).

.. _container-rt-conf-fields:

Platform Definition Fields
--------------------------

**Line 1:** ``PLATFORM_NAME | RT_COMPILER | CONTAINER_IMG | BIND_DIRS``

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Field
     - Description
   * - ``PLATFORM_NAME``
     - Required. Becomes ``MACHINE_ID`` for the run, and is the name used in the
       ``+<PLATFORM_NAME>``/``-<PLATFORM_NAME>`` tags in ``rt.conf``. Use ``container`` for
       container runs: this selects the ``+container`` tests in ``rt.conf`` and the container
       job card templates in ``tests/fv3_conf/``. For a native stack, use a name that matches
       the platform modulefile (see :numref:`Section %s <container-rt-community-mod>`).
   * - ``RT_COMPILER``
     - Required. The platform's only compiler: ``intel`` or ``gnu``. ``rt.conf`` lines that
       are tagged for the platform but use the other compiler are skipped with a notice.
   * - ``CONTAINER_IMG``
     - Absolute path to the Singularity/Apptainer image file (``*.sif``) on the host.
       Leave blank to build and run with a native software stack instead, using the
       modulefile ``modulefiles/ufs_<PLATFORM_NAME>.<compiler>.lua``, which may need to be
       adapted for the user's platform. If the image file does
       not exist, the lines that need it are skipped.
   * - ``BIND_DIRS``
     - Comma-separated list of host directories to bind/mount to the container.
       Include all filesystems containing the source tree, input data, and run directory.
       See :numref:`Section %s <container-rt-binddirs>` for typical values on Tier 1 platforms.
       Ignored for a native stack.

**Line 2:** ``TPN | SCHEDULER | PARTITION | QUEUE | MPI_LAUNCH``

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Field
     - Description
   * - ``TPN``
     - Required. Tasks per node: the number of MPI tasks per node on this platform
       (a numeric value, e.g., ``48`` on Stampede3).
   * - ``SCHEDULER``
     - Required. ``slurm``, ``pbs``, or ``none`` for interactive runs without a scheduler.
       This field cannot be left blank.
   * - ``PARTITION``
     - Slurm partition name (leave blank for PBS or ``none``).
   * - ``QUEUE``
     - Slurm QOS / PBS queue name (leave blank for ``none``).
   * - ``MPI_LAUNCH``
     - MPI launch command used when ``SCHEDULER`` is ``none`` and in PBS job cards:
       ``mpirun`` or ``mpiexec``. Defaults to ``mpirun`` if omitted.

The scheduler account (project) is not set in this file. It is always given on the command
line with ``rt.sh -a <account>``.

**Line 3:** ``RUNDIR_ROOT``

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Field
     - Description
   * - ``RUNDIR_ROOT``
     - Required. Top-level directory where compile and test run directories are created.
       This should be a user-writable path, preferably on a scratch or work filesystem.
       ``rt.sh`` uses this directory exactly as given and does not delete it at the end of
       the run. A symlink ``tests/run_dir`` pointing to it is created.

**Line 4:** ``INPUTDATA_ROOT | INPUTDATA_ROOT_WW3 | INPUTDATA_LM4 | INPUTDATA_GFSv17opn``

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Field
     - Description
   * - ``INPUTDATA_ROOT``
     - Required. Input data directory, e.g., ``.../input-data-20260617``
       (see :numref:`Section %s <container-rt-data>`).
   * - ``INPUTDATA_ROOT_WW3``
     - Optional. WaveWatch III input data directory. Defaults to
       ``${INPUTDATA_ROOT}/WW3_input_data_20260811`` if left blank.
   * - ``INPUTDATA_LM4``
     - Optional. LM4 land model input data directory. Defaults to
       ``${INPUTDATA_ROOT}/LM4_input_data`` if left blank.
   * - ``INPUTDATA_GFSv17opn``
     - Optional. GFS v17 operational input data directory, needed only by tests that
       use it. None of the ``+container`` tests need it.

**Line 5** (optional; needed only for ``-r`` or ``-e``): ``ROCOTO_SCHEDULER | WORKFLOW_MODULE_CMD``

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Field
     - Description
   * - ``ROCOTO_SCHEDULER``
     - Rocoto's name for the scheduler, which is not always the same as ``SCHEDULER``
       (e.g., PBS Professional is ``pbspro``). Required to use ``-r``.
   * - ``WORKFLOW_MODULE_CMD``
     - A shell command that puts the Rocoto or ecFlow command-line tools on ``PATH``
       (e.g., ``module load rocoto``). Leave blank if they are already on ``PATH``.

.. _container-rt-run:

====================
Running the Tests
====================

Tests are launched by running ``rt.sh`` with ``-P`` from the ``tests/`` directory:

.. code-block:: console

   cd ${WM_HOME}/tests
   ./rt.sh -a <account> -P <platform.def> [options]

For example, to run all ``+container`` tests in ``rt.conf`` with the Stampede3 example:

.. code-block:: console

   ./rt.sh -a <account> -P stampede.def -l rt.conf

Command-Line Options
--------------------

The usual ``rt.sh`` options can be combined with ``-P``. The ones most useful for container
and community platform runs are:

.. list-table::
   :widths: 15 60
   :header-rows: 1

   * - Option
     - Description
   * - ``-a <account>``
     - Scheduler account/project. Always required.
   * - ``-P <platform.def>``
     - Platform definition file (see :numref:`Section %s <container-rt-conf>`).
   * - ``-l <file>``
     - Use ``<file>`` instead of ``rt.conf`` to select tests.
   * - ``-n "<test> <compiler>"``
     - Run a single test, for example ``-n "control_c48 intel"``. Its compile still runs.
   * - ``-s <file>``
     - Run only the subset of tests listed in ``<file>``.
   * - ``-o``
     - Compile only; skip all tests.
   * - ``-x``
     - Dry run. With ``-P``, ``rt.sh`` still compiles; for each test it then checks the
       compiled executable, stages the input data, checks the container (if any), and prepares
       the job card, but does not submit the job. Results are reported as
       ``DRY RUN SUCCESS``/``DRY RUN FAIL``.
   * - ``-c``
     - Create a baseline under ``${RUNDIR_ROOT}/REGRESSION_TEST``.
   * - ``-m``
     - Compare against the baseline previously created with ``-c`` under
       ``${RUNDIR_ROOT}/REGRESSION_TEST``.
   * - ``-r`` / ``-e``
     - Use the Rocoto or ecFlow workflow manager instead of running sequentially.
       Requires line 5 of the platform definition file.
   * - ``-d``
     - Delete run directories that are not used by other tests.
   * - ``-v``
     - Verbose output (shell tracing).
   * - ``-h``
     - Print help and exit.

Without ``-c`` or ``-m``, results are not compared against any baseline.

.. note::

   If the executable ``tests/fv3_<compile_name>_<compiler>.exe`` from an earlier run is
   already present, ``rt.sh -P`` reuses it and skips that compile. Remove the executable
   to force a rebuild.

   Because ``-P`` runs do not use the ``rt.sh`` lock, several runs can proceed at the same
   time (e.g., one per compiler). Give each platform definition file its own ``RUNDIR_ROOT``
   so that the runs do not share directories.

Job Script Templates
--------------------

When ``SCHEDULER`` is ``slurm`` or ``pbs``, compile and test jobs are created from job
script templates in ``tests/fv3_conf/``, selected by scheduler and ``PLATFORM_NAME``:

- ``compile_slurm.IN_<PLATFORM_NAME>`` — Slurm compile job card
- ``fv3_slurm.IN_<PLATFORM_NAME>`` — Slurm run job card
- ``compile_qsub.IN_<PLATFORM_NAME>`` — PBS compile job card
- ``fv3_qsub.IN_<PLATFORM_NAME>`` — PBS run job card

Templates for the container option (``PLATFORM_NAME=container``) are provided in the
repository. For a native stack with a scheduler, users must create templates for their
platform name, following the container or Tier 1 templates in ``tests/fv3_conf/``.
When ``SCHEDULER`` is ``none``, no template is used and compile and run steps execute
directly on the current host.

These templates contain scheduler directives and environment setup that may require
platform-specific adjustments, such as partition or queue names, wall-clock limits, node
counts, and any platform-specific environment variables required before launching the model.

Running with a Job Scheduler (Slurm or PBS)
--------------------------------------------

When ``SCHEDULER`` is ``slurm`` or ``pbs``, set ``PARTITION`` and ``QUEUE`` in the platform
definition file and run ``rt.sh`` from a login node:

.. code-block:: console

   ./rt.sh -a <account> -P <platform.def> -l rt.conf

By default, ``rt.sh`` submits each compile and test job to the scheduler and waits for it to
finish before submitting the next one. With ``-r`` or ``-e``, the jobs are managed by Rocoto
or ecFlow instead.

Running Interactively (No Scheduler)
--------------------------------------

When ``SCHEDULER`` is ``none``, jobs run directly on the current host —
suitable for an allocated compute node or single-workstation development.

If login nodes may not be used for running tests, request an interactive compute node
allocation before running ``rt.sh``. On **Slurm** systems the command may look similar to:

.. code-block:: console

   salloc -N 1 -n <cores> -A <account> -t <time> -q <qos> --partition=<partition>

On **PBS** systems the command may look similar to:

.. code-block:: console

   qsub -I -l walltime=<time> -A <account> -q <queue> -l select=1:ncpus=<cores>:mpiprocs=<cores>

After the allocation is granted (and connecting via ``ssh`` to the compute node if required),
run ``rt.sh``. The model is started with the ``MPI_LAUNCH`` command from line 2 of the
platform definition file.

.. code-block:: console

   cd tests
   ./rt.sh -a <account> -P <platform.def> -l rt.conf

.. note::

   For the **container option**, ``rt.sh`` starts the software container first, and then
   runs the MPI tasks entirely inside it, which is the correct approach for single-node
   interactive runs. Interactive container runs require ``PLATFORM_NAME=container``.

   The ``--mpi=pmi2`` flag is a Slurm ``srun``-specific option and should **not** be used
   with ``mpirun`` or ``mpiexec``. When ``SCHEDULER`` is ``none``, ``srun`` is not used.

.. _container-rt-output:

=============================
Run Directory and Log Files
=============================

``rt.sh`` creates the following structure under ``RUNDIR_ROOT``:

**Container option** (``PLATFORM_NAME=container``):

.. code-block:: text

   ${RUNDIR_ROOT}/
   ├── compile_<compile_name>_<compiler>/   # compile working directory
   │   ├── job_card                         # compile job script
   │   ├── modulefiles/                     # modulefiles staged for the build
   │   └── out / err                        # job stdout and stderr files
   ├── <test_name>_<compiler>/              # test working directory
   │   ├── job_card                         # test job script (scheduler runs)
   │   ├── fv3_container_run.sh             # script run inside the container (no scheduler)
   │   ├── modulefiles/                     # modulefiles staged for the run
   │   └── out / err                        # job stdout and stderr files
   └── REGRESSION_TEST/                     # baseline, only with -c

**Community platform option** (native software stack):

.. code-block:: text

   ${RUNDIR_ROOT}/
   ├── compile_<compile_name>_<compiler>/   # compile working directory
   │   ├── job_card                         # compile job script (scheduler runs)
   │   └── out / err                        # job stdout and stderr files
   ├── <test_name>_<compiler>/              # test working directory
   │   ├── job_card                         # test job script (scheduler runs)
   │   ├── fv3_run.sh                       # native run script (no scheduler)
   │   ├── modulefiles/                     # modulefile staged for the run
   │   └── out / err                        # job stdout and stderr files
   └── REGRESSION_TEST/                     # baseline, only with -c

Because ``RUNDIR_ROOT`` is reused from run to run, a compile or test directory left from an
earlier run is not deleted. It is renamed to ``<directory>_old_<YYYYMMDDHHMM>`` before the
new run starts. A symlink ``tests/run_dir`` points to ``RUNDIR_ROOT``.

Log files are written to the ``tests/logs/`` directory, as for Tier 1 runs
(see :numref:`Section %s <log-files>`):

* ``RegressionTests_<PLATFORM_NAME>.log`` — summary of the run. For a default run
  without ``-c`` or ``-m``, the file is named ``RegressionTests_weekly_<PLATFORM_NAME>.log``
  because the comparison step is skipped.
* ``log_<PLATFORM_NAME>/`` — detailed compile and test logs.

For sequential runs, a ``COMPILE/TEST SUMMARY`` with PASS/FAIL for each compile and test
is also printed to the terminal when all tests have finished.
