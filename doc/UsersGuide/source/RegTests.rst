.. _RegTests:
  
******************************************
Regression Tests for Development
******************************************

* Motivation

* Description of the default tests

* Directory structure for inputs and outputs

* How to run a subset of tests

* Baselines: existing baselines, creating new baselines, comparing against a new baseline

* Running regression tests

* Checking the results

* Adding a regression test


System Requirements, Libraries, and Compilers
---------------------------------------------
The build system for the UFS with CCPP relies on the use of the Python scripting language, along with ``cmake``.

The basic requirements for building and running the UFS with CCPP are listed below. The versions listed reflect successful tests and there is no guarantee that the code will work with different versions.

    * FORTRAN 90+ compiler versions
        * ifort 15.1.133, 18.0.1.163 and 19.0.2
        * gfortran 6.2, 8.1, and 9.1
    * C compiler versions
        * icc v18.0.1.163 and 19.0.2
        * gcc 6.2.0 and 8.1
        * AppleClang 10.0
    * MPI job scheduler versions
        * mpt 2.19
        * impi 5.1.1.109 and 5.1.2.150
        * mpich 3.2.1
    * cmake versions 2.8.12.1, 2.8.12.2, and 3.6.2
    * netCDF with HDF5, ZLIB and SZIP versions 4.3.0, 4.4.0, 4.4.1.1, 4.5.0, 4.6.1, and 4.6.3 (not 3.x)
    * Python versions 2.7.5, 2.7.9, and 2.7.13 (not 3.x)

A number of NCEP libraries are required to build and run FV3 and are listed in :numref:`Table %s <NCEP_lib_FV3>`.

.. _NCEP_lib_FV3:

.. table:: *NCEP libraries required to build the UFS Atmosphere*

    +---------------------------+-------------+----------------------------------------------------+
    | Library                   | Version     | Description                                        |
    +===========================+=============+====================================================+
    | bacio                     | 2.0.1       | NCEP binary I/O library                            |
    +---------------------------+-------------+----------------------------------------------------+
    | ip                        | 2.0.0/3.0.0 | NCEP general interpolation library                 |
    +---------------------------+-------------+----------------------------------------------------+
    | nemsio                    | 2.2.3       | NEMS I/O routines                                  |
    +---------------------------+-------------+----------------------------------------------------+
    | sp                        | 2.0.2       | NCEP spectral grid transforms                      |
    +---------------------------+-------------+----------------------------------------------------+
    | w3emc                     | 2.2.0       | NCEP/EMC library for decoding data in GRIB1 format |
    +---------------------------+-------------+----------------------------------------------------+
    | w3nco/v2.0.6              | 2.0.6       | NCEP/NCO library for decoding data in GRIB1 format |
    +---------------------------+-------------+----------------------------------------------------+

These libraries are prebuilt on most NOAA machines using the Intel compiler. For those needing to build the libraries themselves, GMTB recommends using the source code from GitHub at https://github.com/NCAR/NCEPlibs.git, which includes build files for various compilers and machines using OpenMP flags and which are thread-safe.

In addition to the NCEP libraries, some additional external libraries are needed (:numref:`Table %s <ext_lib_FV3>`).

.. _ext_lib_FV3:

.. table:: *External libraries necessary to build the UFS Atmosphere*

    +--------------------+-------------------------+---------------------------------------------------------------------------------------------+
    | Library            | Version                 | Description                                                                                 |
    +====================+=========================+=============================================================================================+
    | ESMF               | V7.1.0r and v8.0.0_bs21 | Earth System Modeling Framework for coupling applications                                   |
    +--------------------+-------------------------+---------------------------------------------------------------------------------------------+
    | netCDF             | 4.3.0 and 4.6.1         | Interface to data access functions for storing and retrieving data arrays                   |
    +--------------------+-------------------------+---------------------------------------------------------------------------------------------+
    | SIONlib (optional) | v1.7.2                  | Parallel I/O library (link) that can be used to read precomputed lookup tables instead of \ |
    |                    |                         | computing them on the fly (or using traditional Fortran binary data files)                  |
    +--------------------+-------------------------+---------------------------------------------------------------------------------------------+

The Earth System Modeling Framework (ESMF), the SIONlib, the NCEPlibs, and the netCDF libraries must be built with the same compiler as the other components of the UFS Atmosphere.

Building the UFS Atmosphere
---------------------------

A complete listing and description of the FV3 build options were discussed elsewhere. This section will describe the commands needed to build the different options using the script ``compile.sh`` provided in the NEMSfv3gfs distribution. This script calls ``ccpp_prebuild.py``, so users do not need to run the *prebuild* step manually. All builds using ``compile.sh`` are made from the ``./tests`` directory of NEMSfv3gfs and follow the basic command:

.. code-block:: console

    ./compile.sh $PWD/../FV3 system.compiler 'MAKEOPTS'

Here, ``system`` stands for the machine on which the code is compiled and can be any of the following machines and compilers: *theia, jet, cheyenne, gaea, stampede, wcoss_cray, wcoss_dell_p3, supermuc_phase2, macosx*, or *linux*.

``compiler`` stands for the compiler to use and depends on the system. For *theia* and *cheyenne*, the available options are ``intel`` and ``gnu``. For *macosx* and *linux*, the only tested compiler is ``gnu``. For all other platforms, ``intel`` is the only option at this time.

The ``MAKEOPTS`` string, enclosed in single or double quotes, allows to specify options for compiling the code. The following options are of interest for building the CCPP version of NEMSfv3gfs:

* **CCPP=Y** - enables :term:`CCPP` (default is ``N``)
* **STATIC=Y** - enables the CCPP static mode; requires ``CCPP=Y`` (default is ``N``) and ``SUITES=...`` (see below)
* **SUITES=XYZ, ABC, DEF, ...** - specify SDF(s) to use when compiling the code in CCPP static mode; SDFs are located in ``ccpp/suites/``, omit the path in the argument; requires ``CCPP=Y STATIC=Y`` (default is ``‘’``)
* **SION=Y** - enables support for the SIONlib I/O library (used by CCPP to read precomputed lookup tables instead of computing them on the fly); available on *Theia, Cheyenne, Jet*; also available on *Mac OS X* and *Linux* if instructions in ``doc/README_{macosx,linux}.txt`` are followed (default is ``N``)
* **32BIT=Y** - compiles FV3 dynamical core in single precision; note that physics are always compiled in double precision; this option is only available on *Theia, Cheyenne*, and *Jet* (default is ``N``)
* **REPRO=Y** - compiles code in :term:`REPRO` mode, i.e. removes certain compiler optimization flags used in the default :term:`PROD` mode to obtain bit-for-bit (b4b) identical results between CCPP and non-CCPP code (default is ``N``)
* **DEBUG=Y** - compiles code in DEBUG mode, i.e. removes all optimization of :term:`PROD` mode and add bound checks; mutually exclusive with ``REPRO=Y`` (default is ``N``)
* **TRANSITION=Y** - applies selective lowering of optimization for selected files to obtain b4b with non-CCPP code in PROD mode (only when using Intel 15 on *Theia*)

Examples:

* Compile non-CCPP code with 32-bit dynamics on *Theia* with the Intel compiler

    .. code-block:: console

        ./compile.sh $PWD/../FV3 theia.intel ‘32BIT=Y’

* Compile dynamic CCPP code in ``DEBUG`` mode on *Jet*

    .. code-block:: console

        ./compile.sh $PWD/../FV3 jet.intel ‘CCPP=Y DEBUG=Y’

* Compile static CCPP code for the CPT suite on *Linux* with the GNU compiler, enable support for the SIONlib I/O library (requires that the library to be installed)

    .. code-block:: console

        ./compile.sh $PWD/../FV3 linux.gnu ‘SION=Y CCPP=Y STATIC=Y SUITES=FV3_CPT_v0’

* *Cheyenne* static build with multiple suites:

    .. code-block:: console

        ./compile.sh $PWD/../FV3 cheyenne.intel ‘CCPP=Y STATIC=Y SUITES=FV3_GFS_v15,FV3_CPT_v0’


Running the UFS Atmosphere Using the Regression Tests (RTs)
------------------------------------------------------------

Regression testing is the process of testing changes to the programs to make sure that the existing functionalities still work when changes are introduced. By running the RTs (or a subset of them by copying a RT configuration file and editing it), the code is compiled, the run directories are set up, and the code is executed. The results are typically compared against a pre-existing baseline, but in certain occasions it is necessary to first create a new baseline (for example, in a new platform where a baseline does not exist or when it is expected that a new development will change the answer). Because the RTs set up the run directories, this is a useful and easy way to get started, since all the model configuration files and necessary input data (initial conditions, fixed data) are copied into the right place.

Overview of the RTs
^^^^^^^^^^^^^^^^^^^

The RT configuration files are located in ``./tests`` relative to the top-level directory of NEMSfv3gfs and have names ``rt*.conf``. The default RT configuration file, supplied with the NEMSfv3gfs master is called ``rt.conf`` and runs four types of configurations: IPD PROD, IPD REPRO, CCPP PROD, and CCPP REPRO. For the IPD configurations, CCPP is not used, that is, the code is compiled with ``CCPP=N``. The PROD configurations use the compiler flags used in NCEP operations for superior performance, while the REPRO configurations remove certain compiler flags to create b4b identical results between CCPP and IPD configurations. Before running the RT script ``rt.sh`` in directory ``./tests``, the user has to set some environment variables on the working shell: ``ACCNR`` (account to be charged for running the RTs), ``NEMS_COMPILER`` (optional for the ``intel`` compiler option, set to ``gnu`` to switch), and potentially ``RUNDIR_ROOT`` (location for the RT run directories), underneath which directories called ``rt_$PID`` are created (``$PID`` is the process identifier of the ``rt.sh`` invocation). This may be required on systems where the user does not have write permissions in the default run directory tree.

.. code-block:: console

    export ACCNR=...
    export NEMS_COMPILER=intel
    export RUNDIR_ROOT=/full/path/under/which/rt_$PID/will/be/created

Running the full default RT suite defined in ``rt.conf`` using the script ``rt.sh``:

.. code-block:: console

    ./rt.sh -f

This command can only be used on a NOAA machine using the Intel compiler, where an *official baseline* is available. For information on testing the CCPP code, or using alternate computational platforms, see the following sections.

This command and all others below produce log output in ``./tests/log_machine.compiler``. These log files contain information on the location of the run directories that can be used as templates for the user. Each ``rt*.conf`` contains one or more compile commands preceding a number of tests.


Baselines
^^^^^^^^^^^^^^^^^^^

Regression testing is only possible on machines for which baselines exist. EMC maintains *official baselines* on *Theia* and *Wcoss* created with the Intel compiler. GMTB maintains additional baselines on *Jet*, *Cheyenne*, and *Gaea*. While GMTB is trying to keep up with changes to the official repositories, baselines maintained by GMTB are not guaranteed to be up-to-date.

When porting the code to a new machine, it is useful to start by establishing a *personal baseline*. Future runs of the RT can then be compared against the *personal baseline* to ascertain that the results have not been inadvertently affected by code developments. The ``rt.sh -c`` option is used to create a *personal baseline*.

.. code-block:: console

    ./rt.sh -l rt.conf -c fv3 # create own reg. test baseline

Once the *personal baseline* has been created, future runs of the RT should be compared against the *personal baseline* using the ``-m`` option.

.. code-block:: console

    ./rt.sh -l rt.conf -m # compare against own baseline

The script rt.sh
^^^^^^^^^^^^^^^^^^^

``rt.sh`` is a bash shell file to run the RT and has the following options:

.. code-block:: console

    Usage: $0 -c <model> | -f | -s | -l <file> | -m | -r | -e | -h
    -c  create new baseline results for <model>
    -f  run full suite of regression tests
    -s  run standard suite of regression tests
    -l  run test specified in <file>
    -m  compare against new baseline results
    -r  use Rocoto workflow manager
    -e  use ecFlow workflow manager
    -h  display this help

The location of the run directories and *personal baseline* directories is controlled in ``rt.sh`` on a per-machine basis. The user is strongly advised to NOT modify the path to the *official baseline* directories.

The *official baseline* directory is defined as:

.. code-block:: console

    RTPWD=$DISKNM/trunk-yyyymmdd/${COMPILER} # on Cheyenne
    RTPWD=$DISKNM/trunk-yyyymmdd             # elsewhere

Note that ``yyyymmdd`` is the year, month and day the baseline was created using top of master code.

.. warning::  Modifying ``$DISKNM`` will break the RTs!

*Personal baseline* results (see below) are stored in

.. code-block:: console

    NEW_BASELINE=${STMP}/${USER}/FV3_RT/REGRESSION_TEST

and RTs are run in ``$RUNDIR_ROOT``.

Example: *Theia*

.. code-block:: console

    ...
    dprefix=/scratch4/NCEPDEV
    DISKNM=$dprefix/nems/noscrub/emc.nemspara/RT
    STMP=$dprefix/stmp4
    PTMP=$dprefix/stmp3
    ..

In case a user does not have write permissions to ``$STMP (/scratch4/NCEPDEV/stmp4/)``, ``$STMP`` must be modified without modifying ``$DISKNM`` (i.e. ``dprefix``). Similarly, if the user does not have write permissions to ``$PTMP``, the user can set the ``$RUNDIR_ROOT`` environment variable to change the location of the run directories as described below.

.. code-block:: console

    # Overwrite default RUNDIR_ROOT if environment variable RUNDIR_ROOT is set
    RUNDIR_ROOT=${RUNDIR_ROOT:-${PTMP}/${USER}/FV3_RT}/rt_$$


Compatibility between the Code Base, the SDF, and the Namelist in the UFS Atmosphere
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The variable ``suite_name`` within the ``namelist.input`` file used in the UFS Atmosphere determines which suite will be employed at run time (e.g., ``suite_name=FV3_GFS_v15``). It is the user’s responsibility to ascertain that the other variables in ``namelist.input`` are compatible with the chosen suite. When runs are executed using the RT framework described in the preceding sections, compatibility is assured. For new experiments, users are responsible for modifying the two files (``SDF`` and ``namelist.input``) consistently, since limited checks are in place.

Information about the UFS Atmosphere physics namelist can be found with the CCPP Scientific Documentation at https://dtcenter.org/GMTB/v4.0/sci_doc/.
