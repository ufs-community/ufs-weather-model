.. _add-hsd-case:

*******************************************************************
General Guidelines for Adding a Test Case to the UFS HSD Framework
*******************************************************************

A user wishing to add a test case to the Unified Forecast System (UFS) Hierarchical System Development (HSD) framework via the ``tests-dev`` directory will generally adhere to the following procedure to ensure the proper integration of the test case inside the UFS Weather Model (WM).

After selecting the case of interest (e.g., idealized case, model improvement, process/phenomenon/component isolation, etc.), the user can either utilize an existing model configuration (e.g., ``ATM``, ``S2S``, ``S2SWA``, ``HAFSW``) if the case can make use of a pre-existing model configuration or bring in their own executable. However, if a new model configuration is needed, the code changes necessary to introduce the new model configuration/executable are the responsibility of the user to add to the UFS-WM. Therefore, the process to introduce a new model configuration is beyond the scope of this document.

For the cases currently in the UFS HSD, only the atmospheric model is utilized, and the ``ATM`` model configuration option is used. This is called at compile time in the CMake build system command with ``-DAPP=ATM``. A tropical cyclone case (currently in development) will utilize the ``HAFSW`` model configuration (i.e., ``-DAPP=HAFSW``). Therefore, cases that can utilize the ``ATM`` or ``HAFSW`` model configurations serve as good starting points.

================================
General Steps to Add a Test Case
================================

To add a test case to be run with ``ufs_test.sh`` in ``tests-dev``, the user must update ``ufs_test.yaml``. The new case should be added under the ``tests:`` section of the model configuration that corresponds to the specific model configuration/executable it will be run with.

For example, to add an “atmosphere-only” case using “FV3” and the “Intel” compiler, users would add their new test case name in the ``ufs_test.yaml`` file by updating the ``tests:`` section under ``atm_dyn32_intel`` configuration. For these new test cases, set the ``project`` option to ``test_case`` and ``baseline`` option to ``True``. For an Intel example, please refer to ``2020_CAPE`` and ``baroclinic_wave`` test cases under ``ufs_test.yaml#L108``. Users can also add a GNU-based test (see where the ``baroclinic_wave`` test was added).

A test case file corresponding to each of these new HSD test cases needs to be added under:

  https://github.com/ufs-community/ufs-weather-model/tree/develop/tests-dev/test_cases/tests/

This test case file exports model, I/O, resource settings, and several other parameters to the ``ufs_test.sh`` workflow and sets the CCPP suite, input namelist, model configuration, diagnostic table, field table, and ``FV3_run`` files that will be used in conjunction with the HSD test cases. Refer to Section 5, *Configuration Files*, for details on updating the configuration files corresponding to these HSD test cases.

=====================
Model Compile Options
=====================

Model compile options are set in ``ufs_test.yaml``. Default compile options are set in the ``options`` section (e.g., model mode such as ``ATM``), CCPP suites, etc. Users can add compile options here. If a physics suite needs to be added for the case, it should be added to the ``CCPP_SUITES`` section of the compile options.

In addition to setting the model compile options, users should consider which compiler they wish their model to be built with. Intel and GNU options are generally supported in the UFS-WM. To enable compilation of the model with both Intel and GNU, ensure that the new HSD test case is added under the appropriate section in both the Intel and GNU components of the ``ufs_test.yaml`` file.

========
Physics
========

UFS-WM supports options from the CCPP. Users can select a suite to use for their case that properly captures the physical processes and parametrizations required for the specific scientific case. The suite should be included in the compile options (see section *Model Compile Options*) and can be exported in the test case file. This is then added to the case's input namelist file, if it is a dynamic variable.

Users can bring in their own physics suites, but would need to submit a PR to the CCPP repository before this can be incorporated into the UFS-WM.

===========
Input Data
===========

The user is responsible for generating, staging, and properly pointing all configuration files to the appropriate input data for the specific HSD test case (e.g., initial/boundary conditions, restart files, etc.). Input data copying from fixed directories on RHDPCS platforms into the test case run directory is handled in the workflow via the ``FV3_run`` file, so this file should include lines to copy the necessary data into the run directory.

- Input data such as ICs/BCs should be copied into the ``INPUT`` directory of the run directory.
- Restart files (for warm starts) should be copied into the ``RESTART`` directory.

Input data can be generated via ``UFS_UTILS`` or the user’s own scripts, as long as it matches the format/variable structure required by the UFS-WM component(s) that the HSD test case will utilize (e.g., for atmosphere-only cases, GFS data such as ``gfs.*.nc``, ``sfc*.nc``, ``orog*.nc``, and gridding files are required).

=====================
Configuration Files
=====================

Input namelist, model configure, UFS configuration, diagnostic table, and field tables should be added in conjunction with the HSD test case. At minimum, a namelist and ``FV3_run`` file should be included that is appropriate for the test case. The ``FV3_run`` file is used by the workflow to read in exported values from the tests file and properly stage input/restart data and fix files.

Default model configurations, UFS configurations, diagnostic table, and field tables may be used, but the user should confirm that these are appropriate for their specific HSD test case. Users may reference other UFS-WM regression test files and namelists, model configurations, UFS configurations, diagnostic/field tables, and ``FV3_run`` files for help in setting up their own files.

- Namelists, model configuration, diagnostic/field tables, and UFS configuration files should be added to:

  https://github.com/ufs-community/ufs-weather-model/tree/develop/tests-dev/test_cases/parm

- The ``FV3_run`` file should go in:

  https://github.com/ufs-community/ufs-weather-model/tree/develop/tests-dev/test_cases/exp_conf

============
Fix Files
============

The user is responsible for staging all fix files needed for the case. The ``FV3_run`` file should include lines to copy the necessary fix files from fixed directories on RHDPCS into the HSD test case run directory.

See ``*.IN`` files in:

  https://github.com/ufs-community/ufs-weather-model/tree/develop/tests-dev/test_cases/exp_conf

or:

  https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/fv3_conf/control_run.IN

for examples of how this is done.

==========================
Running the Test Case
==========================

With the test case added to ``ufs_test.yaml``, all necessary configuration files in place, any source code changes made to relevant subcomponents, and input and fix file data staged, the user should be able to perform a test run of their case from the ``tests-dev`` directory using:

.. code-block:: console

   ./ufs_test.sh -a <account> -s -r -c -n "<test_name> <compiler_option>"
