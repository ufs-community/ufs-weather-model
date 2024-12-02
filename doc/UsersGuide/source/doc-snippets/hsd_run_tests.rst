Launch tests from the ``${UFS_WM}/tests-dev`` directory with the following command:

.. code-block:: console

   cd ${UFS_WM}/tests-dev
   ./ufs_test.sh -a <ACCOUNT> [-s] [-c] -k -r -n "<CASE_NAME> <COMPILER>"

where:

* ``<ACCOUNT>``: Account/project number for batch jobs.
* ``<CASE_NAME>``: Name of the test case (e.g., ``2020_CAPE`` or ``baroclinic_wave``).
* ``<COMPILER>``: Compiler used for the tests (``intel`` or ``gnu``).

**Command-line Options:**

* ``-s``: Syncs scripts from ``./ufs-wm/tests`` to ``./ufs-wm/tests-dev`` (only required on the first run)
* ``-c``: Creates a new baseline (necessary until idealized case baselines are staged in the ``UFS_WM_RT`` directory).  
* ``-k``: Keeps runtime directories after test completion
* ``l``: Runs test cases listed in a YAML file
* ``-m``: Compares against existing baseline results (baseline must exist)
* ``-n``: Runs a single test case
* ``-r``: Uses Rocoto workflow manager

.. note::

   After the initial run of ``ufs_test.sh`` with the ``-s`` option, users do not need to use ``-s`` again unless they subsequently alter files inside of ``tests-dev/test_cases``. If files have changed, the user will need to rerun ``ufs_test.sh`` with ``-s`` or work from a fresh clone so that everything is properly copied via ``-s``.