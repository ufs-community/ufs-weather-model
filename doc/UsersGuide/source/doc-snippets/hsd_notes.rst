.. note::
   
   Once the tests run successfully with the ``-c`` option (baseline created), users can compare future test results with the newly created baseline using ``-m`` instead of ``-c``.

For further test management, users may save the test directory location in an environment variable:

.. code-block:: console

   export UFS_WM_TEST=/path/to/expt_dirs/ufs_test