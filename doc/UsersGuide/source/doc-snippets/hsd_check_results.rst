When the test case finishes running, users should see console output that includes a ``SUCCESS`` message. For example: 

.. code-block:: console
   :emphasize-lines: 2 

   Performing Cleanup...
   REGRESSION TEST RESULT: SUCCESS
   + echo 'ufs_test.sh finished'
   ufs_test.sh finished
   + cleanup
   ++ awk '{print $2}'
   + PID_LOCK=2133541
   + [[ 2133541 == \2\1\3\3\5\4\1 ]]
   + rm -rf /scratch2/NAGAPE/epic/User.Name/ufs-weather-model/tests-dev/lock
   + [[ false == true ]]
   + trap 0
   + exit

Compilation and model run directories can be accessed in the local repository via the ``run_dir`` softlink, which points to the actual ``FV3_RT`` directory. Each test generates ``atm*.nc`` and ``sfc*.nc`` files at specified forecast hour intervals. 

Users can view progress of compile or model run phases by using the ``tail -f <file>`` command or ``vi``/``vim`` on the ``err`` or ``out`` files in the ``run_dir/compile*`` or ``run_dir/<case_name>`` directories. 