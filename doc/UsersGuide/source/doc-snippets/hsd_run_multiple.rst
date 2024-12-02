To run multiple cases at once, modify ``ufs_test.yaml`` to contain only a subset of tests (e.g., ``2020_CAPE`` and ``baroclinic_wave``) and use the ``-l`` argument:

.. code-block:: console

   ./ufs_test.sh -a epic -s -c -k -r -l ufs_test.yaml

Alternatively, users may copy the sections for ``2020_CAPE``/``baroclinic_wave`` tests into a new YAML file (e.g., ``my_test_cases.yaml``) to call via ``ufs_test.sh``:

.. code-block:: console

   ./ufs_test.sh -a epic -s -c -k -r -l my_test_cases.yaml