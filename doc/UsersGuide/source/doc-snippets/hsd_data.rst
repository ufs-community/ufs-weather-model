Data for the HSD cases is already staged on Tier-1 platforms at the ``INPUTDATA_ROOT*`` locations listed in `baseline_setup.yaml <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests-dev/baseline_setup.yaml>`_. However, users on any platform can download the data directly from the `HTF data bucket <https://registry.opendata.aws/noaa-ufs-htf-pds/>`_ using ``wget``. 

.. code-block:: console 

   wget https://noaa-ufs-htf-pds.s3.amazonaws.com/develop-20241115/HSD_INPUT_DATA/HSD_INPUT_DATA.tar.gz
   tar xvfz HSD_INPUT_DATA.tar.gz
