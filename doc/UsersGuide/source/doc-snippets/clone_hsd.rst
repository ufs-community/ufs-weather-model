To start, recursively clone the repository:

.. code-block:: console

   git clone --recursive -b develop https://github.com/ufs-community/ufs-weather-model.git
   cd ufs-weather-model

After cloning, users may save (or "export") the path to the UFS WM in an environment variable:

.. code-block:: console

   export UFS_WM=$PWD

Although this step is optional, users may find it convenient when navigating between directories. This documentation will use ``${UFS_WM}`` to refer to the path to the ``ufs-weather-model`` directory, but users may choose to type out the full path instead. 
