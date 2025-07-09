Users may need to modify the baseline configuration file (``${UFS_WM}/tests-dev/baseline_setup.yaml``), which contains details on the location of staged input data, user-specific output directories, and batch job scheduling. The following variables are of particular importance:

* ``dprefix``: Set this value to an existing directory where the user has write permissions. 
* ``STMP``: Directory for baseline test output (typically ``${dprefix}/stmp4``)
* ``PTMP``: Directory for runtime files (typically ``${dprefix}/stmp2``)
