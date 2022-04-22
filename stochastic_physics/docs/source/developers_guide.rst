Developer's guide
=================

Code is housed on github at https://github.com/noaa-psd/stochastic_physics.  Please read more about the Development process at https://github.com/ufs-community/ufs/wiki/Developing-with-Gitflow.

Please make a fork and checkout the entire ufs-community weather model at https://github.com/ufs-community/ufs-weather-model and point to your fork of the stochastic_physics submodule.

Standalone testing
""""""""""""""""""
If you intend to make modifications to the stochastic physics source code, there is a simplified program that exercises the random pattern generator without needing to run the entire model.  Please see README.standalone in the stochastic_physics directory.

Full model tests
""""""""""""""""
The code updates are not expected to change existing results, so the full model regression tests need to be run.  All of the tests must pass, although only a sub-set of tests are needed to consider adding changes to the stochastic_physics repository: fv3_control, fv3_stochy, fv3_ccpp_control, and fv3_ccpp_stochy.  If the results are expected to change, then there needs to be scientific evidience that the change in results are what is expected.  
