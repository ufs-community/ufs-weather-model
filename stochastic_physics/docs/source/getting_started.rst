Getting Started
==================================================
The stochastic physics currently only works with the UFS-atmosphere model

You should get the full system at https://github.com/ufs-community/ufs-weather-model, which will include the stochastic physics code.

In order to enable stochastic physics in a model run, you will need to turn it on via `namelist options <namelist_options.html>`_

If using the CIME workflow decribed at https://ufs-mrweather-app.readthedocs.io/en/latest/, please add do_sppt=T, etc. to user_nl_ufsatm in the case directory.

