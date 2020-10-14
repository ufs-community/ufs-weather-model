[![Read The Docs Status](https://readthedocs.org/projects/ufs-weather-model/badge/?badge=latest)](http://ufs-weather-model.readthedocs.io/)

# ufs-weather-model

This is the UFS weather model code.

# Where to find information

Start at the [ufs-weather-model wiki](https://github.com/ufs-community/ufs-weather-model/wiki) which has quick start instructions.

[User's reference guide](http://ufs-weather-model.readthedocs.io/) is hosted on read the docs.

# What files are what

The top level directory structure groups source code and input files as follow:

| File/directory            | Purpose |
| --------------            | ------- |
| ```LICENSE.md```          | A copy of the GNU Lesser General Public License, Version 3. |
| ```README.md```           | This file with basic pointers to more information. |
| ```FMS/```                | Contains Flexible Modeling System source code. |
| ```NEMS/```               | Contains NOAA Environmental Modeling System source code and nems compset runi scripts. |
| ```CMEPS-interface/```    | Contains CMEPS mediator |
| ```FV3/```                | Contains FV3 atmosphere model component including FV3 dynamical core, dynamics to physics driver, physics and IO. |
| ```WW3/```                | Contains community wave modeling framework WW3. |
| ```MOM6-interface/```     | Contains MOM6 ocean model component |
| ```CICE-interface/```     | Contains CICE sea-ice model component including CICE6 and Icepack |
| ```stochastic physics/``` | Contains the stochastic physics source code. |
| ```cmake/```              | Contains compile option files on various platforms. |
| ```modulefiles/```        | Contains module files on various platforms.|

# Disclaimer

The United States Department of Commerce (DOC) GitHub project code is provided
on an "as is" basis and the user assumes responsibility for its use. DOC has
relinquished control of the information and no longer has responsibility to
protect the integrity, confidentiality, or availability of the information. Any
claims against the Department of Commerce stemming from the use of its GitHub
project will be governed by all applicable Federal law. Any reference to
specific commercial products, processes, or services by service mark,
trademark, manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce. The
Department of Commerce seal and logo, or the seal and logo of a DOC bureau,
shall not be used in any manner to imply endorsement of any commercial product
or activity by DOC or the United States Government.
