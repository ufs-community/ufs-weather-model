[![Read The Docs Status](https://readthedocs.org/projects/ufs-weather-model/badge/?badge=latest)](http://ufs-weather-model.readthedocs.org/)

# ufs-weather-model

This repository contains the UFS Weather Model source code.

# Where to find information

The [WM User's Guide](https://ufs-weather-model.readthedocs.org/) is hosted on ReadTheDocs.
The [WM wiki](https://github.com/ufs-community/ufs-weather-model/wiki) also has instructions for getting started.

# What files are what

The top level directory structure groups source code and input files as follow:

| File/directory            | Purpose |
| --------------            | ------- |
| ```LICENSE.md```          | A copy of the GNU Lesser General Public License, Version 3. |
| ```README.md```           | This file with basic pointers to more information. |
| ```AQM/```                | Contains air quality modeling component |
| ```CDEPS-interface/```    | Contains CDEPS data components, including data-atmosphere (DATM) and data-ocean (DOCN) model components|
| ```CICE-interface/```     | Contains CICE sea-ice model component including CICE6 and Icepack |
| ```CMakeModules/```       | Contains common cmake modulefiles used by Spack and CMake to find dependencies
| ```CMEPS-interface/```    | Contains CMEPS mediator |
| ```fire_behavior/```      | Contains the Community Fire Behavior Model component |
| ```GOCART/```             | Contains GOCART aerosol model component |
| ```HYCOM-interface/```    | Contains HYCOM ocean model component |
| ```LM4-driver/```         | Contains LM4 land component |
| ```MOM6-interface/```     | Contains MOM6 ocean model component |
| ```NOAHMP-interface/```   | Contains Noah-MP land model component |
| ```stochastic_physics/``` | Contains the stochastic physics source code |
| ```UFSATM/```                | Contains FV3 atmosphere model component including FV3 dynamical core, dynamics to physics driver, physics and IO. |
| ```WW3/```                | Contains community wave modeling framework WW3 |
| ```cmake/```              | Contains compile option files on various platforms. |
| ```modulefiles/```        | Contains module files on various platforms. |
| ```tests/```              | Regression and unit testing framework scripts. |
| ```tests-dev/```          | Developmental testing framework scripts. |
| ```build.sh```            | Script to build the model executable (also used by `tests/`) |

To use `build.sh` to build the coupled model with `FV3_GFS_v17_coupled_p8_ugwpv1` as the CCPP suite, run:
```
$> module use modulefiles
$> module load ufs_<machine>.<compiler>
$> CMAKE_FLAGS="-DAPP=S2S -D32BIT=ON -DHYDRO=ON -DCCPP_SUITES=FV3_GFS_v17_coupled_p8_ugwpv1" ./build.sh
```
where the machine is any Tier 1-4 platform listed in an existing modulefile and the compiler is `intel`, `intelllvm`, or `gnu`.

The build system is regularly tested on [Tier-1 platforms](
https://github.com/ufs-community/ufs-weather-model/wiki/Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers).
Configurations for other platforms that are available with UFS should be used with the understanding that they are not regularly
tested, and users will have to adapt the code to make it work on those platforms.

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
