.. role:: raw-html(raw)
    :format: html

.. _aquaplanet:

*********************
AquaPlanet Test Case
*********************

The AquaPlanet test case is an idealized atmosphere-only forecast configuration designed to study atmospheric dynamics in a simplified Earth-like setting where all land is replaced with ocean. This configuration removes the complexity of land-surface interactions, topography, and regional variations, allowing researchers to focus on fundamental atmospheric processes such as tropical convection, Hadley circulation, jet stream dynamics, and the global energy budget.

The test case runs at C48 resolution with the ``FV3_GFS_v17_p8_ugwpv1`` physics suite. Initial conditions are created by modifying standard GFS data to represent an aquaplanet configuration: sea surface temperatures (SST) follow an idealized latitudinal profile, all land is converted to ocean, topography is set to zero, and sea ice is removed. The atmospheric initial state is configured with idealized vertical profiles of temperature, humidity, and winds appropriate for an aquaplanet simulation.

A key feature of this configuration is the 90-day spin-up period required to allow the model to adjust from its initial state to a balanced aquaplanet climate. After spin-up, the simulation can be run for extended periods to study seasonal variations, climate statistics, and atmospheric phenomena in the absence of land-surface influences.

============================
Obtaining Data for HSD Cases
============================

.. include:: ./doc-snippets/hsd_data.rst

.. _setup-aquaplanet:

==========================================
Setting Up the AquaPlanet Experiment
==========================================

This section describes how to set up the AquaPlanet experiment from scratch. The process involves modifying input data files, running a 90-day spin-up simulation, and then running experiments.

Clone Required Repositories
----------------------------

First, clone the UFS Weather Model and UFS_UTILS repositories:

.. code-block:: console

   git clone https://github.com/ufs-community/ufs-weather-model.git
   cd ufs-weather-model/
   git submodule update --init --recursive

Then clone the UFS_UTILS repository:

.. code-block:: console

   git clone https://github.com/ufs-community/UFS_UTILS.git
   cd UFS_UTILS/
   git submodule update --init --recursive

Finally, clone the aquaplanet tools repository:

.. code-block:: console

   git clone https://github.com/RatkoVasic-NOAA/Aquaplanet

Build UFS Weather Model
------------------------

Navigate to the UFS Weather Model test directory and configure the regression test system to compile the model:

.. code-block:: console

   cd ufs-weather-model/tests/

Edit ``rt.conf`` to include only the following two lines:

.. code-block:: console

   COMPILE | atm_dyn32 | intel | -DAPP=ATM -DCCPP_SUITES=FV3_GFS_v17_p8_ugwpv1 -D32BIT=ON | | fv3 |
   RUN | control_c48 | | baseline |

Execute the regression test to compile the model and create a baseline run directory:

.. code-block:: console

   ./rt.sh -a <account> -k

Replace ``<account>`` with your project code (e.g., ``epic``). Save the location of the ``control_c48_intel`` run directory for later use.

Build UFS_UTILS
---------------

Compile the UFS_UTILS tools:

.. code-block:: console

   cd UFS_UTILS/
   module purge
   module use $PWD/modulefiles
   module load build.<platform>.intelllvm
   ./build_all.sh

Replace ``<platform>`` with your system name (e.g., ``ursa``).

Link the fix directories:

.. code-block:: console

   cd fix/
   ./link_fixdirs.sh emc <platform>

Since orography files need to be edited, remove the link and copy the files:

.. code-block:: console

   rm orog
   mkdir orog
   cp /scratch3/NCEPDEV/global/role.glopara/fix/orog/20240917/*.nc orog/.
   cp /scratch3/NCEPDEV/global/role.glopara/fix/orog/20240917/*.dat orog/.
   cp -r /scratch3/NCEPDEV/global/role.glopara/fix/orog/20240917/C48/ orog/.