.. _Introduction:
  
*************************
Introduction
*************************

The Unified Forecast System (:term:`UFS`) :term:`Weather Model` (WM) is a prognostic model that can be
used for short- and medium-range research and operational forecasts, as exemplified by
its use in the operational Global Forecast System (GFS) of the National Oceanic and 
Atmospheric Administration (NOAA). The UFS WM v1.0 is the first public release of this
software and represents a snapshot of a continuously evolving system undergoing open
development. More information about the UFS can be found in its portal at https://ufscommunity.org/. 

Key architectural elements of the UFS WM, along with links to external detailed documentation
for those elements, are listed below:

- The Finite-Volume Cubed-Sphere (FV3) dynamical core (https://noaa-emc.github.io/FV3_Dycore/html/index.html).

- The Flexible Modeling System (:term:`FMS`), a software infrastructure used for functions such as
  parallelization (https://www.gfdl.noaa.gov/fms/).

- The Common-Community Physics Package (:term:`CCPP`) library of 
  physical parameterizations (:term:`Parameterization or physics scheme`) and the 
  framework to use it with the model
  (https://dtcenter.org/community-code/common-community-physics-package-ccpp).

- The stochastic physics capability, including the Stochastic Kinetic Backscatter Scheme (SKEBS),
  the Stochastically Perturbed Parameterization Tendencies (SPPT) scheme, the perturbed boundary
  layer humidity (SHUM) scheme, and the cellular automata method 
  (https://stochastic-physics.readthedocs.io/en/ufs_public_release/).

- The NOAA Environmental Modeling System (:term:`NEMS`) model driver used to create the main program
  (https://docs.google.com/document/d/1-kFhPBf7GBTUd5SaB5D_3OUGX_93pWKP21QgTh4y6ok/edit#heading=h.dah4y9bxn10l).

- The libraries needed to build the system, such as: 
    - National Centers for Environmental Prediction (NCEP) Libraries 
    - Earth System Modeling Framework (ESMF; https://www.earthsystemcog.org/projects/esmf/)
    - System libraries

- The build system used to compile the code and generate the executable.

- The regression tests used to maintain software integrity as innovations are added.

For the UFS WM v1.0 release, the following aspects are supported:

- Global configuration with resolutions of C96 (~100 km), C192 (~50 km), C384 (25 km), and C768 (~13 km)

- Sixty-four vertical levels at predetermined locations.

- Two physics suites (:term:`suite`), corresponding to GFS v15.2 (operational at the time of the release) and
  GFS v16beta (October 2019 version, in preparation for operational implementation in 2021).
 
- Ability to run with or without SKEBS, SPPT, and SHUM.
 
- Ability to initialize from GFS files in Gridded Binary v2 (GRIB2) format for past dates, 
  starting January 1, 2018, when the preprocessing utility chgres is employed. Dates before
  that may work, but are not guaranteed.
 
- Output files in Network Common Data Form (NetCDF) format.

The GFS v15.2 physics suite uses the following physical parameterization: the Simplified Arakawa Schubert shallow and deep convective schemes, the Geophysical Fluid Dynamics Laboratory (GFDL) microphysics scheme, the Noah Land Surface Model (LSM), the Rapid Radiative Transfer Model for Global Circulation Models (RRTMG) radiation scheme, the hybrid eddy-diffusivity mass-flux (EDMF) planetary boundary layer (PBL) scheme based on the Smagorinsky K theory, and an orographic gravity wave drag (GWD) parameterization. In the GFS v16beta suite, a moist TKE-based EDMF scheme replaces the K-based one and a non-stationary GWD parameterization is added. It should be noted that the public release code reflects the state of development as of the fall of 2019, and therefore the parameterizations contains innovations beyond what is currently in operations. In other words, the GFDL microphysics distributed for use in GFS v15.2 and GFS v16beta is the same scheme and contains development beyond what was transitioned to operations for GFS v15 in June 2019. It should also be noted that further changes are expected to the GFS v16 suite before it is implemented in operations in 2021.

The UFS WM v1 code is portable and can be used with Linux and Mac operating systems with Intel and GNU compilers. It has been tested in a variety of platforms widely used by atmospheric scientists, such as the NOAA research Hera system, the National Center for Atmospheric Research (NCAR) Cheyenne system, the National Science Foundation Stampede system, and Mac laptops.

.. note::

   At this time, the following aspects are unsupported:  standalone regional domains, configurations in which a mediator is used to couple the atmospheric model to models of other earth domains (such as ocean, ice, and waves), horizontal resolutions other than the supported ones, different number or placement of vertical levels, physics suites other than GFS v15.2 and GFS v16beta the *cellular automata* stochastic scheme, initialization from sources other than GFS, the use of different file formats for input and output, and the use of the model in different computational platforms. It is expected that the UFS WM supported capabilities will be expanded in future releases.

It should be noted that the UFS WM is a component of the UFS Medium-Range Weather Application, which also contains pre- and post-processing components, a comprehensive build system, and workflows for configuration and execution of the application. At this time, the UFS WM is only supported to the general community for use as part of the UFS Medium-Range Weather Application. However, those wishing to contribute development to the UFS WM should become familiar with the procedures for running the model as a standalone component and for executing the regression tests to make sure no inadvertent changes to the results have been introduced during the development process.

Support for the UFS WM is provided through the UFS Forum by the Developmental Testbed Center (DTC) and other groups involved in UFS development, such as NOAA’s Environmental Modeling Center (EMC), NOAA research laboratories (GFDL, NSSL, ESRL, and AOML), and NCAR. UFS users and developers are encouraged not only to post questions, but also to help address questions posted by other members of the community. 

This WM User’s Guide is organized as follows:

- :numref:`Chapter %s <CodeOverview>` (Code Overview) provides a description of the various
  code repositories from which source code is pulled and an overview of the directory structure. 

- :numref:`Chapter %s <InputsOutputs>` (Inputs and Outputs) lists the model inputs and outputs
  and has a description of the key files.

- :numref:`Chapter %s <SDFandNamelistExamplePractices>` (SDF and namelist samples and best practices)
  contains a description of the :term:`Suite Definition File (SDF)` and namelists needed to configure the model
  for running with the GFS v15.2 and GFS v16beta physics suites. 

- :numref:`Chapter %s <FAQforModelConfiguration>` (FAQ on model configuration) contains information on
  miscellaneous topics pertaining to using the model in configurations that differ from the default. 

The next three chapters:

- :numref:`Chapter %s <ContributingDevelopment>` (Contributing development)
- :numref:`Chapter %s <CompilingCodeWithoutApp>` (Compiling the WM code without an application)
- :numref:`Chapter %s <RegTests>` (Regression tests for development)

go beyond the capabilities supported in the public release to cover code management for conducting
development and proposing contributions back to the authoritative code repositories. It should be noted that the regression tests described here are mandatory for committing code back to the ufs-weather-model authoritative code repository. These regressions tests differ from those distributed with the workflows for UFS applications, which are intended for application users and developers to assess the quality of their installations and the impact of their code changes. Finally,
:numref:`Chapters %s <Acronyms>` and :numref:`%s <Glossary>` contain a list of acronyms and a glossary, respectively.

.. This is how you cite a reference :cite:`Bernardet2018`.

.. bibliography:: references.bib
