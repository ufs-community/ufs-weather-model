.. _Introduction:

*************************
Introduction
*************************

The Unified Forecast System (:term:`UFS`) :term:`Weather Model` (WM) is a prognostic model that can be
used for short- and medium-range research and operational forecasts, as exemplified by
its use in the operational Global Forecast System (GFS) of the National Oceanic and
Atmospheric Administration (NOAA). The UFS WM v1.1 is the latest public release of this
software and represents a snapshot of a continuously evolving system undergoing open
development. More information about the UFS can be found in its portal at https://ufscommunity.org/.

Key architectural elements of the UFS WM, along with links to external detailed documentation
for those elements, are listed below:

- `The Finite-Volume Cubed-Sphere (FV3) dynamical core <https://noaa-emc.github.io/FV3_Dycore_ufs-v1.1.0/html/index.html>`_.

- `The Flexible Modeling System <https://www.gfdl.noaa.gov/fms/>`_ (:term:`FMS`), a software infrastructure used for functions such as
  parallelization.

- `The Common-Community Physics Package <https://dtcenter.org/community-code/common-community-physics-package-ccpp>`_ (:term:`CCPP`), a library of
  physical parameterizations and the framework to use it with the model. :term:`Parameterization or physics scheme` is defined here.

- `The stochastic physics capability <https://stochastic-physics.readthedocs.io/en/release-ufs-v1.1.0/>`_, including the Stochastic Kinetic Backscatter Scheme (SKEBS),
  the Stochastically Perturbed Parameterization Tendencies (SPPT) scheme, the perturbed boundary
  layer humidity (SHUM) scheme, and the cellular automata method.

- `The NOAA Environmental Modeling System <https://noaa-emc.github.io/NEMS_doc_ufs-v1.1.0/index.html>`_ (:term:`NEMS`) model driver used to create the main program.

- The libraries needed to build the system, such as:
    - `National Centers for Environmental Prediction (NCEP) Libraries <https://github.com/NOAA-EMC/NCEPLIBS/wiki>`_
    - `Earth System Modeling Framework (ESMF) <https://www.earthsystemcog.org/projects/esmf/>`_
    - `External libraries <https://github.com/NOAA-EMC/NCEPLIBS-external/wiki>`_

- The build system used to compile the code and generate the executable.

- The regression tests used to maintain software integrity as innovations are added.

For the UFS WM v1.1 release, the following aspects are supported:

- Global configuration with resolutions of C96 (~100 km), C192 (~50 km), C384 (25 km), and C768 (~13 km)

- Sixty-four vertical levels at predetermined locations.

- Four physics suites (:term:`suite`), corresponding to GFS v15.2 (operational at the time of the release) and
  GFS v16beta (October 2019 version, in preparation for operational implementation in 2021). Variants
  with and without prediction of Sea Surface Temperature (SST) are included.

- Ability to run with or without SKEBS, SPPT, and SHUM.

- Ability to initialize from GFS files in Gridded Binary v2 (GRIB2), NEMS Input/Output (NEMSIO), or
  Network Common Data Form (netCDF) format for past dates, starting January 1, 2018, when the 
  preprocessing utility chgres_cube is employed.  Dates before that may work, but are not guaranteed.

- Output files in Network Common Data Form (NetCDF) format.

The GFS_v15p2 physics suite uses the following physical parameterizations: the
Simplified Arakawa Schubert shallow and deep convective schemes, the Geophysical
Fluid Dynamics Laboratory (GFDL) microphysics scheme, the Noah Land Surface Model (LSM),
the Rapid Radiative Transfer Model for Global Circulation Models (RRTMG) radiation scheme,
the hybrid eddy-diffusivity mass-flux (EDMF) planetary boundary layer (PBL) scheme based on the Smagorinsky K theory,
an orographic gravity wave drag (GWD) parameterization, and the Near SST (NSST) ocean scheme to predict SST.
In the GFS_v16beta suite, a moist TKE-based EDMF scheme replaces the K-based one and a non-stationary GWD parameterization is added.
The GFS_v15p2_no_nsst and the GFS_v16beta_no_nsst suites use a simple ocean scheme instead of the NSST scheme.
This simple ocean scheme keeps the SST constant throughout the forecast and is recommended for use when the initial
conditions do not contain all fields needed to initialize the NSST scheme.


Even when using physics suite GFS_v15p2, the UFS WM v1.1 differs from the operational GFS v15.2 in a few ways. First, the public release code
reflects the state of development as of the fall of 2019,
and therefore the parameterizations contain innovations beyond what is in GFSv15.2 operations.
For example, the GFDL microphysics distributed for use in GFS v15.2 and GFS v16beta
is the same scheme and contains development beyond what was transitioned to operations
for GFS v15 in June 2019. Second, the public release code uses the CCPP as the
interface for calling physics, while in operations the Interoperable Physics Driver
(IPD) is used. NOAA is currently working toward phasing out the IPD from UFS applications.
Validation tests demonstrated that CCPP and IPD give bit-for-bit identical results
when the same physics is employed and selected performance flags are excluded at
compilation time. When performance compiler flags employed in operational production are used, runs with
CCPP and IPD for the same physics suite yield differences comparable to running
the model in different computational platforms. Finally, the operational GFS
runs in NOAA Central Operations computational platforms. When users run the model
in different platforms, the results will differ.

It should also be noted that further changes are expected to the GFS v16 suite before it is implemented in operations in 2021.

The UFS WM v1 code is portable and can be used with Linux and Mac operating systems with Intel and GNU compilers. It has been tested in a variety of platforms widely used by atmospheric scientists, such as the NOAA research Hera system, the National Center for Atmospheric Research (NCAR) Cheyenne system, the National Science Foundation Stampede system, and Mac laptops.

.. note::

   At this time, the following aspects are unsupported:  standalone regional domains, configurations in which a mediator is used to couple the atmospheric model to models of other earth domains (such as ocean, ice, and waves), horizontal resolutions other than the supported ones, different number or placement of vertical levels, physics suites other than GFS v15.2 and GFS v16beta, the *cellular automata* stochastic scheme, initialization from sources other than GFS, the use of different file formats for input and output, and the use of the model in different computational platforms. It is expected that the UFS WM supported capabilities will be expanded in future releases.

It should be noted that the UFS WM is a component of the UFS Medium-Range (MR) Weather Application (App), which also contains pre- and post-processing components, a comprehensive build system, and workflows for configuration and execution of the application. At this time, the UFS WM is only supported to the general community for use as part of the UFS MR Weather App. However, those wishing to contribute development to the UFS WM should become familiar with the procedures for running the model as a standalone component and for executing the regression tests described in the UFS WM GitHub `wiki <https://github.com/ufs-community/ufs-weather-model/wiki/Making-code-changes-in-the-UFS-weather-model-and-its-subcomponents>`_ to make sure no inadvertent changes to the results have been introduced during the development process.

Support for the UFS WM is provided through the `UFS Forum <https://forums.ufscommunity.org/forum/ufs-weather-model>`_ by the Developmental Testbed Center (DTC) and other groups involved in UFS development, such as NOAA’s Environmental Modeling Center (EMC), NOAA research laboratories (GFDL, NSSL, ESRL, and AOML), and NCAR. UFS users and developers are encouraged not only to post questions, but also to help address questions posted by other members of the community.

This WM User’s Guide is organized as follows:

- :numref:`Chapter %s <CodeOverview>` (Code Overview) provides a description of the various
  code repositories from which source code is pulled and an overview of the directory structure.

- :numref:`Chapter %s <BuildingAndRunning>` (Building and Running the WM) explains how to use the WM without an application.

- :numref:`Chapter %s <InputsOutputs>` (Inputs and Outputs) lists the model inputs and outputs
  and has a description of the key files.

- :numref:`Chapter %s <SDFandNamelistExamplePractices>` (SDF and namelist samples and best practices)
  contains a description of the :term:`Suite Definition File (SDF)` and namelists needed to configure the model
  for running with the GFS v15.2 and GFS v16beta physics suites.

- :numref:`Chapter %s <FAQ>` (FAQ) lists frequently asked questions and answers.

Finally, :numref:`Chapters %s <Acronyms>` and :numref:`%s <Glossary>` contain a list of acronyms and a glossary, respectively.

.. This is how you cite a reference :cite:`Bernardet2018`.

.. bibliography:: references.bib
