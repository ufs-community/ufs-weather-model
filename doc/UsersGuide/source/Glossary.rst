.. _Glossary:

*************************
Glossary
*************************

.. glossary::

   CCPP
      Model agnostic, vetted, collection of codes containing atmospheric physical parameterizations
      and suites for use in NWP along with a framework that connects the physics to host models

   CCPP-Framework
     The infrastructure that connects physics schemes with a host model; also refers to a software
     repository of the same name 

   CCPP-Physics
      The pool of CCPP-compliant physics schemes; also refers to a software repository of the same name

   chgres_cube
      The preprocessing software used to create initial and boundary condition files to "coldstart" the forecast model. It is part of :term:`UFS_UTILS`. 

   dycore
   dynamical core
      Global atmospheric model based on fluid dynamics principles, including Euler's equations of motion.

   EMC
   Environmental Modeling Center
      The `Environmental Modeling Center <https://www.emc.ncep.noaa.gov/emc_new.php>`__ is one of :term:`NCEP`'s nine centers and leads the :term:`National Weather Service`'s modeling efforts.

   FMS
     The Flexible Modeling System (FMS) is a software framework for supporting the efficient
     development, construction, execution, and scientific interpretation of atmospheric, 
     oceanic, and climate system models.

   FV3
   FV3 dycore
   FV3 dynamical core
      The Finite-Volume Cubed-Sphere :term:`dynamical core` (dycore). Developed at NOAA's `Geophysical 
      Fluid Dynamics Laboratory <https://www.gfdl.noaa.gov/>`__ (GFDL), it is a scalable and flexible dycore capable of both hydrostatic and non-hydrostatic atmospheric simulations. It is the dycore used in the UFS Weather Model.

   HPC-Stack
      The `HPC-Stack <https://github.com/NOAA-EMC/hpc-stack>`__ is a repository that provides a unified, shell script-based build system for building the software stack required for numerical weather prediction (NWP) tools such as the `Unified Forecast System (UFS) <https://ufscommunity.org/>`__ and the `Joint Effort for Data assimilation Integration (JEDI) <https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/>`__ framework.

   NWS
   National Weather Service
      The `National Weather Service <https://www.weather.gov/>`__ (NWS) is an agency of the United States government that is tasked with providing weather forecasts, warnings of hazardous weather, and other weather-related products to organizations and the public for the purposes of protection, safety, and general information. It is a part of the National Oceanic and Atmospheric Administration (NOAA) branch of the Department of Commerce.

   NWP
   Numerical Weather Prediction
      Numerical Weather Prediction (NWP) takes current observations of weather and processes them with computer models to forecast the future state of the weather. 

   NCEP
   National Centers for Environmental Prediction
      National Centers for Environmental Prediction (NCEP) is a branch of the :term: `National Weather Service` and consists of nine centers, including the :term:`Environmental Modeling Center`. More information can be found at https://www.ncep.noaa.gov.

   NCEPLIBS
      The software libraries created and maintained by :term:`NCEP` that are required for running 
      :term:`chgres_cube`, the UFS Weather Model, and the :term:`UPP`. They are included in the `HPC-Stack <https://github.com/NOAA-EMC/hpc-stack>`__ and in `spack-stack <https://github.com/NOAA-EMC/spack-stack>`__. 

   NCEPLIBS-external
      A collection of third-party libraries required to build :term:`NCEPLIBS`, :term:`chgres_cube`, 
      the UFS Weather Model, and the :term:`UPP`. They are included in the :term:`HPC-Stack` and in :term:`spack-stack`.  

   NEMS
      The NOAA Environmental Modeling System - a software infrastructure that supports 
      NCEP/EMC’s forecast products.

   NUOPC
      The National Unified Operational Prediction Capability is a consortium of Navy, NOAA,
      and Air Force modelers and their research partners. It aims to advance the weather
      modeling systems used by meteorologists, mission planners, and decision makers. NUOPC
      partners are working toward a common model architecture - a standard way of building
      models - in order to make it easier to collaboratively build modeling systems.

   Parameterization or physics scheme
      The representation, in a dynamic model, of physical effects in terms of admittedly
      oversimplified parameters, rather than realistically requiring such effects to be 
      consequences of the dynamics of the system (AMS Glossary)

   Post-processor
      Software that enhances the value of the raw forecasts produced by the modeling application to make them more useful. At :term:`NCEP`, the :term:`UPP` (Unified Post Processor) software is used to convert data from spectral to gridded format, de-stagger grids, interpolate data vertically (e.g., to isobaric levels) and horizontally (to various predefined grids), and to compute derived variables. Some types of post-processors, such as statistical post-processors, use historical information of previous runs and observations to de-bias and calibrate its output.

   spack-stack
      The `spack-stack <https://github.com/NOAA-EMC/spack-stack>`__ is a collaborative effort between the NOAA Environmental Modeling Center (EMC), the UCAR Joint Center for Satellite Data Assimilation (JCSDA), and the Earth Prediction Innovation Center (EPIC). *spack-stack* is a repository that provides a Spack-based method for building the software stack required for numerical weather prediction (NWP) tools such as the `Unified Forecast System (UFS) <https://ufscommunity.org/>`__ and the `Joint Effort for Data assimilation Integration (JEDI) <https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/>`__ framework. *spack-stack* uses the Spack package manager along with custom Spack configuration files and Python scripts to simplify installation of the libraries required to run various applications. The *spack-stack* can be installed on a range of platforms and comes pre-configured for many systems. Users can install the necessary packages for a particular application and later add the missing packages for another application without having to rebuild the entire stack.

   Suite Definition File (SDF)
     An external file containing information about the 
     construction of a physics suite. It describes the schemes that are called, in which
     order they are called, whether they are subcycled, and whether they are assembled
     into groups to be called together

   Suite
      A collection of primary physics schemes and interstitial schemes that are known to work
      well together

   UFS
      A Unified Forecast System (UFS) is a community-based, coupled comprehensive Earth
      system modeling system. The UFS numerical applications span local to global domains
      and predictive time scales from sub-hourly analyses to seasonal predictions. It is
      designed to support the Weather Enterprise and to be the source system for NOAA's
      operational numerical weather prediction applications

   UFS_UTILS
      The :term:`UFS` Utilities repository contains a collection of pre-processing programs for use with the UFS Weather Model and UFS applications. These programs set up the model grid and create coldstart initial conditions. The code is publicly available on the `UFS_UTILS <https://github.com/ufs-community/UFS_UTILS>`__ Github repository. 

   UPP
   Unified Post Processor
      The `Unified Post Processor <https://dtcenter.org/community-code/unified-post-processor-upp>`__ is the :term:`post-processor` software developed at :term:`NCEP`. It is used operationally to 
      convert the raw output from a variety of :term:`NCEP`'s :term:`NWP` models, including the :term:`FV3 dycore`, to a more useful form.

   Weather Model
      A prognostic model that can be used for short- and medium-range research and 
      operational forecasts. It can be an atmosphere-only model or be an atmospheric
      model coupled with one or more additional components, such as a wave or ocean model.
      The UFS Weather Model repository is publicly available on `GitHub <https://github.com/ufs-community/ufs-weather-model>`__. 

