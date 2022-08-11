.. _Glossary:

*************************
Glossary
*************************

.. glossary::

   AQM
   Air Quality Model
      The AQM is a UFS Application that dynamically couples the Community Multiscale Air Quality (:term:`CMAQ`) model with the UFS Weather Model through the :term:`NUOPC` Layer to simulate temporal and spatial variations of atmospheric compositions (e.g., ozone and aerosol compositions). The CMAQ, treated as a column chemistry model, updates concentrations of chemical species (e.g., ozone and aerosol compositions) at each integration time step. The transport terms (e.g., advection and diffusion) of all chemical species are handled by the UFS Weather Model as tracers.

   ATM
      The Weather Model configuration that runs the standalone atmospheric model only. 

   CCPP
      The `Common Community Physics Package <https://dtcenter.org/community-code/common-community-physics-package-ccpp>`__ is a forecast-model agnostic, vetted collection of code containing atmospheric physical parameterizations and suites of parameterizations for use in Numerical Weather Prediction (:term:`NWP`) along with a framework that connects the physics to the host forecast model.

   CCPP-Framework
     The infrastructure that connects physics schemes with a host model; also refers to a software
     repository of the same name 

   CCPP-Physics
      The pool of CCPP-compliant physics schemes; also refers to a software repository of the same name

   chgres_cube
      The preprocessing software used to create initial and boundary condition files to "coldstart" the forecast model. It is part of :term:`UFS_UTILS`. 

   CICE
   Sea Ice Model
      `CICE <https://github.com/CICE-Consortium/CICE>`__ is a computationally efficient model for simulating the growth, melting, and movement of polar sea ice. It was designed as one component of coupled atmosphere-ocean-land-ice global climate models. CICE has several interacting components, including a model of ice dynamics, a transport model that describes advection of different state variables; and a vertical physics package called "Icepack". When coupled with other earth system model components, routines external to the CICE model prepare and execute data exchanges with an external "flux coupler".

      ..
         COMMENT: Clarify definition!!!

   CDEPS
   Community Data Models for Earth Predictive Systems
      The Community Data Models for Earth Predictive Systems repository (`CDEPS <https://github.com/ESCOMP/CDEPS>`__) contains a set of :term:`NUOPC`-compliant data components along with :term:`ESMF`-based "stream" code that enables new capabilities in selectively removing feedbacks in coupled model systems. The CDEPS data models perform the basic function of reading external data files, modifying those data, and then sending the data back to the :term:`CMEPS` mediator. The fields sent to the :term:`mediator` are the same as those that would be sent by an active component. This takes advantage of the fact that the mediator and other CMEPS-compliant model components have no fundamental knowledge of whether another component is fully active or just a data component.

      ..
         COMMENT: Clarify definition!!!

   CESM
   Community Earth System Model
      The `Community Earth System Model <https://www.cesm.ucar.edu/>`__ is a community climate model centered at the National Center for Atmospheric Research (:term:`NCAR`). 

   CMAQ
   Community Multiscale Air Quality Model
      The Community Multiscale Air Quality Model (`CMAQ <https://www.epa.gov/cmaq/cmaq-models-0>`__, pronounced cee-mak) is a numerical air quality model that predicts the concentration of airborne gases and particles and the deposition of these pollutants back to Earth's surface. The purpose of CMAQ is to provide fast, technically sound estimates of ozone, particulates, toxics, and acid deposition. CMAQ is an active open-source development project of the U.S. Environmental Protection Agency (EPA). Code is publicly availably at https://github.com/USEPA/CMAQ. 

   CMEPS
      The Community Mediator for Earth Prediction Systems (`CMEPS <https://github.com/ESCOMP/CMEPS>`__) is a :term:`NUOPC`-compliant :term:`mediator` used for coupling Earth system model components. It is currently being used in NCAR's Community Earth System Model (CESM) and NOAA's subseasonal-to-seasonal (S2S) coupled system.

   DATM
      ..
         COMMENT: Add definition!!!

   dycore
   dynamical core
      Global atmospheric model based on fluid dynamics principles, including Euler's equations of motion.

   EMC
   Environmental Modeling Center
      The `Environmental Modeling Center <https://www.emc.ncep.noaa.gov/emc_new.php>`__ is one of :term:`NCEP`'s nine centers and leads the :term:`National Weather Service`'s modeling efforts.

   ESMF
      `Earth System Modeling Framework <https://earthsystemmodeling.org/docs/release/latest/ESMF_usrdoc/>`__. The ESMF defines itself as "a suite of software tools for developing high-performance, multi-component Earth science modeling applications." It is a community-developed software infrastructure for building and coupling models. 

   FMS
     The Flexible Modeling System (FMS) is a software framework for supporting the efficient
     development, construction, execution, and scientific interpretation of atmospheric, 
     oceanic, and climate system models.

   FV3
   FV3 dycore
   FV3 dynamical core
      The Finite-Volume Cubed-Sphere :term:`dynamical core` (dycore). Developed at NOAA's `Geophysical 
      Fluid Dynamics Laboratory <https://www.gfdl.noaa.gov/>`__ (GFDL), it is a scalable and flexible dycore capable of both hydrostatic and non-hydrostatic atmospheric simulations. It is the dycore used in the UFS Weather Model.

   GOCART
      NASA's Goddard Chemistry Aerosol Radiation and Transport (GOCART) model simulates the distribution of major tropospheric aerosol types, including sulfate, dust, organic carbon (OC), black carbon (BC), and sea salt aerosols. The UFS Weather Model integrates a prognostic aerosol component using GOCART. The code is publicly available on GitHub at https://github.com/GEOS-ESM/GOCART.   
         
      ..
         COMMENT: Check definition!!! 
         
   HPC-Stack
      The `HPC-Stack <https://github.com/NOAA-EMC/hpc-stack>`__ is a repository that provides a unified, shell script-based build system for building the software stack required for numerical weather prediction (NWP) tools such as the `Unified Forecast System (UFS) <https://ufscommunity.org/>`__ and the `Joint Effort for Data assimilation Integration (JEDI) <https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/>`__ framework.

   HYCOM
   Hybrid Coordinate Ocean Model
      The HYbrid Coordinate Ocean Model (`HYCOM <https://www.hycom.org/>`__) was developed to address known shortcomings in the vertical coordinate scheme of the Miami Isopycnic-Coordinate Ocean Model (MICOM) developed by Rainer Bleck and colleagues. HYCOM is a primitive equation, general circulation model with vertical coordinates that remain isopycnic in the open, stratified ocean. However, the isopycnal vertical coordinates smoothly transition to z-coordinates in the weakly stratified upper-ocean mixed layer, to terrain-following sigma coordinates in shallow water regions, and back to z-level coordinates in very shallow water. The latter transition prevents layers from becoming too thin where the water is very shallow. See the `HYCOM User's Guide <https://www.hycom.org/attachments/063_hycom_users_guide.pdf>`__ for more information.
   
   Mediator
      A mediator, sometimes called a coupler, is a software component that includes code for representing component interactions. Typical operations include merging data fields, ensuring consistent treatment of coastlines, computing fluxes, and temporal averaging.

   MOM6
   Modular Ocean Model
      MOM6 is the latest generation of the Modular Ocean Model. It is numerical model code for simulating the ocean general circulation. MOM6 was originally developed by the `Geophysical Fluid Dynamics Laboratory <https://www.gfdl.noaa.gov/mom-ocean-model/>`__. Currently, `MOM6 code <https://github.com/mom-ocean/MOM6>`__ and an `extensive suite of test cases <https://github.com/NOAA-GFDL/MOM6-examples/wiki>`__ are available under an open-development software framework. Although there are many public forks of MOM6, the `NOAA EMC fork <https://github.com/NOAA-EMC/MOM6>`__ is used in the UFS Weather Model. 

   NWS
   National Weather Service
      The `National Weather Service <https://www.weather.gov/>`__ (NWS) is an agency of the United States government that is tasked with providing weather forecasts, warnings of hazardous weather, and other weather-related products to organizations and the public for the purposes of protection, safety, and general information. It is a part of the National Oceanic and Atmospheric Administration (NOAA) branch of the Department of Commerce.

   NWP
   Numerical Weather Prediction
      Numerical Weather Prediction (NWP) takes current observations of weather and processes them with computer models to forecast the future state of the weather. 

   NCAR
      The `National Center for Atmospheric Research <https://ncar.ucar.edu/>`__. 

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
      The NOAA Environmental Modeling System is a common modeling framework whose purpose is 
      to streamline components of operational modeling suites at :term:`NCEP`.

   NG-GODAS
      Next Generation-Global Ocean Data Assimilation System.

   .. ADD!!!

   NUOPC
   National Unified Operational Prediction Capability
      The `National Unified Operational Prediction Capability <https://earthsystemmodeling.org/nuopc/>`__ is a consortium of Navy, NOAA, and Air Force modelers and their research partners. It aims to advance the weather modeling systems used by meteorologists, mission planners, and decision makers. NUOPC partners are working toward a common model architecture --- a standard way of building models --- in order to make it easier to collaboratively build modeling systems.

   NUOPC Layer
      The :term:`NUOPC` Layer "defines conventions and a set of generic components for building coupled models using the Earth System Modeling Framework (:term:`ESMF`)." 
      NUOPC applications are built on four generic components: driver, model, :term:`mediator`, and connector. For more information, visit the `NUOPC website <https://earthsystemmodeling.org/nuopc/>`__.

   Parameterization
   Parameterizations
      Simplified functions that approximate the effects of small-scale processes (e.g., microphysics, gravity wave drag) that cannot be explicitly resolved by a model grid's representation of the earth. Common categories of parameterizations include radiation, surface layer, planetary boundary layer and vertical mixing, deep and shallow cumulus, and microphysics. Parameterizations can be grouped together into physics suites (such as the :term:`CCPP` physics suites), which are sets of parameterizations known to work well together. 

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

   WW3
   WWIII
   WaveWatch III
      WAVEWATCH III (WW3) is a community wave modeling framework that includes the latest scientific advancements in the field of wind-wave modeling and dynamics. The core of the framework consists of the WAVEWATCH III third-generation wave model (WAVE-height, WATer depth and Current Hindcasting), developed at NOAA/:term:`NCEP`. WAVEWATCH III differs from its predecessors in many important points such as governing equations, model structure, numerical methods and physical parameterizations. The model code is publicly available on GitHub at https://github.com/NOAA-EMC/WW3. 

   WM
   Weather Model
      A prognostic model that can be used for short- and medium-range research and 
      operational forecasts. It can be an atmosphere-only model or be an atmospheric
      model coupled with one or more additional components, such as a wave or ocean model.
      The UFS Weather Model repository is publicly available on `GitHub <https://github.com/ufs-community/ufs-weather-model>`__. 

