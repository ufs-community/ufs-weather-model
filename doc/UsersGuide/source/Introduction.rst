.. _Introduction:

*************************
Introduction
*************************

The Unified Forecast System (:term:`UFS`) Weather Model (:term:`WM`) is a prognostic model that can be
used for short- and medium-range research and operational forecasts, as exemplified by
its use in the operational Global Forecast System (GFS) of the National Oceanic and
Atmospheric Administration (NOAA). In addition to its use in NOAA's operational forecast systems, the UFS WM is the atmospheric model used in public UFS application releases, such as the Short-Range Weather (SRW) Application v2.2.0 release. These releases represent a snapshot of a continuously evolving system undergoing open
development. More information about the UFS can be found on the UFS Community Portal at https://ufscommunity.org/ and on the Earth Prediction Innovation Center (EPIC) website at https://epic.noaa.gov/get-code/ufs-weather-model/.

Key architectural elements of the UFS WM, along with links to external detailed documentation
for those elements, are listed below:

   * The `Finite-Volume Cubed-Sphere (FV3) dynamical core <https://noaa-emc.github.io/FV3_Dycore_ufs-v2.0.0/html/index.html>`__ is the computational part of an atmospheric model that solves the equations of fluid motion.

   * The `Flexible Modeling System <https://www.gfdl.noaa.gov/fms/>`__ (:term:`FMS`), is a software framework for supporting the efficient development, construction, execution, and scientific interpretation of atmospheric, oceanic, and climate system models. It is used for functions such as parallelization. 

   * `The Common-Community Physics Package <https://dtcenter.org/community-code/common-community-physics-package-ccpp>`__ (:term:`CCPP`), provides a framework and library of physics schemes, or :term:`parameterizations`, that support interoperable atmospheric physics. Atmospheric physics is a set of numerical methods approximating the effects of small-scale processes such as clouds, turbulence, radiation, and their interactions. 

   * `Stochastic physics <https://stochastic-physics.readthedocs.io/en/latest/>`__ schemes apply randomized perturbations to the physical tendencies, or physical parameters, of a model in order to compensate for model uncertainty. They include the Stochastic Kinetic Backscatter Scheme (SKEBS), the Stochastically Perturbed Parameterization Tendencies (SPPT) scheme, the perturbed boundary layer humidity (SHUM) scheme, the Stochastically Perturbed Parameterizations (SPP) scheme, Land Surface Model SPP (LSM-SPP), and the cellular automata method (:cite:t:`BengtssonEtAl2020`).

   * The libraries needed to build the system, which are bundled together via `spack-stack <https://spack-stack.readthedocs.io/en/latest/>`__ and include:
   
      * `National Centers for Environmental Prediction (NCEP) Libraries <https://github.com/NOAA-EMC/NCEPLIBS/wiki>`__
      * `Earth System Modeling Framework (ESMF) <https://earthsystemmodeling.org/>`__
      * `External libraries <https://github.com/NOAA-EMC/NCEPLIBS-external/wiki>`__

   * The build system used to compile the code and generate the executable.

   * The regression tests used to maintain software integrity as innovations are added.

.. COMMENT: Should NCEP, ESMF, and external libraries be grouped as part of HPC-Stack? Or is this a different set of libraries?

The UFS Weather Model is currently included in two UFS Application releases: The UFS Short-Range Weather (:term:`SRW`) Application v2.2.0 release (October 2023) and the UFS Medium Range Weather Application (:term:`MRW`) v1.1.0 release (October 2020). These UFS Apps also contain pre- and post-processing components, a comprehensive build system, and workflows for configuration and execution of the application. The SRW App v2.2.0 documentation and details can be found `here <https://ufs-srweather-app.readthedocs.io/en/release-public-v2.2.0/>`__. The MRW App v1.1.0 documentation and details can be found `here <https://ufs-mrweather-app.readthedocs.io/en/ufs-v1.1.0>`__.

The UFS WM code is portable and can be used with Linux or Mac operating systems and with Intel or GNU compilers. It has been tested on a variety of platforms widely used by atmospheric scientists, such as the NOAA Research Hera system, the National Center for Atmospheric Research (:term:`NCAR`) Derecho system, the National Science Foundation Stampede system, and Mac laptops.

.. note::

   At this time, the following aspects are unsupported: configurations in which a mediator is used to couple the atmospheric model to models of other earth domains (such as ocean, ice, and waves), horizontal resolutions other than the supported ones, different number or placement of vertical levels, the *cellular automata* stochastic scheme, and the use of different file formats for input and output.  It is expected that the UFS WM supported capabilities will be expanded in future releases.

.. COMMENT: Are coupled versions of the WM now supported? With 12 configurations it would seem that perhaps some are? 
.. COMMENT: Is the cellular automata stochastic scheme now supported?
.. COMMENT: Which horizontal/vertical levels & placements are supported? Just the default ones? 

Those wishing to contribute development to the UFS WM should become familiar with the procedures for running the model as a standalone component and for executing the regression tests described in the UFS WM GitHub `wiki <https://github.com/ufs-community/ufs-weather-model/wiki/Making-code-changes-in-the-UFS-weather-model-and-its-subcomponents>`__ to make sure no inadvertent changes to the results have been introduced during the development process.

Support for the UFS WM is provided through the `UFS Forum <https://forums.ufscommunity.org/forum/ufs-weather-model>`__ by the Developmental Testbed Center (DTC) and other groups involved in UFS development, such as NOAA's Environmental Modeling Center (:term:`EMC`), NOAA research laboratories (GFDL, NSSL, ESRL, and AOML), and :term:`NCAR`. UFS users and developers are encouraged not only to post questions, but also to help address questions posted by other members of the community.

This WM User's Guide is organized as follows:

   * :numref:`Chapter %s <CodeOverview>` (Code Overview) provides a description of the various code repositories from which source code is pulled and an overview of the directory structure.

   * :numref:`Chapter %s <BuildingAndRunning>` (Building and Running the WM) explains how to use the WM without an application.

   * :numref:`Chapter %s <InputsOutputs>` (Data: Input, Model Configuration, and Output Files) lists the model inputs and outputs and has a description of the key files.

   * :numref:`Chapter %s <ConfigParams>` (Configuration Parameters) lists the purpose and valid values for various configuration parameters.

   * :numref:`Chapter %s <AutomatedTesting>` (Automated Testing) describes UFS WM automated testing options.
   
   * :numref:`Chapter %s <FAQ>` (FAQ) lists frequently asked questions and answers.

Finally, :numref:`Chapters %s <Acronyms>` and :numref:`%s <Glossary>` contain a list of acronyms and a glossary, respectively.

.. bibliography:: references.bib
