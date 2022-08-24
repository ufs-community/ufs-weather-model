.. _Introduction:

*************************
Introduction
*************************

The Unified Forecast System (:term:`UFS`) Weather Model (:term:`WM`) is a prognostic model that can be
used for short- and medium-range research and operational forecasts, as exemplified by
its use in the operational Global Forecast System (GFS) of the National Oceanic and
Atmospheric Administration (NOAA). The UFS WM v2.0 is the latest public release of this
software and represents a snapshot of a continuously evolving system undergoing open
development. More information about the UFS can be found in its portal at https://ufscommunity.org/.

Key architectural elements of the UFS WM, along with links to external detailed documentation
for those elements, are listed below:

   * The `Finite-Volume Cubed-Sphere (FV3) dynamical core <https://noaa-emc.github.io/FV3_Dycore_ufs-v2.0.0/html/index.html>`__ is the computational part of an atmospheric model that solves the equations of fluid motion..

   * The `Flexible Modeling System <https://www.gfdl.noaa.gov/fms/>`__ (:term:`FMS`), is a software framework for supporting the efficient development, construction, execution, and scientific interpretation of atmospheric, oceanic, and climate system models. It is used for functions such as parallelization. 

   * `The Common-Community Physics Package <https://dtcenter.org/community-code/common-community-physics-package-ccpp>`__ (:term:`CCPP`), provides a framework and library of physics schemes, or :term:`parameterizations`, that support interoperable atmospheric physics. Atmospheric physics is a set of numerical methods approximating the effects of small-scale processes such as clouds, turbulence, radiation, and their interactions. 

   * `Stochastic physics <https://stochastic-physics.readthedocs.io/en/ufs-v2.0.0/>`__ schemes apply randomized perturbations to the physical tendencies, or physical parameters, of a model in order to compensate for model uncertainty. They include Stochastic Kinetic Backscatter Scheme (SKEBS), the Stochastically Perturbed Parameterization Tendencies (SPPT) scheme, the perturbed boundary layer humidity (SHUM) scheme, and the cellular automata method.

   ..
      COMMENT: Add SPP and SPP LSM? Remove cellular automata method?

   .. COMMENT: Remove?: 
      - `The NOAA Environmental Modeling System <https://noaa-emc.github.io/NEMS_doc_ufs-v2.0.0/index.html>`__ (:term:`NEMS`) is a software infrastructure that supports NCEP/EMC's forecast products. In particular, the model driver is used to organize calls to various WM components hierarchically.

   * The libraries needed to build the system, such as:
   
      * `National Centers for Environmental Prediction (NCEP) Libraries <https://github.com/NOAA-EMC/NCEPLIBS/wiki>`__
      * `Earth System Modeling Framework (ESMF) <https://www.earthsystemcog.org/projects/esmf/>`__
      * `External libraries <https://github.com/NOAA-EMC/NCEPLIBS-external/wiki>`__

   * The build system used to compile the code and generate the executable.

   * The regression tests used to maintain software integrity as innovations are added.

The UFS Weather Model is currently included in two UFS Application releases.  These UFS Apps also contain pre- and post-processing components, a comprehensive build system, and workflows for configuration and execution of the application.

The UFS WM is included as part of the UFS Short Range Weather App, and details can be found `here <https://ufs-srweather-app.readthedocs.io/en/release-public-v2/>`__.

The UFS WM v1.1 and v1.0 is included as part of the UFS Medium Range Weather App, and details can be found `here <https://ufs-mrweather-app.readthedocs.io/en/ufs-v1.1.0>`__.

The UFS WM v2.0 code is portable and can be used with Linux or Mac operating systems and with Intel or GNU compilers. It has been tested in a variety of platforms widely used by atmospheric scientists, such as the NOAA Research Hera system, the National Center for Atmospheric Research (NCAR) Cheyenne system, the National Science Foundation Stampede system, and Mac laptops.

.. note::

   At this time, the following aspects are unsupported: configurations in which a mediator is used to couple the atmospheric model to models of other earth domains (such as ocean, ice, and waves), horizontal resolutions other than the supported ones, different number or placement of vertical levels, the *cellular automata* stochastic scheme, and the use of different file formats for input and output.  It is expected that the UFS WM supported capabilities will be expanded in future releases.

Those wishing to contribute development to the UFS WM should become familiar with the procedures for running the model as a standalone component and for executing the regression tests described in the UFS WM GitHub `wiki <https://github.com/ufs-community/ufs-weather-model/wiki/Making-code-changes-in-the-UFS-weather-model-and-its-subcomponents>`_ to make sure no inadvertent changes to the results have been introduced during the development process.

Support for the UFS WM is provided through the `UFS Forum <https://forums.ufscommunity.org/forum/ufs-weather-model>`_ by the Developmental Testbed Center (DTC) and other groups involved in UFS development, such as NOAA’s Environmental Modeling Center (EMC), NOAA research laboratories (GFDL, NSSL, ESRL, and AOML), and NCAR. UFS users and developers are encouraged not only to post questions, but also to help address questions posted by other members of the community.

This WM User's Guide is organized as follows:

   * :numref:`Chapter %s <CodeOverview>` (Code Overview) provides a description of the various code repositories from which source code is pulled and an overview of the directory structure.

   * :numref:`Chapter %s <BuildingAndRunning>` (Building and Running the WM) explains how to use the WM without an application.

   * :numref:`Chapter %s <InputsOutputs>` (Inputs and Outputs) lists the model inputs and outputs and has a description of the key files.

   * :numref:`Chapter %s <FAQ>` (FAQ) lists frequently asked questions and answers.

Finally, :numref:`Chapters %s <Acronyms>` and :numref:`%s <Glossary>` contain a list of acronyms and a glossary, respectively.

.. This is how you cite a reference :cite:`Bernardet2018`.

