# Requirements Document

## Introduction

This document specifies the requirements for integrating CECE (Chemistry/Emissions Computation Engine) into the UFS (Unified Forecast System) weather model as a NUOPC component. CECE is a chemistry and emissions component that already exists as a git submodule with a NUOPC cap (`cece_cap_mod`) implementing the IPDv01 initialization protocol. The integration follows the established UFS pattern used by existing components such as AQM, CDEPS, NOAHMP, and GOCART — covering the CMake build system, driver registration, application configuration, and runtime configuration.

## Glossary

- **UFS_Build_System**: The top-level CMake build system defined in `CMakeLists.txt` that declares component options, finds dependencies, builds component subdirectories, and links them into the `ufs` library with preprocessor definitions.
- **UFS_Driver**: The NUOPC driver component defined in `driver/UFSDriver.F90` that imports component SetServices routines via preprocessor guards, registers components in `SetModelServices`, and reads runtime configuration from `ufs.configure`.
- **App_Configurator**: The CMake application configuration module at `cmake/configure_apps.cmake` that maps application names (e.g., ATM, ATMAQ, S2S) to sets of enabled components.
- **CECE_Interface**: A CMake build directory (`CECE-interface/`) or equivalent build integration that compiles the CECE submodule sources and produces a linkable library target for the UFS build system.
- **CECE_Cap**: The existing NUOPC Model cap at `CECE/src/cece_cap.F90` (module `cece_cap_mod`) that implements IPDv01 with `CECE_SetServices`, `CECE_InitializeAdvertise`, `CECE_InitializeRealize`, `CECE_Run`, and `CECE_Finalize`.
- **Runtime_Config**: The `ufs.configure` template files (e.g., `tests/parm/ufs.configure.*.IN`) that define `EARTH_component_list`, per-component model/petlist/thread settings, attributes, and the `runSeq::` execution sequence.
- **NUOPC**: The National Unified Operational Prediction Capability interoperability layer built on ESMF that standardizes component coupling via SetServices, Advertise, Realize, Run, and Finalize phases.
- **IPDv01**: Initialize Phase Definition version 01, the NUOPC initialization protocol used by CECE and DATM, consisting of Advertise (phase 1) and Realize (phase 3) steps.
- **FRONT_CECE**: The C preprocessor macro that, when defined, enables conditional compilation of CECE import and registration code in `UFSDriver.F90`, set to the Fortran module name providing `SetServices`.
- **CECE_Application**: A UFS application configuration (e.g., ATME or an extension of ATMAQ) that enables the CECE component alongside required companion components.
- **Kokkos**: A C++ performance portability library used by CECE for parallel execution on CPUs and GPUs.
- **TIDE**: The data ingestor library used by CECE for reading input data streams, built as a subdirectory of the CECE source tree.

## Requirements

### Requirement 1: CMake Component Option Declaration

**User Story:** As a UFS developer, I want CECE declared as a CMake build option in the top-level `CMakeLists.txt`, so that CECE can be conditionally enabled or disabled like other UFS components.

#### Acceptance Criteria

1. THE UFS_Build_System SHALL declare a `CECE` CMake cache variable of type BOOL with a default value of `OFF`.
2. WHEN the `CECE` option is `ON`, THE UFS_Build_System SHALL include the CECE build directory via `add_subdirectory`.
3. THE UFS_Build_System SHALL print the CECE option status in the component status messages alongside existing components (AQM, GOCART, MOM6, etc.).

### Requirement 2: CECE Build Integration

**User Story:** As a UFS developer, I want CECE compiled and linked into the UFS executable through the CMake build system, so that the CECE NUOPC cap and its C++ core are available at link time.

#### Acceptance Criteria

1. WHEN the `CECE` option is `ON`, THE UFS_Build_System SHALL build the CECE library target from the CECE submodule sources.
2. WHEN the `CECE` option is `ON`, THE UFS_Build_System SHALL add a dependency from the `ufs` library target to the CECE library target.
3. WHEN the `CECE` option is `ON`, THE UFS_Build_System SHALL define the preprocessor macro `FRONT_CECE=cece_cap_mod` for the `ufs` library target.
4. WHEN the `CECE` option is `ON`, THE UFS_Build_System SHALL link the CECE library target as a public dependency of the `ufs` library target.
5. THE CECE_Interface SHALL resolve CECE's C++ dependencies (Kokkos, yaml-cpp) either through CECE's own FetchContent declarations or through pre-installed system packages.
6. THE CECE_Interface SHALL compile the Fortran NUOPC cap (`cece_cap.F90`) and the TIDE data ingestor alongside the C++ core sources.

### Requirement 3: Driver Registration

**User Story:** As a UFS developer, I want CECE registered in the UFS driver so that the NUOPC framework can initialize, run, and finalize CECE as part of a coupled simulation.

#### Acceptance Criteria

1. WHEN `FRONT_CECE` is defined, THE UFS_Driver SHALL import `CECE_SS` (aliased from `CECE_SetServices`) from the `FRONT_CECE` module.
2. WHEN `FRONT_CECE` is defined and the runtime configuration specifies `model` equal to `"cece"`, THE UFS_Driver SHALL call `NUOPC_DriverAddComp` to register CECE with the driver using `CECE_SS` and the configured `petList`.
3. WHEN `FRONT_CECE` is defined and the runtime configuration specifies `model` equal to `"cece"` and `ompNumThreads` is greater than 1, THE UFS_Driver SHALL reject the configuration with an error message indicating that ESMF-aware threading is not implemented for CECE.
4. WHEN `FRONT_CECE` is not defined and the runtime configuration specifies `model` equal to `"cece"`, THE UFS_Driver SHALL report an error indicating that no component named `"cece"` was found.

### Requirement 4: Application Configuration

**User Story:** As a UFS developer, I want at least one UFS application configuration that enables CECE, so that CECE can be built and tested as part of a recognized UFS application.

#### Acceptance Criteria

1. THE App_Configurator SHALL define at least one application name that enables the `CECE` component.
2. WHEN the CECE-enabled application is selected, THE App_Configurator SHALL set `CECE` to `ON`.
3. WHEN the CECE-enabled application is selected, THE App_Configurator SHALL also enable all companion components required by CECE (at minimum: `FMS`, `FV3`, `STOCH_PHYS` for an atmosphere-coupled configuration).
4. THE App_Configurator SHALL add the CECE-enabled application name to the `VALID_APPS` list in the top-level `CMakeLists.txt`.

### Requirement 5: Runtime Configuration Template

**User Story:** As a UFS developer, I want a `ufs.configure` template file for CECE-enabled runs, so that the NUOPC driver can read CECE's component settings and execute it in the correct run sequence.

#### Acceptance Criteria

1. THE Runtime_Config SHALL include `CECE` in the `EARTH_component_list` alongside `ATM` and any other required components.
2. THE Runtime_Config SHALL define `CECE_model:`, `CECE_petlist_bounds:`, and `CECE_omp_num_threads:` configuration entries with template placeholders.
3. THE Runtime_Config SHALL define a `CECE_attributes::` section for component-specific attributes.
4. THE Runtime_Config SHALL include CECE in the `runSeq::` block, defining the execution order and data exchange pattern between ATM and CECE.
5. WHEN CECE is coupled with ATM, THE Runtime_Config SHALL define the run sequence such that ATM executes first, ATM transfers data to CECE, CECE executes, and CECE transfers data back to ATM within each coupling timestep.

### Requirement 6: CECE Cap NUOPC Compliance

**User Story:** As a UFS developer, I want the CECE NUOPC cap to be compatible with the UFS driver's component lifecycle, so that CECE initializes, runs, and finalizes correctly within the coupled system.

#### Acceptance Criteria

1. THE CECE_Cap SHALL export a `CECE_SetServices` subroutine that accepts an `ESMF_GridComp` and an integer return code, matching the signature expected by `NUOPC_DriverAddComp`.
2. THE CECE_Cap SHALL implement the IPDv01 initialization protocol with an Advertise phase (IPDv01p1) and a Realize phase (IPDv01p3).
3. THE CECE_Cap SHALL advertise import fields for meteorological data and export fields for chemical species during the Advertise phase.
4. THE CECE_Cap SHALL create and allocate ESMF fields on a valid ESMF Grid or Mesh during the Realize phase.
5. THE CECE_Cap SHALL implement a Run phase registered under the `label_Advance` specialization that reads the ESMF clock, processes import fields, executes the CECE computation, and populates export fields.
6. THE CECE_Cap SHALL implement a Finalize phase registered under the `label_Finalize` specialization that releases all allocated resources.
7. IF the CECE core initialization fails during any phase, THEN THE CECE_Cap SHALL return a non-success ESMF return code to the driver.

### Requirement 7: Preprocessor Guard Consistency

**User Story:** As a UFS developer, I want all CECE-related code paths guarded by the `FRONT_CECE` preprocessor macro, so that the UFS builds correctly both with and without CECE enabled.

#### Acceptance Criteria

1. THE UFS_Driver SHALL compile without errors WHEN `FRONT_CECE` is not defined.
2. THE UFS_Driver SHALL compile without errors WHEN `FRONT_CECE` is defined and the CECE library target is available.
3. THE UFS_Build_System SHALL define `FRONT_CECE` only WHEN the `CECE` CMake option is `ON`.

### Requirement 8: CECE Build Isolation

**User Story:** As a UFS developer, I want CECE's build to not interfere with other UFS components, so that enabling CECE does not break existing application builds.

#### Acceptance Criteria

1. WHEN the `CECE` option is `OFF`, THE UFS_Build_System SHALL not attempt to find or build any CECE dependencies (Kokkos, yaml-cpp).
2. WHEN the `CECE` option is `ON`, THE CECE_Interface SHALL not modify global CMake compiler flags or definitions that affect other UFS component targets.
3. WHEN the `CECE` option is `ON`, THE CECE_Interface SHALL place Fortran module files in a CECE-specific module directory to avoid name collisions with other components.
