# Implementation Plan: CECE UFS Integration

## Overview

Integrate CECE into the UFS weather model by wiring it into the CMake build system, registering it in the NUOPC driver, adding the ATME application configuration, and providing a runtime configuration template. Each task builds incrementally — build system first, then driver, then app config, then runtime config, then final wiring verification.

## Tasks

- [x] 1. Add CECE CMake option and build integration to top-level CMakeLists.txt
  - [x] 1.1 Declare the `CECE` cache variable and add status message
    - Add `set(CECE OFF CACHE BOOL "Enable CECE")` alongside existing component options in `CMakeLists.txt`
    - Add `ATME` to the `VALID_APPS` list
    - Add `message("CECE ............. ${CECE}")` in the component status messages block
    - _Requirements: 1.1, 1.3, 4.4_
  - [x] 1.2 Add conditional `add_subdirectory(CECE)` block
    - Add a new `### CECE ###` section after the AQM section with `if(CECE) add_subdirectory(CECE) endif()`
    - CECE builds directly from its submodule — no wrapper directory needed
    - _Requirements: 1.2, 2.1_
  - [x] 1.3 Wire CECE into the `ufs` library target
    - In the UFS Library section, add the CECE block following the AQM pattern:
      - `add_dependencies(ufs cece)`
      - `list(APPEND _ufs_defs_private FRONT_CECE=cece_cap_mod)`
      - `list(APPEND _ufs_libs_public cece)`
    - _Requirements: 2.2, 2.3, 2.4, 7.3_

- [x] 2. Checkpoint - Verify CMake changes
  - Ensure all CMake changes are syntactically correct, ask the user if questions arise.

- [x] 3. Add ATME application configuration in configure_apps.cmake
  - [x] 3.1 Add the ATME application block
    - Add a new `if(APP MATCHES "^(ATME)$")` block in `cmake/configure_apps.cmake`
    - Enable `FMS`, `FV3`, `STOCH_PHYS`, and `CECE` with `ON CACHE BOOL ... FORCE`
    - Add message: `"Configuring UFS app in Atmosphere with Emissions mode"`
    - Place the block after the existing ATM-family blocks (after the `ATMF` / `ATMMPAS` handling)
    - _Requirements: 4.1, 4.2, 4.3_

- [x] 4. Add CECE driver registration in UFSDriver.F90
  - [x] 4.1 Add the CECE import block
    - Add `#ifdef FRONT_CECE` / `use FRONT_CECE, only: CECE_SS => SetServices` / `#endif` in the module import section of `driver/UFSDriver.F90`, after the AQM import block
    - _Requirements: 3.1, 7.1, 7.2_
  - [x] 4.2 Add the CECE registration block in SetModelServices
    - Add a `#ifdef FRONT_CECE` guarded block in the component registration loop, after the AQM block
    - Match on `trim(model) == "cece"`
    - Include threading bail-out for `ompNumThreads > 1` (matching AQM/DATM pattern)
    - Call `NUOPC_DriverAddComp(driver, trim(prefix), CECE_SS, petList=petList, comp=comp, rc=rc)`
    - Set `found_comp = .true.`
    - _Requirements: 3.2, 3.3, 3.4_

- [x] 5. Checkpoint - Verify driver changes compile
  - Ensure all Fortran preprocessor changes are syntactically correct, ask the user if questions arise.

- [x] 6. Create runtime configuration template
  - [x] 6.1 Create `tests/parm/ufs.configure.atme.IN`
    - Create the file with ESMF, EARTH, ATM, and CECE sections
    - Set `EARTH_component_list: ATM CECE`
    - Define `CECE_model:`, `CECE_petlist_bounds:`, `CECE_omp_num_threads:` with `@[...]` template placeholders
    - Define `CECE_attributes::` section with `Verbosity = 0` and `ConfigFile = @[cece_config_file]`
    - Define `runSeq::` with ATM first, then `ATM -> CECE`, then `CECE`, then `CECE -> ATM` inside a coupling interval
    - _Requirements: 5.1, 5.2, 5.3, 5.4, 5.5_

- [x] 7. Final checkpoint - Verify all integration files
  - Ensure all files are consistent and correctly cross-reference each other, ask the user if questions arise.

## Notes

- No changes to the CECE cap (`CECE/src/cece_cap.F90`) are required — it already implements IPDv01 with the expected interface
- CECE builds directly from `add_subdirectory(CECE)` using its own CMakeLists.txt which handles Kokkos, yaml-cpp, and ESMF dependencies via FetchContent
- Threading is disabled initially (matching AQM pattern) because CECE manages its own parallelism via Kokkos
- The design explicitly notes that property-based testing does not apply to this build/driver integration task
- Each task references specific requirements for traceability
