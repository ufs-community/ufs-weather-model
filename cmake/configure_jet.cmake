message("")
message("Setting configuration for $ENV{CMAKE_Platform}")
message("")

get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("C       compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION} (${C_COMPILER_NAME})")
message("CXX     compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} (${CXX_COMPILER_NAME})")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} (${Fortran_COMPILER_NAME})")
message("")

option(DEBUG   "Enable DEBUG mode" OFF)
option(REPRO   "Enable REPRO mode" OFF)
option(VERBOSE "Enable VERBOSE mode" OFF)
option(32BIT   "Enable 32BIT (single precision arithmetic in dycore)" OFF)
option(OPENMP  "Enable OpenMP threading" ON)
option(AVX2    "Enable AVX2 instruction set" OFF)

option(INLINE_POST "Enable inline post" OFF)

include( cmake/${CMAKE_Fortran_COMPILER_ID}.cmake )

message("AVX2 is   ENABLED on Jet (multi-tagret executable)")
string (REPLACE "-xHOST" "-axSSE4.2,AVX,CORE-AVX2" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
string (REPLACE "-xHOST" "-axSSE4.2,AVX,CORE-AVX2" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

string(REPLACE "-i_dynamic" "-shared-intel"
       CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS
       "${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS}")

set(NEMSIO_INC $ENV{NEMSIO_INC})
set(NCEP_LIBS $ENV{NEMSIO_LIB} $ENV{BACIO_LIB4} $ENV{SP_LIBd} $ENV{W3EMC_LIBd} $ENV{W3NCO_LIBd})

set(ESMF_MOD ${ESMF_F90COMPILEPATHS})
set(ESMF_LIBS "${ESMF_F90ESMFLINKRPATHS} ${ESMF_F90ESMFLINKPATHS} ${ESMF_F90ESMFLINKLIBS}")

set(NETCDF_INC_DIR $ENV{NETCDF}/include)
set(NETCDF_LIBDIR $ENV{NETCDF}/lib)
set(NETCDF_LIBS -L$ENV{NETCDF}/lib -lnetcdff -lnetcdf)

message("")
