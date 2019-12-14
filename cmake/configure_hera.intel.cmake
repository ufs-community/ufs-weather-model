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
option(AVX2    "Enable AVX2 instruction set" ON)

option(INLINE_POST "Enable inline post" ON)

include( cmake/${CMAKE_Fortran_COMPILER_ID}.cmake )

set(NEMSIO_INC $ENV{NEMSIO_INC})
set(POST_INC $ENV{POST_INC})
set(NCEP_LIBS $ENV{POST_LIB} $ENV{NEMSIO_LIB} $ENV{G2_LIB4} $ENV{G2TMPL_LIB} $ENV{BACIO_LIB4} $ENV{SP_LIBd} $ENV{W3EMC_LIBd} $ENV{W3NCO_LIBd} $ENV{CRTM_LIB} $ENV{PNG_LIB} $ENV{JASPER_LIB} $ENV{Z_LIB})

set(ESMF_MOD ${ESMF_F90COMPILEPATHS})
set(ESMF_LIBS "${ESMF_F90ESMFLINKRPATHS} ${ESMF_F90ESMFLINKPATHS} ${ESMF_F90ESMFLINKLIBS}")

set(NETCDF_INC_DIR $ENV{NETCDF}/include)
set(NETCDF_LIBDIR $ENV{NETCDF}/lib)
set(NETCDF_LIBS -L$ENV{NETCDF}/lib -lnetcdff -lnetcdf)

message("")
