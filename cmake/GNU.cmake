# https://github.com/ufs-community/ufs-weather-model/issues/2159
if(MOVING_NEST)
  message(FATAL_ERROR "Option MOVING_NEST not compatible with ${CMAKE_Fortran_COMPILER_ID}}")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ggdb -fbacktrace -cpp -fcray-pointer -ffree-line-length-none -fno-range-check")

if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
endif()

if(NOT 32BIT)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
endif()

if(DEBUG)
    add_definitions(-DDEBUG)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check")
    # https://github.com/ufs-community/ufs-weather-model/issues/2155
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND ${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "arm64")
      set( CMAKE_Fortran_FLAGS_DEBUG     "${CMAKE_Fortran_FLAGS_DEBUG} -mcmodel=small" )
    else()
      set( CMAKE_Fortran_FLAGS_DEBUG     "${CMAKE_Fortran_FLAGS_DEBUG} -mcmodel=medium" )
    endif()
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
else()
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set(CMAKE_C_FLAGS_RELEASE "-O2")
endif()
