
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mpreprocess")

if(DEBUG)
    message("DEBUG is       ENABLED")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -Ktrap=fp -Mbounds -traceback")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -traceback")
else()
    message("DEBUG is       disabled (optimized build)")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
    set(FAST "-fast-transcendentals")
endif()


if(REPRO)
    message("REPRO is       ENABLED")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O1 -g -traceback")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
else()
    message("REPRO is       disabled")
endif()


if(32BIT)
    message("32BIT is       ENABLED")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i4 -r4")
    add_definitions(-DOVERLOAD_R4)
    add_definitions(-DOVERLOAD_R8)
else()
    message("32BIT is       disabled")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i4 -r8 -Mfprelaxed=div -Mfprelaxed=sqrt")
endif()


if(OPENMP)
    message("OPENMP is      ENABLED")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mp")
    add_definitions(-DOPENMP)
else()
    message("OPENMP is      disabled")
endif()


if(AVX2)
    message("AVX2 is        ENABLED")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -tp=haswell")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -tp=haswell")
else()
    message("AVX2 is        disabled")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -tp=x64")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -tp=x64")
endif()


if(INLINE_POST)
    message("INLINE_POST is ENABLED")
else()
    message("INLINE_POST is disabled")
endif()
