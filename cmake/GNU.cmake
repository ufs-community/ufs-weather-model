
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcray-pointer -ffree-line-length-none -fno-range-check -fbacktrace")

if(DEBUG)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O0 -ggdb")
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -fno-range-check")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
endif()


if(REPRO)
    if(APPLE)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -ggdb")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -ggdb")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
    endif()
else()
endif()


if(32BIT)
    add_definitions(-DOVERLOAD_R4)
    add_definitions(-DOVERLOAD_R8)
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
endif()


if(OPENMP)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
    add_definitions(-DOPENMP)
endif()

if(AVX2)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
endif()

if(QUAD_PRECISION)
    add_definitions(-DENABLE_QUAD_PRECISION)
endif()

if(MULTI_GASES)
    add_definitions(-DMULTI_GASES)
endif()
