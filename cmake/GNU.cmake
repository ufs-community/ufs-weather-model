set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ggdb -fbacktrace -cpp -fcray-pointer -ffree-line-length-none -fno-range-check")

if(NOT 32BIT)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
endif()

if(DEBUG)
    add_definitions(-DDEBUG)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check -mcmodel=medium")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
else()
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set(CMAKE_C_FLAGS_RELEASE "-O2")
endif()
