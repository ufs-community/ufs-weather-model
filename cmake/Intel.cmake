
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -nowarn -sox -align array64byte")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qno-opt-dynamic-align")


if(32BIT)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i4 -real-size 32")
    add_definitions(-DOVERLOAD_R4)
    add_definitions(-DOVERLOAD_R8)
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i4 -real-size 64")
    if(NOT REPRO)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-prec-div -no-prec-sqrt")
    endif()
endif()

if(REPRO)
    add_definitions(-DREPRO)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qno-opt-dynamic-align")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qno-opt-dynamic-align")
elseif(DEBUG)
    add_definitions(-DDEBUG)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qno-opt-dynamic-align")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qno-opt-dynamic-align")
else()
    if(AVX2)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xCORE-AVX2 -qno-opt-dynamic-align")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xCORE-AVX2 -qno-opt-dynamic-align")
    elseif(SIMDMULTIARCH)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -axSSE4.2,AVX,CORE-AVX2,CORE-AVX512 -qno-opt-dynamic-align")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -axSSE4.2,AVX,CORE-AVX2,CORE-AVX512 -qno-opt-dynamic-align")
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qno-opt-dynamic-align")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qno-opt-dynamic-align")
    endif()
endif()

if(REPRO)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -debug minimal -fp-model consistent -qoverride-limits -g -traceback")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -debug minimal")
elseif(DEBUG)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -traceback -ftrapuv")
    if(DEBUG_LINKMPI)
      if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -link_mpi=dbg_mt")
      else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -link_mpi=dbg")
      endif()
    endif()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -ftrapuv -traceback")
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -debug minimal")
    set(FAST "-fast-transcendentals")
endif()

if(QUAD_PRECISION)
    add_definitions(-DENABLE_QUAD_PRECISION)
endif()

if(MULTI_GASES)
    add_definitions(-DMULTI_GASES)
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__IFC -sox -fp-model source")
