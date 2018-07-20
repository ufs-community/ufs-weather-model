#-----------------------------------------------------------------------------
# find_package(OpenMP) does not support Fortran and has other problems, too

macro(detect_openmp)

    # C
    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
        set (OpenMP_C_FLAGS "-qopenmp")
    elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
        set (OpenMP_C_FLAGS "-mp")
    elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
        set (OpenMP_C_FLAGS "-fopenmp")
    elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
        set (OpenMP_C_FLAGS "-fopenmp")
    else ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
        message (FATAL_ERROR "C compiler ${CMAKE_C_COMPILER_ID} not configured")
    endif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    message (STATUS "Detecting OpenMP flags for ${CMAKE_C_COMPILER_ID} C compiler: ${OpenMP_C_FLAGS}")

    # C++ not used
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set (OpenMP_CXX_FLAGS "-qopenmp")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
        set (OpenMP_CXX_FLAGS "-mp")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set (OpenMP_CXX_FLAGS "-fopenmp")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set (OpenMP_CXX_FLAGS "-fopenmp")
    else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        message (FATAL_ERROR "C compiler ${CMAKE_CXX_COMPILER_ID} not configured")
    endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    message (STATUS "Detecting OpenMP flags for ${CMAKE_CXX_COMPILER_ID} C++ compiler: ${OpenMP_CXX_FLAGS}")

    # Fortran
    if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
        set (OpenMP_Fortran_FLAGS "-qopenmp")
    elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
        set (OpenMP_Fortran_FLAGS "-mp")
    elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
        set (OpenMP_Fortran_FLAGS "-fopenmp")
    elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Clang")
        set (OpenMP_Fortran_FLAGS "-fopenmp")
    else ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
        message (FATAL_ERROR "Fortran compiler ${CMAKE_Fortran_COMPILER_ID} not configured")
    endif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    message (STATUS "Detecting OpenMP flags for ${CMAKE_Fortran_COMPILER_ID} Fortran compiler: ${OpenMP_Fortran_FLAGS}")

endmacro()
