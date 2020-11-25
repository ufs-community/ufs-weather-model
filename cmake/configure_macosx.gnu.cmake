set(INLINE_POST     ON  CACHE BOOL "Enable inline post" FORCE)

# OpenMP broken for clang compiler
if(${CMAKE_C_COMPILER_ID} MATCHES "^(Clang|AppleClang)$")
  set(OPENMP        OFF CACHE BOOL "Enable OpenMP threading" FORCE)
endif()
