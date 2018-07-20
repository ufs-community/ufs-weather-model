#-----------------------------------------------------------------------------
# Darwin uses Mach-O binaries, so we don't need to export our own functions.
if (NOT (${CMAKE_SYSTEM_NAME} MATCHES "Darwin"))
    set(LINK_EXPORT_DYNAMIC "-Wl,--export-dynamic")
    set(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} ${LINK_EXPORT_DYNAMIC}")
endif()
