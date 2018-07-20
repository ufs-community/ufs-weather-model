#------------------------------------------------------------------------------
# Define a macro to find and add sources
macro(add_sources)
    file(RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}/ccpp/src" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach(_src ${ARGN})
        if (_relPath)
            list (APPEND SOURCES "${_relPath}/${_src}")
        else()
            list (APPEND SOURCES "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SOURCES to parent directory
        set(SOURCES ${SOURCES} PARENT_SCOPE)
    endif()
endmacro()
