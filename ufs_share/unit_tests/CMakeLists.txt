# ------------------
# Add dependencies
# ------------------
find_package(OpenMP REQUIRED COMPONENTS Fortran)

# ----------------------------------------------------------------------------
# Build the Unit Test Executable
# ----------------------------------------------------------------------------
add_executable(test_qc_exe test_qc_coupling.F90)

# Link it to the parent library. 
# Because qc_coupling linked ESMF publicly, test_qc_exe gets ESMF automatically!
target_link_libraries(test_qc_exe PRIVATE ufs_share)

# ----------------------------------------------------------------------------
# Stage the Test Data
# ----------------------------------------------------------------------------
# Copies test_fd.yaml from the source unit_tests/ folder directly into the 
# build/unit_tests/ folder so the executable can read it locally.
file(COPY test_qc.nml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# ----------------------------------------------------------------------------
# Register the Test
# ----------------------------------------------------------------------------
add_test(NAME ctest_qc_coupling COMMAND mpiexec -n 3 $<TARGET_FILE:test_qc_exe>)
# add_test(NAME ctest_qc_coupling COMMAND test_qc_exe)
set_tests_properties(ctest_qc_coupling PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")

