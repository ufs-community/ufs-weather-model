-- ufs_container.runtime.lua
-- Host-side runtime module for container-based UFS-WM regression tests.
--
-- This file is loaded on the HOST (not inside the container) when the test
-- job card runs.  Use it to load any host modules that must be present before
-- launching apptainer/singularity, such as:
--   * the apptainer or singularity module itself (if not in PATH by default)
--   * a host MPI module that srun/mpirun needs
--   * any other system modules required at submission time
--
-- Platform-specific examples are shown below (all commented out).
-- Uncomment and adapt the lines appropriate for your system.

-- -- Orion / Hercules — load the singularity module.
-- load("singularity")

-- -- Gaea-c6 (GFDL) — apptainer is in PATH by default; no module needed.

-- -- Derecho (NCAR) — load apptainer, host GNU compilers, and host OpenMPI
-- load("apptainer")
-- load("gcc/14.3.0")
-- load("openmpi/5.0.9")

-- -- NOAA cloud — singularity is in PATH by default; no module needed

-- -- Generic placeholder: load host MPI if needed for srun
-- load("openmpi")
