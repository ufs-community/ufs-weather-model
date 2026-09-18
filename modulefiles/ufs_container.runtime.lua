-- ufs_container.runtime.lua
-- Host-side runtime module for container-based UFS-WM regression tests.
--
-- This file is loaded on the HOST (not inside the container) before
-- apptainer/singularity is invoked, from within each Tier 1 platform's own
-- native job-card template ("module use .../modulefiles; module load
-- ufs_container.runtime") whenever a container image was requested. Which
-- container runtime module (if any) is needed depends on the platform the
-- job is actually running on -- identified via the MACHINE_ID environment
-- variable the calling job script already exports for its own use.

local host = os.getenv("MACHINE_ID") or ""

if (host == "hercules" or host == "orion") then
  -- Hercules / Orion: apptainer/singularity is only available via modules.
  load("singularity")
elseif (host == "derecho") then
  -- Derecho (NCAR): apptainer is only available via modules.
  load("apptainer")
elseif (host == "container") then
  -- community.sh's own generic platform (MACHINE_ID=container), targeting
  -- one fixed, site-configured cluster rather than a specific Tier 1 host
  -- above. Placeholder -- uncomment and adapt for your system, e.g.:
  -- load("apptainer")
  -- load("singularity")
  -- load("openmpi")
end
-- Gaea-c6, Ursa, NOAA Cloud: apptainer/singularity is already in PATH by
-- default -- nothing to load. Likewise when MACHINE_ID is otherwise
-- unrecognized or unset.
