###############################################################################
# Disable all Components except applications
###############################################################################
set(FMS        OFF CACHE BOOL "Disable FMS"                FORCE)
set(CDEPS      OFF CACHE BOOL "Disable CDEPS"              FORCE)
set(FV3        OFF CACHE BOOL "Disable FV3"                FORCE)
set(MOM6       OFF CACHE BOOL "Disable MOM6"               FORCE)
set(CICE6      OFF CACHE BOOL "Disable CICE6"              FORCE)
set(WW3        OFF CACHE BOOL "Disable WaveWatch III"      FORCE)
set(STOCH_PHYS OFF CACHE BOOL "Disable Stochastic Physics" FORCE)

###############################################################################
### Checks
###############################################################################
if(DATM AND S2S)
  message(FATAL_ERROR "DATM=ON and S2S=ON are incompatible, ABORT!")
endif()

###############################################################################
### Enable Application Specific components
###############################################################################
if(DATM)
  set(FMS        ON  CACHE BOOL "Enable MOM6"                FORCE)
  set(FV3        OFF CACHE BOOL "Disable FV3"                FORCE)
  set(MOM6       ON  CACHE BOOL "Enable MOM6"                FORCE)
  set(CICE6      ON  CACHE BOOL "Enable CICE6"               FORCE)
  set(STOCH_PHYS OFF CACHE BOOL "Disable Stochastic Physics" FORCE)
  set(CDEPS      ON  CACHE BOOL "Enable CDEPS"               FORCE)
  message("Configuring UFS in Data Atmosphere application mode")
  return()
endif()

if(S2S)
  set(FMS        ON  CACHE BOOL "Enable MOM6"                FORCE)
  set(FV3        ON  CACHE BOOL "Enable FV3"                 FORCE)
  set(MOM6       ON  CACHE BOOL "Enable MOM6"                FORCE)
  set(CICE6      ON  CACHE BOOL "Enable CICE6"               FORCE)
  set(WW3        ON  CACHE BOOL "Enable WaveWatch III"       FORCE)
  set(STOCH_PHYS ON  CACHE BOOL "Enable Stochastic Physics"  FORCE)
  set(32BIT      OFF CACHE BOOL "Disable 32IT"               FORCE)
  message("Configuring UFS in S2S application mode")
  return()
endif()

