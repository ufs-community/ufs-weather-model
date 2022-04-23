###############################################################################
### Doc
# This file turns ON (only ON, not OFF) a component in a valid application
# If a user wishes to add a new application, they should add an
# Application handle (APP_NAME) in the list of `VALID_APPS` in the
# top-level CMakeLists.txt
# Next, they can define an if-endif block with that APP_NAME and
# turn ON the components specific to this new application.
# Note, only the components required for the application should be turned ON.
# It is forbidden to turn OFF or change any other CMake option in this file.
###############################################################################

###############################################################################
### Configure Application Components
###############################################################################
if(APP MATCHES "^(ATM|ATMW|ATMAQ)$")
  set(FMS        ON  CACHE BOOL "Enable FMS"                 FORCE)
  set(FV3        ON  CACHE BOOL "Enable FV3"                 FORCE)
  set(STOCH_PHYS ON  CACHE BOOL "Enable Stochastic Physics"  FORCE)
  if(APP MATCHES "ATMW")
    set(WW3      ON  CACHE BOOL "Enable WAVEWATCH III"       FORCE)
    message("Configuring UFS app in Atmosphere with Waves mode")
  elseif(APP MATCHES "ATMAQ")
    set(AQM      ON  CACHE BOOL "Enable AQM"                 FORCE)
    message("Configuring UFS app in Atmosphere with Air Quality mode")
  else()
    message("Configuring UFS app in Atmosphere Only mode")
  endif()
endif()

if(APP MATCHES "^(NG-GODAS)$")
  set(CMEPS      ON  CACHE BOOL "Enable CMEPS"               FORCE)
  set(STOCH_PHYS ON  CACHE BOOL "Enable Stochastic Physics"  FORCE)
  set(FMS        ON  CACHE BOOL "Enable FMS"                 FORCE)
  set(MOM6       ON  CACHE BOOL "Enable MOM6"                FORCE)
  set(CICE6      ON  CACHE BOOL "Enable CICE6"               FORCE)
  set(CDEPS      ON  CACHE BOOL "Enable CDEPS"               FORCE)
  message("Configuring UFS app in (CDEPS) Data Atmosphere mode")
endif()

if(APP MATCHES "^(S2S|S2SA|S2SW|S2SWA)$")
  set(APP_MSG "Configuring UFS app in S2S")
  set(CMEPS      ON  CACHE BOOL "Enable CMEPS"               FORCE)
  set(FMS        ON  CACHE BOOL "Enable FMS"                 FORCE)
  set(FV3        ON  CACHE BOOL "Enable FV3"                 FORCE)
  set(MOM6       ON  CACHE BOOL "Enable MOM6"                FORCE)
  set(CICE6      ON  CACHE BOOL "Enable CICE6"               FORCE)
  set(STOCH_PHYS ON  CACHE BOOL "Enable Stochastic Physics"  FORCE)
  if(APP MATCHES "^S2SW")
    set(WW3      ON  CACHE BOOL "Enable WAVEWATCH III"       FORCE)
    string(CONCAT APP_MSG ${APP_MSG} " with Waves")
  endif()
  if(APP MATCHES "A$")
    set(UFS_GOCART ON  CACHE BOOL "Enable GOCART"            FORCE)
    if(WW3)
      string(CONCAT APP_MSG ${APP_MSG} " and Aerosols")
    else()
      string(CONCAT APP_MSG ${APP_MSG} " with Aerosols")
    endif()
  endif()
  message("${APP_MSG} mode")
endif()

if(APP MATCHES "^(HAFS|HAFSW|HAFS-ALL)$")
  set(CMEPS      ON  CACHE BOOL "Enable CMEPS"               FORCE)
  if(APP MATCHES "^(HAFS-ALL)$")
    set(CDEPS    ON  CACHE BOOL "Enable CDEPS"               FORCE)
    message("Configuring UFS app in HAFS with CDEPS mode")
  endif()
  set(FMS        ON  CACHE BOOL "Enable FMS"                 FORCE)
  set(FV3        ON  CACHE BOOL "Enable FV3"                 FORCE)
  set(STOCH_PHYS ON  CACHE BOOL "Enable Stochastic Physics"  FORCE)
  set(HYCOM      ON  CACHE BOOL "Enable HYCOM"               FORCE)
  if(APP MATCHES "^(HAFSW|HAFS-ALL)$")
    set(WW3      ON  CACHE BOOL "Enable WAVEWATCH III"       FORCE)
    message("Configuring UFS app in HAFS with Waves mode")
  endif()
endif()

if(APP MATCHES "^(ATMAERO)$")
  set(FMS        ON  CACHE BOOL "Enable FMS"                 FORCE)
  set(FV3        ON  CACHE BOOL "Enable FV3"                 FORCE)
  set(STOCH_PHYS ON  CACHE BOOL "Enable Stochastic Physics"  FORCE)
  set(UFS_GOCART ON  CACHE BOOL "Enable GOCART"              FORCE)
  message("Configuring UFS app in Atmosphere with Aerosols mode")
endif()
