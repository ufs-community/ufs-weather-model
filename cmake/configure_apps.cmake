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
if(APP MATCHES "^(ATM|ATMW)$")
  set(FMS        ON  CACHE BOOL "Enable FMS"                 FORCE)
  set(FV3        ON  CACHE BOOL "Enable FV3"                 FORCE)
  set(STOCH_PHYS ON  CACHE BOOL "Enable Stochastic Physics"  FORCE)
  if(APP MATCHES "ATMW")
    set(WW3      ON  CACHE BOOL "Enable WAVEWATCH III"       FORCE)
    message("Configuring UFS app in Atmosphere with Waves mode")
  else()
    message("Configuring UFS app in Atmosphere Only mode")
  endif()
endif()

if(APP MATCHES "^(NG-GODAS|NG-GODAS-NEMSDATM)$")
  set(CMEPS      ON  CACHE BOOL "Enable CMEPS"               FORCE)
  set(FMS        ON  CACHE BOOL "Enable FMS"                 FORCE)
  set(MOM6       ON  CACHE BOOL "Enable MOM6"                FORCE)
  set(CICE6      ON  CACHE BOOL "Enable CICE6"               FORCE)
  if(APP MATCHES "NG-GODAS-NEMSDATM")
    set(NEMSdatm ON  CACHE BOOL "Enable NEMS DataAtm"        FORCE)
    message("Configuring UFS app in (NEMS) Data Atmosphere mode")
  elseif(APP MATCHES "NG-GODAS")
    set(CDEPS    ON  CACHE BOOL "Enable CDEPS"               FORCE)
    message("Configuring UFS app in (CDEPS) Data Atmosphere mode")
  endif()
endif()

if(APP MATCHES "^(S2S|S2SW)$")
  set(CMEPS      ON  CACHE BOOL "Enable CMEPS"               FORCE)
  set(FMS        ON  CACHE BOOL "Enable FMS"                 FORCE)
  set(FV3        ON  CACHE BOOL "Enable FV3"                 FORCE)
  set(MOM6       ON  CACHE BOOL "Enable MOM6"                FORCE)
  set(CICE6      ON  CACHE BOOL "Enable CICE6"               FORCE)
  set(STOCH_PHYS ON  CACHE BOOL "Enable Stochastic Physics"  FORCE)
  if(APP MATCHES "S2SW")
    set(WW3      ON  CACHE BOOL "Enable WAVEWATCH III"       FORCE)
    message("Configuring UFS app in S2S with Waves mode")
  else()
    message("Configuring UFS app in S2S mode")
  endif()
endif()
