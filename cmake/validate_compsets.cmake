###############################################################################
### Application Checks
###############################################################################
if((ATM OR ATMW) AND (S2S OR S2SW))
  message(FATAL_ERROR "ATM|ATMW=ON and S2S|S2SW=ON are incompatible, ABORT!")
endif()

if(DATM AND (S2S OR S2SW))
  message(FATAL_ERROR "DATM=ON and S2S|S2SW=ON are incompatible, ABORT!")
endif()

if(DATM AND (ATM OR ATMW))
  message(FATAL_ERROR "DATM=ON and ATM|ATMW=ON are incompatible, ABORT!")
endif()

###############################################################################
### Enable Components in Supported Applications
###############################################################################
if(ATM OR ATMW)
  set(FMS        ON  CACHE BOOL "Enable FMS"                 FORCE)
  set(FV3        ON  CACHE BOOL "Enable FV3"                 FORCE)
  set(STOCH_PHYS ON  CACHE BOOL "Enable Stochastic Physics"  FORCE)
  if(ATMW)
    set(WW3      ON  CACHE BOOL "Enable WAVEWATCH III"       FORCE)
    message("Configuring UFS app in Atmosphere with Waves mode")
  else()
    message("Configuring UFS app in Atmosphere Only mode")
  endif()
  return()
endif()

if(DATM_NEMS)
  set(CMEPS      ON  CACHE BOOL "Enable CMEPS"               FORCE)
  set(NEMSdatm   ON  CACHE BOOL "Enable NEMS DataAtm"        FORCE)
  set(FMS        ON  CACHE BOOL "Enable FMS"                 FORCE)
  set(MOM6       ON  CACHE BOOL "Enable MOM6"                FORCE)
  set(CICE6      ON  CACHE BOOL "Enable CICE6"               FORCE)
  message("Configuring UFS app in Data Atmosphere mode")
  return()
endif()

if(DATM)
  set(CMEPS      ON  CACHE BOOL "Enable CMEPS"               FORCE)
  set(CDEPS      ON  CACHE BOOL "Enable CDEPS"               FORCE)
  set(FMS        ON  CACHE BOOL "Enable FMS"                 FORCE)
  set(MOM6       ON  CACHE BOOL "Enable MOM6"                FORCE)
  set(CICE6      ON  CACHE BOOL "Enable CICE6"               FORCE)
  message("Configuring UFS app in Data Atmosphere mode")
  return()
endif()

if(S2S OR S2SW)
  set(CMEPS      ON  CACHE BOOL "Enable CMEPS"               FORCE)
  set(FMS        ON  CACHE BOOL "Enable FMS"                 FORCE)
  set(FV3        ON  CACHE BOOL "Enable FV3"                 FORCE)
  set(MOM6       ON  CACHE BOOL "Enable MOM6"                FORCE)
  set(CICE6      ON  CACHE BOOL "Enable CICE6"               FORCE)
  set(STOCH_PHYS ON  CACHE BOOL "Enable Stochastic Physics"  FORCE)
  set(32BIT      OFF CACHE BOOL "Disable 32BIT"              FORCE)
  if(S2SW)
    set(WW3      ON  CACHE BOOL "Enable WAVEWATCH III"       FORCE)
    message("Configuring UFS app in S2S with Waves mode")
  else()
    message("Configuring UFS app in S2S mode")
  endif()
  return()
endif()

###############################################################################
### The application you are building is unsupported, proceed with caution!
###############################################################################
message("Configuring UFS app in an unsupported configuration!")
