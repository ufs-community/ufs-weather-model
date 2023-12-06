#!/bin/bash

# Do not set ACCNR here or it will break the per-system defaults in rt.sh.

case $(hostname -f) in

  clogin01.cactus.wcoss2.ncep.noaa.gov)	  MACHINE_ID=wcoss2 ;; ### cactus
  clogin02.cactus.wcoss2.ncep.noaa.gov)	  MACHINE_ID=wcoss2 ;; ### cactus
  clogin03.cactus.wcoss2.ncep.noaa.gov)	  MACHINE_ID=wcoss2 ;; ### cactus
  clogin04.cactus.wcoss2.ncep.noaa.gov)	  MACHINE_ID=wcoss2 ;; ### cactus
  clogin05.cactus.wcoss2.ncep.noaa.gov)	  MACHINE_ID=wcoss2 ;; ### cactus
  clogin06.cactus.wcoss2.ncep.noaa.gov)	  MACHINE_ID=wcoss2 ;; ### cactus
  clogin07.cactus.wcoss2.ncep.noaa.gov)	  MACHINE_ID=wcoss2 ;; ### cactus
  clogin08.cactus.wcoss2.ncep.noaa.gov)	  MACHINE_ID=wcoss2 ;; ### cactus
  clogin09.cactus.wcoss2.ncep.noaa.gov)	  MACHINE_ID=wcoss2 ;; ### cactus

  dlogin01.dogwood.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2 ;; ### dogwood
  dlogin02.dogwood.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2 ;; ### dogwood
  dlogin03.dogwood.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2 ;; ### dogwood
  dlogin04.dogwood.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2 ;; ### dogwood
  dlogin05.dogwood.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2 ;; ### dogwood
  dlogin06.dogwood.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2 ;; ### dogwood
  dlogin07.dogwood.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2 ;; ### dogwood
  dlogin08.dogwood.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2 ;; ### dogwood
  dlogin09.dogwood.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2 ;; ### dogwood

  alogin01.acorn.wcoss2.ncep.noaa.gov)  MACHINE_ID=acorn ;; ### acorn
  alogin02.acorn.wcoss2.ncep.noaa.gov)  MACHINE_ID=acorn ;; ### acorn
  alogin03.acorn.wcoss2.ncep.noaa.gov)  MACHINE_ID=acorn ;; ### acorn

  gaea10)                  MACHINE_ID=gaea ;; ### gaea10
  gaea11)                  MACHINE_ID=gaea ;; ### gaea11
  gaea12)                  MACHINE_ID=gaea ;; ### gaea12
  gaea13)                  MACHINE_ID=gaea ;; ### gaea13
  gaea14)                  MACHINE_ID=gaea ;; ### gaea14
  gaea15)                  MACHINE_ID=gaea ;; ### gaea15
  gaea16)                  MACHINE_ID=gaea ;; ### gaea16
  gaea10.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea10
  gaea11.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea11
  gaea12.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea12
  gaea13.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea13
  gaea14.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea14
  gaea15.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea15
  gaea16.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea16

  gaea51.ncrc.gov)         MACHINE_ID=gaea-c5 ;; ### gaea51
  gaea52.ncrc.gov)         MACHINE_ID=gaea-c5 ;; ### gaea52
  gaea53.ncrc.gov)         MACHINE_ID=gaea-c5 ;; ### gaea53
  gaea54.ncrc.gov)         MACHINE_ID=gaea-c5 ;; ### gaea54
  gaea55.ncrc.gov)         MACHINE_ID=gaea-c5 ;; ### gaea55
  gaea56.ncrc.gov)         MACHINE_ID=gaea-c5 ;; ### gaea56
  gaea57.ncrc.gov)         MACHINE_ID=gaea-c5 ;; ### gaea57
  gaea58.ncrc.gov)         MACHINE_ID=gaea-c5 ;; ### gaea58

  hfe01)                   MACHINE_ID=hera ;; ### hera01
  hfe02)                   MACHINE_ID=hera ;; ### hera02
  hfe03)                   MACHINE_ID=hera ;; ### hera03
  hfe04)                   MACHINE_ID=hera ;; ### hera04
  hfe05)                   MACHINE_ID=hera ;; ### hera05
  hfe06)                   MACHINE_ID=hera ;; ### hera06
  hfe07)                   MACHINE_ID=hera ;; ### hera07
  hfe08)                   MACHINE_ID=hera ;; ### hera08
  hfe09)                   MACHINE_ID=hera ;; ### hera09
  hfe10)                   MACHINE_ID=hera ;; ### hera10
  hfe11)                   MACHINE_ID=hera ;; ### hera11
  hfe12)                   MACHINE_ID=hera ;; ### hera12
  hecflow01)               MACHINE_ID=hera ;; ### heraecflow01

  s4-submit.ssec.wisc.edu) MACHINE_ID=s4 ;; ### s4

  fe1)                     MACHINE_ID=jet ;; ### jet01
  fe2)                     MACHINE_ID=jet ;; ### jet02
  fe3)                     MACHINE_ID=jet ;; ### jet03
  fe4)                     MACHINE_ID=jet ;; ### jet04
  fe5)                     MACHINE_ID=jet ;; ### jet05
  fe6)                     MACHINE_ID=jet ;; ### jet06
  fe7)                     MACHINE_ID=jet ;; ### jet07
  fe8)                     MACHINE_ID=jet ;; ### jet08
  tfe1)                    MACHINE_ID=jet ;; ### jet09
  tfe2)                    MACHINE_ID=jet ;; ### jet10

  Orion-login-1.HPC.MsState.Edu) MACHINE_ID=orion ;; ### orion1
  Orion-login-2.HPC.MsState.Edu) MACHINE_ID=orion ;; ### orion2
  Orion-login-3.HPC.MsState.Edu) MACHINE_ID=orion ;; ### orion3
  Orion-login-4.HPC.MsState.Edu) MACHINE_ID=orion ;; ### orion4

  derecho1.hsn.de.hpc.ucar.edu) MACHINE_ID=derecho ;; ### derecho1
  derecho2.hsn.de.hpc.ucar.edu) MACHINE_ID=derecho ;; ### derecho2
  derecho3.hsn.de.hpc.ucar.edu) MACHINE_ID=derecho ;; ### derecho3
  derecho4.hsn.de.hpc.ucar.edu) MACHINE_ID=derecho ;; ### derecho4
  derecho5.hsn.de.hpc.ucar.edu) MACHINE_ID=derecho ;; ### derecho5
  derecho6.hsn.de.hpc.ucar.edu) MACHINE_ID=derecho ;; ### derecho6
  derecho7.hsn.de.hpc.ucar.edu) MACHINE_ID=derecho ;; ### derecho7
  derecho8.hsn.de.hpc.ucar.edu) MACHINE_ID=derecho ;; ### derecho8

  Hercules-login-1.HPC.MsState.Edu) MACHINE_ID=hercules;; ### hercules1 
  Hercules-login-2.HPC.MsState.Edu) MACHINE_ID=hercules;; ### hercules2 
  Hercules-login-3.HPC.MsState.Edu) MACHINE_ID=hercules;; ### hercules3
  Hercules-login-4.HPC.MsState.Edu) MACHINE_ID=hercules;; ### hercules4

  login1.stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede1
  login2.stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede2
  login3.stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede3
  login4.stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede4
  
  
  login01.expanse.sdsc.edu) MACHINE_ID=expanse ;; ### expanse1
  login02.expanse.sdsc.edu) MACHINE_ID=expanse ;; ### expanse2
  
esac

case $(echo ${PW_CSP:-nono}) in

  aws) MACHINE_ID=aws ;; ### parallelworks aws
  google)  MACHINE_ID=gcp ;; ### parallelworks gcp
  azure)  MACHINE_ID=azure ;; ### parallelworks azure
  
esac
[[ ${MACHINE_ID} =~ "aws" || ${MACHINE_ID} =~ "gcp" || ${MACHINE_ID} =~ "azure" ]] && MACHINE_ID=noaacloud

# Overwrite auto-detect with RT_MACHINE if set
MACHINE_ID=${RT_MACHINE:-${MACHINE_ID}}
