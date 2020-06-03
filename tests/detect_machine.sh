#!/bin/bash

# Default account "nems"
export ACCNR=${ACCNR:-nems}

case $(hostname -f) in

  llogin1)                 MACHINE_ID=wcoss_cray ;; ### luna
  llogin2)                 MACHINE_ID=wcoss_cray ;; ### luna
  llogin3)                 MACHINE_ID=wcoss_cray ;; ### luna

  slogin1)                 MACHINE_ID=wcoss_cray ;; ### surge
  slogin2)                 MACHINE_ID=wcoss_cray ;; ### surge
  slogin3)                 MACHINE_ID=wcoss_cray ;; ### surge

  v71a1.ncep.noaa.gov)     MACHINE_ID=wcoss_dell_p3 ;; ### venus
  v71a2.ncep.noaa.gov)     MACHINE_ID=wcoss_dell_p3 ;; ### venus
  v71a3.ncep.noaa.gov)     MACHINE_ID=wcoss_dell_p3 ;; ### venus
  v72a1.ncep.noaa.gov)     MACHINE_ID=wcoss_dell_p3 ;; ### venus
  v72a2.ncep.noaa.gov)     MACHINE_ID=wcoss_dell_p3 ;; ### venus
  v72a3.ncep.noaa.gov)     MACHINE_ID=wcoss_dell_p3 ;; ### venus

  m71a1.ncep.noaa.gov)     MACHINE_ID=wcoss_dell_p3 ;; ### mars
  m71a2.ncep.noaa.gov)     MACHINE_ID=wcoss_dell_p3 ;; ### mars
  m71a3.ncep.noaa.gov)     MACHINE_ID=wcoss_dell_p3 ;; ### mars
  m72a1.ncep.noaa.gov)     MACHINE_ID=wcoss_dell_p3 ;; ### mars
  m72a2.ncep.noaa.gov)     MACHINE_ID=wcoss_dell_p3 ;; ### mars
  m72a3.ncep.noaa.gov)     MACHINE_ID=wcoss_dell_p3 ;; ### mars

  gaea9)                   MACHINE_ID=gaea ;; ### gaea9
  gaea10)                  MACHINE_ID=gaea ;; ### gaea10
  gaea11)                  MACHINE_ID=gaea ;; ### gaea11
  gaea12)                  MACHINE_ID=gaea ;; ### gaea12
  gaea13)                  MACHINE_ID=gaea ;; ### gaea13
  gaea14)                  MACHINE_ID=gaea ;; ### gaea14
  gaea15)                  MACHINE_ID=gaea ;; ### gaea15
  gaea16)                  MACHINE_ID=gaea ;; ### gaea16
  gaea9.ncrc.gov)          MACHINE_ID=gaea ;; ### gaea9
  gaea10.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea10
  gaea11.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea11
  gaea12.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea12
  gaea13.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea13
  gaea14.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea14
  gaea15.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea15
  gaea16.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea16

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

  cheyenne1.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne1
  cheyenne2.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne2
  cheyenne3.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne3
  cheyenne4.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne4
  cheyenne5.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne5
  cheyenne6.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne6
  cheyenne1.ib0.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne1
  cheyenne2.ib0.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne2
  cheyenne3.ib0.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne3
  cheyenne4.ib0.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne4
  cheyenne5.ib0.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne5
  cheyenne6.ib0.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne6

  login1.stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede1
  login2.stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede2
  login3.stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede3
  login4.stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede4
esac

# Overwrite auto-detect with NEMS_MACHINE if set
MACHINE_ID=${NEMS_MACHINE:-${MACHINE_ID}}

# Append compiler
if [ $MACHINE_ID = orion ] || [ $MACHINE_ID = hera ] || [ $MACHINE_ID = cheyenne ] || [ $MACHINE_ID = jet ] || [ $MACHINE_ID = gaea ] || [ $MACHINE_ID = stampede ] ; then
    MACHINE_ID=${MACHINE_ID}.${COMPILER}
fi

echo "Machine: " $MACHINE_ID "    Account: " $ACCNR
