#!/bin/bash

export ACCNR=${ACCNR:-nems}

case $(hostname -f) in

  g10a1.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=1;; ### gyre 1
  g10a2.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=1;; ### gyre 2
  g14a1.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=1;; ### gyre 3
  g14a2.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=1;; ### gyre 4

  t10a1.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=1;; ### tide 1
  t10a2.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=1;; ### tide 2
  t14a1.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=1;; ### tide 3
  t14a2.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=1;; ### tide 4

  g20a1.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=2;; ### gyre phase2
  g20a2.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=2;; ### gyre phase2
  g20a3.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=2;; ### gyre phase2
  g21a1.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=2;; ### gyre phase2
  g21a2.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=2;; ### gyre phase2
  g21a3.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=2;; ### gyre phase2

  t20a1.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=2;; ### tide phase2
  t20a2.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=2;; ### tide phase2
  t20a3.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=2;; ### tide phase2
  t21a1.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=2;; ### tide phase2
  t21a2.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=2;; ### tide phase2
  t21a3.ncep.noaa.gov)     MACHINE_ID=wcoss ; export pex=2;; ### tide phase2

  llogin1)                 MACHINE_ID=wcoss_cray ;; ### luna
  llogin2)                 MACHINE_ID=wcoss_cray ;; ### luna
  llogin3)                 MACHINE_ID=wcoss_cray ;; ### luna

  slogin1)                 MACHINE_ID=wcoss_cray ;; ### surge
  slogin2)                 MACHINE_ID=wcoss_cray ;; ### surge
  slogin3)                 MACHINE_ID=wcoss_cray ;; ### surge

  gaea1.ncrc.gov)          MACHINE_ID=gaea ;; ### gaea1
  gaea2.ncrc.gov)          MACHINE_ID=gaea ;; ### gaea2
  gaea3.ncrc.gov)          MACHINE_ID=gaea ;; ### gaea3
  gaea4.ncrc.gov)          MACHINE_ID=gaea ;; ### gaea4
  gaea5.ncrc.gov)          MACHINE_ID=gaea ;; ### gaea5
  gaea6.ncrc.gov)          MACHINE_ID=gaea ;; ### gaea6
  gaea7.ncrc.gov)          MACHINE_ID=gaea ;; ### gaea7
  gaea8.ncrc.gov)          MACHINE_ID=gaea ;; ### gaea8
  gaea9.ncrc.gov)          MACHINE_ID=gaea ;; ### gaea9
  gaea10.ncrc.gov)         MACHINE_ID=gaea ;; ### gaea10

  tfe01)                   MACHINE_ID=theia ;; ### theia01
  tfe02)                   MACHINE_ID=theia ;; ### theia02
  tfe03)                   MACHINE_ID=theia ;; ### theia03
  tfe04)                   MACHINE_ID=theia ;; ### theia04
  tfe05)                   MACHINE_ID=theia ;; ### theia05
  tfe06)                   MACHINE_ID=theia ;; ### theia06
  tfe07)                   MACHINE_ID=theia ;; ### theia07
  tfe08)                   MACHINE_ID=theia ;; ### theia08
  tfe09)                   MACHINE_ID=theia ;; ### theia09
  tfe10)                   MACHINE_ID=theia ;; ### theia10

  yslogin1)                MACHINE_ID=yellowstone ;;
  yslogin2)                MACHINE_ID=yellowstone ;;
  yslogin3)                MACHINE_ID=yellowstone ;;
  yslogin4)                MACHINE_ID=yellowstone ;;
  yslogin5)                MACHINE_ID=yellowstone ;;
  yslogin6)                MACHINE_ID=yellowstone ;;
  yslogin7)                MACHINE_ID=yellowstone ;;
  yslogin8)                MACHINE_ID=yellowstone ;;
  yslogin9)                MACHINE_ID=yellowstone ;;
  yslogin10)               MACHINE_ID=yellowstone ;;

esac

echo "Machine: " $MACHINE_ID "    Account: " $ACCNR

# --- for Theia, find available account ID
#if [[ $1"" != "machineonly" ]] ; then
  if [ ${MACHINE_ID} = theia ]; then
    AP=account_params          # Account info
    if [ ${ACCNR:-null} = null ]; then

      ac=`$AP 2>&1 | grep '^\s*Allocation: [0-9]' | awk '$4>100{print $3}'| head -1`
      nr=`echo $ac|wc -w`

      if [ $nr -eq 1 ]; then
        ACCNR=$ac
        echo "Found a valid account: using $ac"
      else
        ac=`$AP 2>&1 | grep '^\s*Allocation: [0-9]' | awk '{print $3}'| head -1`
        nr=`echo $ac|wc -w`
        if [ $nr -eq 1 ]; then
          ACCNR=$ac
          echo "Could not an find account with positive balance: using $ac"
          echo "NOTE: Will run in windfall; longer wait times, be patient!"
        else
          echo "Check your account ID; No compute allocations found"
        fi
      fi
    else
      cphr=`$AP 2>&1 | grep '^\s*Allocation: [0-9]' | grep $ACCNR | awk '{print $4}'`
      nr=`echo $cphr|wc -w`
      if [ $nr -eq 0 ]; then
        echo 'Wrong account choice: ' $ACCNR
      else
        echo "Account: " $ACCNR", available: " $cphr " CPU hrs"
      fi
    fi
  fi
#fi
