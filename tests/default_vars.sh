
###############################################################################
#
# Export variables to the default values
#  - first common variables, then model specific ones
#  - different machines, different defaults:
#
###############################################################################

if [ $MACHINE_ID = wcoss ]; then

  TASKS_dflt=144 ; TPN_dflt=16 ; INPES_dflt=3 ; JNPES_dflt=8

elif [ $MACHINE_ID = wcoss_cray ]; then

  TASKS_dflt=144 ; TPN_dflt=24 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=72  ; TPN_thrd=12 ; INPES_thrd=3 ; JNPES_thrd=4

elif [ $MACHINE_ID = theia ]; then

  TASKS_dflt=144 ; TPN_dflt=24 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=72  ; TPN_thrd=12 ; INPES_thrd=3 ; JNPES_thrd=4


fi

export_fv3 ()
{
export THRD=1
export WLCLK=15
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export TASKS=$TASKS_dflt
export TPN=$TPN_dflt
export WARM_START=.F.
export NGGPS_IC=.T.
export EXTERNAL_IC=.T.
export MAKE_NH=.T.
export MOUNTAIN=.F.
export NA_INIT=1
export DAYS=1
export FDIAG=0,1,2,3,4,5,6,7,8,9,10,11,12,18,24


export ENS_NUM=1
export SYEAR=2016
export SMONTH=10
export SDAY=03
export SHOUR=00
export FHMAX=`expr $DAYS \* 24`
export DT_ATMOS=225

}
