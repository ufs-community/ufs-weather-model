
###############################################################################
#
# Export variables to the default values
#  - first common variables, then model specific ones
#  - different machines, different defaults:
#
###############################################################################

if [[ $MACHINE_ID = wcoss_cray ]]; then

  TASKS_dflt=150 ; TPN_dflt=24 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=78  ; TPN_thrd=12 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_c384=336 ; TPN_c384=20 ; INPES_c384=8 ; JNPES_c384=6
  TASKS_stretch=48 ; TPN_stretch=24 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=24 ; INPES_strnest=2 ; JNPES_strnest=4

elif [[ $MACHINE_ID = wcoss_dell_p3 ]]; then

  TASKS_dflt=150 ; TPN_dflt=28 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=78  ; TPN_thrd=14 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_c384=336 ; TPN_c384=20 ; INPES_c384=8 ; JNPES_c384=6
  TASKS_stretch=48 ; TPN_stretch=28 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=28 ; INPES_strnest=2 ; JNPES_strnest=4

  TASKS_cpl_atmw=196; TPN_cpl_atmw=28; INPES_cpl_atmw=3; JNPES_cpl_atmw=8
  THRD_cpl_atmw=1; WPG_cpl_atmw=6; APB_cpl_atmw="0 149"; WPB_cpl_atmw="150 195"

  TASKS_cpl_atmw_gdas=560; TPN_cpl_atmw_gdas=14; INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8
  THRD_cpl_atmw_gdas=2; WPG_cpl_atmw_gdas=24; APB_cpl_atmw_gdas="0 311"; WPB_cpl_atmw_gdas="312 559"

  TASKS_cpl_c96=192; TPN_cpl_c96=28; INPES_cpl_c96=3; JNPES_cpl_c96=8
  THRD_cpl_c96=1; WPG_cpl_c96=6;  MPB_cpl_c96="0 143"; APB_cpl_c96="0 149"
  OPB_cpl_c96="150 179"; IPB_cpl_c96="180 191"
  NPROC_ICE_cpl_c96=12

  TASKS_cpl_dflt=196; TPN_cpl_dflt=28; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 169"; IPB_cpl_dflt="170 177"; WPB_cpl_dflt="178 195"  
  NPROC_ICE_cpl_dflt=8

  TASKS_cpl_thrd=120; TPN_cpl_thrd=14; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 71"; APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 97"; IPB_cpl_thrd="98 107"; WPB_cpl_thrd="108 119"
  NPROC_ICE_cpl_thrd=10

  TASKS_cpl_dcmp=196; TPN_cpl_dcmp=28; INPES_cpl_dcmp=4; JNPES_cpl_dcmp=6
  THRD_cpl_dcmp=1; WPG_cpl_dcmp=6;  MPB_cpl_dcmp="0 143"; APB_cpl_dcmp="0 149"
  OPB_cpl_dcmp="150 169"; IPB_cpl_dcmp="170 177"; WPB_cpl_dcmp="178 195"
  NPROC_ICE_cpl_dcmp=8

  TASKS_cpl_mpi=280; TPN_cpl_mpi=28; INPES_cpl_mpi=4; JNPES_cpl_mpi=8
  THRD_cpl_mpi=1; WPG_cpl_mpi=6;  MPB_cpl_mpi="0 191"; APB_cpl_mpi="0 197"
  OPB_cpl_mpi="198 231"; IPB_cpl_mpi="232 251"; WPB_cpl_mpi="252 279"
  NPROC_ICE_cpl_mpi=20

  TASKS_cpl_c384=480; TPN_cpl_c384=28; INPES_cpl_c384=6; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=24; MPB_cpl_c384="0 287"; APB_cpl_c384="0 311"
  OPB_cpl_c384="312 431"; IPB_cpl_c384="432 479"
  NPROC_ICE_cpl_c384=48

  TASKS_cpl_bmrk=560; TPN_cpl_bmrk=28; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"; WPB_cpl_bmrk="480 559"
  NPROC_ICE_cpl_bmrk=48

  TASKS_cpl_bmrk_mpi=600; TPN_cpl_bmrk_mpi=28; INPES_cpl_bmrk_mpi=6; JNPES_cpl_bmrk_mpi=8
  THRD_cpl_bmrk_mpi=1; WPG_cpl_bmrk_mpi=24; MPB_cpl_bmrk_mpi="0 287"; APB_cpl_bmrk_mpi="0 311"
  OPB_cpl_bmrk_mpi="312 431"; IPB_cpl_bmrk_mpi="432 479"; WPB_cpl_bmrk_mpi="480 599"
  NPROC_ICE_cpl_bmrk_mpi=48

  TASKS_cpl_c192=288; TPN_cpl_c192=28; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"
  NPROC_ICE_cpl_c192=24

  TASKS_cdeps_100=40; TPN_cdeps_100=28
  MPB_cdeps_100="0 11"; APB_cdeps_100="0 11"
  OPB_cdeps_100="12 27"; IPB_cdeps_100="28 39"

  TASKS_cdeps_025=208; TPN_cdeps_025=28
  MPB_cdeps_025="0 39"; APB_cdeps_025="0 39"
  OPB_cdeps_025="40 159"; IPB_cdeps_025="160 207"

elif [[ $MACHINE_ID = wcoss2 ]]; then

  TASKS_dflt=150 ; TPN_dflt=128 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=78  ; TPN_thrd=64 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_c384=336 ; TPN_c384=128 ; INPES_c384=8 ; JNPES_c384=6
  TASKS_stretch=48 ; TPN_stretch=128 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=128 ; INPES_strnest=2 ; JNPES_strnest=4

  TASKS_cpl_atmw=200; TPN_cpl_atmw=128; INPES_cpl_atmw=3; JNPES_cpl_atmw=8
  THRD_cpl_atmw=1; WPG_cpl_atmw=6; APB_cpl_atmw="0 149"; WPB_cpl_atmw="150 199"

  TASKS_cpl_atmw_gdas=560; TPN_cpl_atmw_gdas=32; INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8
  THRD_cpl_atmw_gdas=4; WPG_cpl_atmw_gdas=24; APB_cpl_atmw_gdas="0 311"; WPB_cpl_atmw_gdas="312 559"

  TASKS_cpl_c96=192; TPN_cpl_c96=128; INPES_cpl_c96=3; JNPES_cpl_c96=8
  THRD_cpl_c96=1; WPG_cpl_c96=6;  MPB_cpl_c96="0 143"; APB_cpl_c96="0 149"
  OPB_cpl_c96="150 179"; IPB_cpl_c96="180 191"
  NPROC_ICE_cpl_c96=12

  TASKS_cpl_dflt=200; TPN_cpl_dflt=128; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 169"; IPB_cpl_dflt="170 179"; WPB_cpl_dflt="180 199"
  NPROC_ICE_cpl_dflt=10

  TASKS_cpl_thrd=120; TPN_cpl_thrd=64; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 71"; APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 97"; IPB_cpl_thrd="98 107"; WPB_cpl_thrd="108 119"
  NPROC_ICE_cpl_thrd=10

  TASKS_cpl_dcmp=200; TPN_cpl_dcmp=128; INPES_cpl_dcmp=4; JNPES_cpl_dcmp=6
  THRD_cpl_dcmp=1; WPG_cpl_dcmp=6;  MPB_cpl_dcmp="0 143"; APB_cpl_dcmp="0 149"
  OPB_cpl_dcmp="150 169"; IPB_cpl_dcmp="170 179"; WPB_cpl_dcmp="180 199"
  NPROC_ICE_cpl_dcmp=10

  TASKS_cpl_mpi=280; TPN_cpl_mpi=128; INPES_cpl_mpi=4; JNPES_cpl_mpi=8
  THRD_cpl_mpi=1; WPG_cpl_mpi=6;  MPB_cpl_mpi="0 191"; APB_cpl_mpi="0 197"
  OPB_cpl_mpi="198 231"; IPB_cpl_mpi="232 251"; WPB_cpl_mpi="252 279"
  NPROC_ICE_cpl_mpi=20

  TASKS_cpl_c384=480; TPN_cpl_c384=128; INPES_cpl_c384=6; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=24; MPB_cpl_c384="0 287"; APB_cpl_c384="0 311"
  OPB_cpl_c384="312 431"; IPB_cpl_c384="432 479"
  NPROC_ICE_cpl_c384=48

  TASKS_cpl_bmrk=560; TPN_cpl_bmrk=128; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"; WPB_cpl_bmrk="480 559"
  NPROC_ICE_cpl_bmrk=48

  TASKS_cpl_bmrk_mpi=600; TPN_cpl_bmrk_mpi=40; INPES_cpl_bmrk_mpi=6; JNPES_cpl_bmrk_mpi=8
  THRD_cpl_bmrk_mpi=1; WPG_cpl_bmrk_mpi=24; MPB_cpl_bmrk_mpi="0 287"; APB_cpl_bmrk_mpi="0 311"
  OPB_cpl_bmrk_mpi="312 431"; IPB_cpl_bmrk_mpi="432 479"; WPB_cpl_bmrk_mpi="480 599"
  NPROC_ICE_cpl_bmrk_mpi=48

  TASKS_cpl_c192=288; TPN_cpl_c192=128; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"
  NPROC_ICE_cpl_c192=24

  TASKS_cdeps_100=40; TPN_cdeps_100=128
  MPB_cdeps_100="0 11"; APB_cdeps_100="0 11"
  OPB_cdeps_100="12 27"; IPB_cdeps_100="28 39"

  TASKS_cdeps_025=208; TPN_cdeps_025=128
  MPB_cdeps_025="0 39"; APB_cdeps_025="0 39"
  OPB_cdeps_025="40 159"; IPB_cdeps_025="160 207"

elif [[ $MACHINE_ID = orion.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=40 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=78  ; TPN_thrd=20 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_c384=336 ; TPN_c384=20 ; INPES_c384=8 ; JNPES_c384=6
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=12 ; INPES_strnest=2 ; JNPES_strnest=4

  TASKS_cpl_atmw=200; TPN_cpl_atmw=40; INPES_cpl_atmw=3; JNPES_cpl_atmw=8
  THRD_cpl_atmw=1; WPG_cpl_atmw=6; APB_cpl_atmw="0 149"; WPB_cpl_atmw="150 199"

  TASKS_cpl_atmw_gdas=560; TPN_cpl_atmw_gdas=20; INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8
  THRD_cpl_atmw_gdas=2; WPG_cpl_atmw_gdas=24; APB_cpl_atmw_gdas="0 311"; WPB_cpl_atmw_gdas="312 559"

  TASKS_cpl_c96=192; TPN_cpl_c96=40; INPES_cpl_c96=3; JNPES_cpl_c96=8
  THRD_cpl_c96=1; WPG_cpl_c96=6;  MPB_cpl_c96="0 143"; APB_cpl_c96="0 149"
  OPB_cpl_c96="150 179"; IPB_cpl_c96="180 191"
  NPROC_ICE_cpl_c96=12

  TASKS_cpl_dflt=200; TPN_cpl_dflt=40; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 169"; IPB_cpl_dflt="170 179"; WPB_cpl_dflt="180 199"  
  NPROC_ICE_cpl_dflt=10

  TASKS_cpl_thrd=120; TPN_cpl_thrd=20; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 71"; APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 97"; IPB_cpl_thrd="98 107"; WPB_cpl_thrd="108 119"
  NPROC_ICE_cpl_thrd=10

  TASKS_cpl_dcmp=200; TPN_cpl_dcmp=40; INPES_cpl_dcmp=4; JNPES_cpl_dcmp=6
  THRD_cpl_dcmp=1; WPG_cpl_dcmp=6;  MPB_cpl_dcmp="0 143"; APB_cpl_dcmp="0 149"
  OPB_cpl_dcmp="150 169"; IPB_cpl_dcmp="170 179"; WPB_cpl_dcmp="180 199"
  NPROC_ICE_cpl_dcmp=10

  TASKS_cpl_mpi=280; TPN_cpl_mpi=40; INPES_cpl_mpi=4; JNPES_cpl_mpi=8
  THRD_cpl_mpi=1; WPG_cpl_mpi=6;  MPB_cpl_mpi="0 191"; APB_cpl_mpi="0 197"
  OPB_cpl_mpi="198 231"; IPB_cpl_mpi="232 251"; WPB_cpl_mpi="252 279"
  NPROC_ICE_cpl_mpi=20

  TASKS_cpl_c384=480; TPN_cpl_c384=40; INPES_cpl_c384=6; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=24; MPB_cpl_c384="0 287"; APB_cpl_c384="0 311"
  OPB_cpl_c384="312 431"; IPB_cpl_c384="432 479"
  NPROC_ICE_cpl_c384=48

  TASKS_cpl_bmrk=560; TPN_cpl_bmrk=40; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"; WPB_cpl_bmrk="480 559"
  NPROC_ICE_cpl_bmrk=48

  TASKS_cpl_bmrk_mpi=600; TPN_cpl_bmrk_mpi=40; INPES_cpl_bmrk_mpi=6; JNPES_cpl_bmrk_mpi=8
  THRD_cpl_bmrk_mpi=1; WPG_cpl_bmrk_mpi=24; MPB_cpl_bmrk_mpi="0 287"; APB_cpl_bmrk_mpi="0 311"
  OPB_cpl_bmrk_mpi="312 431"; IPB_cpl_bmrk_mpi="432 479"; WPB_cpl_bmrk_mpi="480 599"
  NPROC_ICE_cpl_bmrk_mpi=48

  TASKS_cpl_c192=288; TPN_cpl_c192=40; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"
  NPROC_ICE_cpl_c192=24

  TASKS_cdeps_100=40; TPN_cdeps_100=40
  MPB_cdeps_100="0 11"; APB_cdeps_100="0 11"
  OPB_cdeps_100="12 27"; IPB_cdeps_100="28 39"

  TASKS_cdeps_025=208; TPN_cdeps_025=40
  MPB_cdeps_025="0 39"; APB_cdeps_025="0 39"
  OPB_cdeps_025="40 159"; IPB_cdeps_025="160 207"

elif [[ $MACHINE_ID = hera.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=40 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=78  ; TPN_thrd=20 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_c384=336 ; TPN_c384=20 ; INPES_c384=6 ; JNPES_c384=8
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=12 ; INPES_strnest=2 ; JNPES_strnest=4

  TASKS_cpl_atmw=200; TPN_cpl_atmw=40; INPES_cpl_atmw=3; JNPES_cpl_atmw=8
  THRD_cpl_atmw=1; WPG_cpl_atmw=6; APB_cpl_atmw="0 149"; WPB_cpl_atmw="150 199"

  TASKS_cpl_atmw_gdas=560; TPN_cpl_atmw_gdas=20; INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8
  THRD_cpl_atmw_gdas=2; WPG_cpl_atmw_gdas=24; APB_cpl_atmw_gdas="0 311"; WPB_cpl_atmw_gdas="312 559"

  TASKS_cpl_c96=192; TPN_cpl_c96=40; INPES_cpl_c96=3; JNPES_cpl_c96=8
  THRD_cpl_c96=1; WPG_cpl_c96=6;  MPB_cpl_c96="0 143"; APB_cpl_c96="0 149"
  OPB_cpl_c96="150 179"; IPB_cpl_c96="180 191"
  NPROC_ICE_cpl_c96=12

  TASKS_cpl_dflt=200; TPN_cpl_dflt=40; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 169"; IPB_cpl_dflt="170 179"; WPB_cpl_dflt="180 199"
  NPROC_ICE_cpl_dflt=10

  TASKS_cpl_thrd=120; TPN_cpl_thrd=20; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 71"; APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 97"; IPB_cpl_thrd="98 107"; WPB_cpl_thrd="108 119"
  NPROC_ICE_cpl_thrd=10

  TASKS_cpl_dcmp=200; TPN_cpl_dcmp=40; INPES_cpl_dcmp=4; JNPES_cpl_dcmp=6
  THRD_cpl_dcmp=1; WPG_cpl_dcmp=6;  MPB_cpl_dcmp="0 143"; APB_cpl_dcmp="0 149"
  OPB_cpl_dcmp="150 169"; IPB_cpl_dcmp="170 179"; WPB_cpl_dcmp="180 199"
  NPROC_ICE_cpl_dcmp=10

  TASKS_cpl_mpi=280; TPN_cpl_mpi=40; INPES_cpl_mpi=4; JNPES_cpl_mpi=8
  THRD_cpl_mpi=1; WPG_cpl_mpi=6;  MPB_cpl_mpi="0 191"; APB_cpl_mpi="0 197"
  OPB_cpl_mpi="198 231"; IPB_cpl_mpi="232 251"; WPB_cpl_mpi="252 279"
  NPROC_ICE_cpl_mpi=20

  TASKS_cpl_c384=480; TPN_cpl_c384=40; INPES_cpl_c384=6; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=24; MPB_cpl_c384="0 287"; APB_cpl_c384="0 311"
  OPB_cpl_c384="312 431"; IPB_cpl_c384="432 479"
  NPROC_ICE_cpl_c384=48

  TASKS_cpl_bmrk=560; TPN_cpl_bmrk=40; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"; WPB_cpl_bmrk="480 559"
  NPROC_ICE_cpl_bmrk=48

  TASKS_cpl_bmrk_mpi=600; TPN_cpl_bmrk_mpi=40; INPES_cpl_bmrk_mpi=6; JNPES_cpl_bmrk_mpi=8
  THRD_cpl_bmrk_mpi=1; WPG_cpl_bmrk_mpi=24; MPB_cpl_bmrk_mpi="0 287"; APB_cpl_bmrk_mpi="0 311"
  OPB_cpl_bmrk_mpi="312 431"; IPB_cpl_bmrk_mpi="432 479"; WPB_cpl_bmrk_mpi="480 599"
  NPROC_ICE_cpl_bmrk_mpi=48

  TASKS_cpl_c192=288; TPN_cpl_c192=40; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"
  NPROC_ICE_cpl_c192=24

  TASKS_cdeps_100=40; TPN_cdeps_100=40
  MPB_cdeps_100="0 11"; APB_cdeps_100="0 11"
  OPB_cdeps_100="12 27"; IPB_cdeps_100="28 39"

  TASKS_cdeps_025=208; TPN_cdeps_025=40
  MPB_cdeps_025="0 39"; APB_cdeps_025="0 39"
  OPB_cdeps_025="40 159"; IPB_cdeps_025="160 207"

elif [[ $MACHINE_ID = linux.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=40 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=78  ; TPN_thrd=20 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=12 ; INPES_strnest=2 ; JNPES_strnest=4

  TASKS_cpl_c96=192; TPN_cpl_c96=40; INPES_cpl_c96=3; JNPES_cpl_c96=8
  THRD_cpl_c96=1; WPG_cpl_c96=6;  MPB_cpl_c96="0 143"; APB_cpl_c96="0 149"
  OPB_cpl_c96="150 179"; IPB_cpl_c96="180 191"
  NPROC_ICE_cpl_c96=12

  TASKS_cpl_dflt=200; TPN_cpl_dflt=40; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 169"; IPB_cpl_dflt="170 179"; WPB_cpl_dflt="180 199"
  NPROC_ICE_cpl_dflt=10

  TASKS_cpl_thrd=120; TPN_cpl_thrd=40; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 71";  APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 97";  IPB_cpl_thrd="98 107"; WPB_cpl_thrd="108 119"
  NPROC_ICE_cpl_thrd=10

elif [[ $MACHINE_ID = jet.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=24 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=78  ; TPN_thrd=12 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_c384=480 ; TPN_c384=12 ; INPES_c384=12 ; JNPES_c384=6
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=12 ; INPES_strnest=2 ; JNPES_strnest=4

  TASKS_cpl_atmw=192; TPN_cpl_atmw=24; INPES_cpl_atmw=3; JNPES_cpl_atmw=8
  THRD_cpl_atmw=1; WPG_cpl_atmw=6; APB_cpl_atmw="0 149"; WPB_cpl_atmw="150 191"

  TASKS_cpl_atmw_gdas=552; TPN_cpl_atmw_gdas=12; INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8
  THRD_cpl_atmw_gdas=2; WPG_cpl_atmw_gdas=24; APB_cpl_atmw_gdas="0 311"; WPB_cpl_atmw_gdas="312 551"

  TASKS_cpl_c96=192; TPN_cpl_c96=24; INPES_cpl_c96=3; JNPES_cpl_c96=8
  THRD_cpl_c96=1; WPG_cpl_c96=6;  MPB_cpl_c96="0 143"; APB_cpl_c96="0 149"
  OPB_cpl_c96="150 179"; IPB_cpl_c96="180 191"
  NPROC_ICE_cpl_c96=12

  TASKS_cpl_dflt=200; TPN_cpl_dflt=24; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 169"; IPB_cpl_dflt="170 179"; WPB_cpl_dflt="180 199"
  NPROC_ICE_cpl_dflt=10

  TASKS_cpl_thrd=120; TPN_cpl_thrd=12; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 71"; APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 97"; IPB_cpl_thrd="98 107"; WPB_cpl_thrd="108 119"
  NPROC_ICE_cpl_thrd=10

  TASKS_cpl_dcmp=200; TPN_cpl_dcmp=24; INPES_cpl_dcmp=4; JNPES_cpl_dcmp=6
  THRD_cpl_dcmp=1; WPG_cpl_dcmp=6;  MPB_cpl_dcmp="0 143"; APB_cpl_dcmp="0 149"
  OPB_cpl_dcmp="150 169"; IPB_cpl_dcmp="170 179"; WPB_cpl_dcmp="180 199"
  NPROC_ICE_cpl_dcmp=10

  TASKS_cpl_mpi=280; TPN_cpl_mpi=24; INPES_cpl_mpi=4; JNPES_cpl_mpi=8
  THRD_cpl_mpi=1; WPG_cpl_mpi=6;  MPB_cpl_mpi="0 191"; APB_cpl_mpi="0 197"
  OPB_cpl_mpi="198 231"; IPB_cpl_mpi="232 251"; WPB_cpl_mpi="252 279"
  NPROC_ICE_cpl_mpi=20

  TASKS_cpl_c384=480; TPN_cpl_c384=24; INPES_cpl_c384=6; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=24; MPB_cpl_c384="0 287"; APB_cpl_c384="0 311"
  OPB_cpl_c384="312 431"; IPB_cpl_c384="432 479"
  NPROC_ICE_cpl_c384=48

  TASKS_cpl_bmrk=560; TPN_cpl_bmrk=24; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"; WPB_cpl_bmrk="480 559"
  NPROC_ICE_cpl_bmrk=48

  TASKS_cpl_bmrk_mpi=600; TPN_cpl_bmrk_mpi=24; INPES_cpl_bmrk_mpi=6; JNPES_cpl_bmrk_mpi=8
  THRD_cpl_bmrk_mpi=1; WPG_cpl_bmrk_mpi=24; MPB_cpl_bmrk_mpi="0 287"; APB_cpl_bmrk_mpi="0 311"
  OPB_cpl_bmrk_mpi="312 431"; IPB_cpl_bmrk_mpi="432 479"; WPB_cpl_bmrk_mpi="480 599"
  NPROC_ICE_cpl_bmrk_mpi=48

  TASKS_cpl_c192=288; TPN_cpl_c192=24; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"
  NPROC_ICE_cpl_c192=24

  TASKS_cdeps_100=40; TPN_cdeps_100=24
  MPB_cdeps_100="0 11"; APB_cdeps_100="0 11"
  OPB_cdeps_100="12 27"; IPB_cdeps_100="28 39"

  TASKS_cdeps_025=208; TPN_cdeps_025=24
  MPB_cdeps_025="0 39"; APB_cdeps_025="0 39"
  OPB_cdeps_025="40 159"; IPB_cdeps_025="160 207"

elif [[ $MACHINE_ID = s4.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=32 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=78  ; TPN_thrd=16 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_c384=336 ; TPN_c384=16 ; INPES_c384=6 ; JNPES_c384=8
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=12 ; INPES_strnest=2 ; JNPES_strnest=4

  TASKS_cpl_atmw=200; TPN_cpl_atmw=16; INPES_cpl_atmw=3; JNPES_cpl_atmw=8
  THRD_cpl_atmw=2; WPG_cpl_atmw=6; APB_cpl_atmw="0 149"; WPB_cpl_atmw="150 199"

  TASKS_cpl_atmw_gdas=560; TPN_cpl_atmw_gdas=16; INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8
  THRD_cpl_atmw_gdas=2; WPG_cpl_atmw_gdas=24; APB_cpl_atmw_gdas="0 311"; WPB_cpl_atmw_gdas="312 559"

  TASKS_cpl_c96=192; TPN_cpl_c96=32; INPES_cpl_c96=3; JNPES_cpl_c96=8
  THRD_cpl_c96=1; WPG_cpl_c96=6;  MPB_cpl_c96="0 143"; APB_cpl_c96="0 149"
  OPB_cpl_c96="150 179"; IPB_cpl_c96="180 191"
  NPROC_ICE_cpl_c96=12

  TASKS_cpl_dflt=200; TPN_cpl_dflt=32; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 169"; IPB_cpl_dflt="170 179"; WPB_cpl_dflt="180 199"
  NPROC_ICE_cpl_dflt=10

  TASKS_cpl_thrd=120; TPN_cpl_thrd=16; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 71"; APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 97"; IPB_cpl_thrd="98 107"; WPB_cpl_thrd="108 119"
  NPROC_ICE_cpl_thrd=10

  TASKS_cpl_dcmp=200; TPN_cpl_dcmp=32; INPES_cpl_dcmp=4; JNPES_cpl_dcmp=6
  THRD_cpl_dcmp=1; WPG_cpl_dcmp=6;  MPB_cpl_dcmp="0 143"; APB_cpl_dcmp="0 149"
  OPB_cpl_dcmp="150 169"; IPB_cpl_dcmp="170 179"; WPB_cpl_dcmp="180 199"
  NPROC_ICE_cpl_dcmp=10

  TASKS_cpl_mpi=280; TPN_cpl_mpi=32; INPES_cpl_mpi=4; JNPES_cpl_mpi=8
  THRD_cpl_mpi=1; WPG_cpl_mpi=6;  MPB_cpl_mpi="0 191"; APB_cpl_mpi="0 197"
  OPB_cpl_mpi="198 231"; IPB_cpl_mpi="232 251"; WPB_cpl_mpi="252 279"
  NPROC_ICE_cpl_mpi=20

  TASKS_cpl_c384=480; TPN_cpl_c384=32; INPES_cpl_c384=6; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=24; MPB_cpl_c384="0 287"; APB_cpl_c384="0 311"
  OPB_cpl_c384="312 431"; IPB_cpl_c384="432 479"
  NPROC_ICE_cpl_c384=48

  TASKS_cpl_bmrk=560; TPN_cpl_bmrk=32; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"; WPB_cpl_bmrk="480 559"
  NPROC_ICE_cpl_bmrk=48

  TASKS_cpl_bmrk_mpi=600; TPN_cpl_bmrk_mpi=32; INPES_cpl_bmrk_mpi=6; JNPES_cpl_bmrk_mpi=8
  THRD_cpl_bmrk_mpi=1; WPG_cpl_bmrk_mpi=24; MPB_cpl_bmrk_mpi="0 287"; APB_cpl_bmrk_mpi="0 311"
  OPB_cpl_bmrk_mpi="312 431"; IPB_cpl_bmrk_mpi="432 479"; WPB_cpl_bmrk_mpi="480 599"
  NPROC_ICE_cpl_bmrk_mpi=48

  TASKS_cpl_c192=288; TPN_cpl_c192=32; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"
  NPROC_ICE_cpl_c192=24

  TASKS_cdeps_100=40; TPN_cdeps_100=32
  MPB_cdeps_100="0 11"; APB_cdeps_100="0 11"
  OPB_cdeps_100="12 27"; IPB_cdeps_100="28 39"

  TASKS_cdeps_025=208; TPN_cdeps_025=32
  MPB_cdeps_025="0 39"; APB_cdeps_025="0 39"
  OPB_cdeps_025="40 159"; IPB_cdeps_025="160 207"

elif [[ $MACHINE_ID = gaea.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=36 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=78  ; TPN_thrd=18 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_c384=336 ; TPN_c384=18 ; INPES_c384=6 ; JNPES_c384=8
  TASKS_stretch=48 ; TPN_stretch=18 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=18 ; INPES_strnest=2 ; JNPES_strnest=4

  TASKS_cpl_atmw=180; TPN_cpl_atmw=36; INPES_cpl_atmw=3; JNPES_cpl_atmw=8
  THRD_cpl_atmw=1; WPG_cpl_atmw=6; APB_cpl_atmw="0 149"; WPB_cpl_atmw="150 179"

  TASKS_cpl_atmw_gdas=576; TPN_cpl_atmw_gdas=12; INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8
  THRD_cpl_atmw_gdas=3; WPG_cpl_atmw_gdas=24; APB_cpl_atmw_gdas="0 311"; WPB_cpl_atmw_gdas="312 575"

  TASKS_cpl_c96=192; TPN_cpl_c96=36; INPES_cpl_c96=3; JNPES_cpl_c96=8
  THRD_cpl_c96=1; WPG_cpl_c96=6;  MPB_cpl_c96="0 143"; APB_cpl_c96="0 149"
  OPB_cpl_c96="150 179"; IPB_cpl_c96="180 191"
  NPROC_ICE_cpl_c96=12

  TASKS_cpl_dflt=216; TPN_cpl_dflt=36; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 169"; IPB_cpl_dflt="170 181"; WPB_cpl_dflt="182 215"
  NPROC_ICE_cpl_dflt=12

  TASKS_cpl_thrd=120; TPN_cpl_thrd=18; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 71"; APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 97"; IPB_cpl_thrd="98 107"; WPB_cpl_thrd="108 119"
  NPROC_ICE_cpl_thrd=10

  TASKS_cpl_dcmp=200; TPN_cpl_dcmp=36; INPES_cpl_dcmp=4; JNPES_cpl_dcmp=6
  THRD_cpl_dcmp=1; WPG_cpl_dcmp=6;  MPB_cpl_dcmp="0 143"; APB_cpl_dcmp="0 149"
  OPB_cpl_dcmp="150 169"; IPB_cpl_dcmp="170 179"; WPB_cpl_dcmp="180 199"
  NPROC_ICE_cpl_dcmp=10

  TASKS_cpl_mpi=280; TPN_cpl_mpi=36; INPES_cpl_mpi=4; JNPES_cpl_mpi=8
  THRD_cpl_mpi=1; WPG_cpl_mpi=6;  MPB_cpl_mpi="0 191"; APB_cpl_mpi="0 197"
  OPB_cpl_mpi="198 231"; IPB_cpl_mpi="232 251"; WPB_cpl_mpi="252 279"
  NPROC_ICE_cpl_mpi=20

  TASKS_cpl_c384=480; TPN_cpl_c384=36; INPES_cpl_c384=6; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=24; MPB_cpl_c384="0 287"; APB_cpl_c384="0 311"
  OPB_cpl_c384="312 431"; IPB_cpl_c384="432 479"
  NPROC_ICE_cpl_c384=48

  TASKS_cpl_bmrk=560; TPN_cpl_bmrk=36; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"; WPB_cpl_bmrk="480 559"
  NPROC_ICE_cpl_bmrk=48

  TASKS_cpl_bmrk_mpi=600; TPN_cpl_bmrk_mpi=36; INPES_cpl_bmrk_mpi=6; JNPES_cpl_bmrk_mpi=8
  THRD_cpl_bmrk_mpi=1; WPG_cpl_bmrk_mpi=24; MPB_cpl_bmrk_mpi="0 287"; APB_cpl_bmrk_mpi="0 311"
  OPB_cpl_bmrk_mpi="312 431"; IPB_cpl_bmrk_mpi="432 479"; WPB_cpl_bmrk_mpi="480 599"
  NPROC_ICE_cpl_bmrk_mpi=48

  TASKS_cpl_c192=288; TPN_cpl_c192=36; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"
  NPROC_ICE_cpl_c192=24

  TASKS_cdeps_100=40; TPN_cdeps_100=36
  MPB_cdeps_100="0 11"; APB_cdeps_100="0 11"
  OPB_cdeps_100="12 27"; IPB_cdeps_100="28 39"

  TASKS_cdeps_025=208; TPN_cdeps_025=36
  MPB_cdeps_025="0 39"; APB_cdeps_025="0 39"
  OPB_cdeps_025="40 159"; IPB_cdeps_025="160 207"

elif [[ $MACHINE_ID = cheyenne.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=36 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=78  ; TPN_thrd=18 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_c384=336 ; TPN_c384=18 ; INPES_c384=8 ; JNPES_c384=6
  TASKS_stretch=48 ; TPN_stretch=18 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=18 ; INPES_strnest=2 ; JNPES_strnest=4

  TASKS_cpl_atmw=180; TPN_cpl_atmw=40; INPES_cpl_atmw=3; JNPES_cpl_atmw=8
  THRD_cpl_atmw=1; WPG_cpl_atmw=6; APB_cpl_atmw="0 149"; WPB_cpl_atmw="150 179"

  TASKS_cpl_atmw_gdas=576; TPN_cpl_atmw_gdas=12; INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8
  THRD_cpl_atmw_gdas=3; WPG_cpl_atmw_gdas=24; APB_cpl_atmw_gdas="0 311"; WPB_cpl_atmw_gdas="312 575"

  TASKS_cpl_c96=192; TPN_cpl_c96=36; INPES_cpl_c96=3; JNPES_cpl_c96=8
  THRD_cpl_c96=1; WPG_cpl_c96=6;  MPB_cpl_c96="0 143"; APB_cpl_c96="0 149"
  OPB_cpl_c96="150 179"; IPB_cpl_c96="180 191"
  NPROC_ICE_cpl_c96=12

  TASKS_cpl_dflt=200; TPN_cpl_dflt=36; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 169"; IPB_cpl_dflt="170 179"; WPB_cpl_dflt="180 199"
  NPROC_ICE_cpl_dflt=10

  TASKS_cpl_thrd=120; TPN_cpl_thrd=18; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 71"; APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 97"; IPB_cpl_thrd="98 107"; WPB_cpl_thrd="108 119"
  NPROC_ICE_cpl_thrd=10

  TASKS_cpl_dcmp=200; TPN_cpl_dcmp=36; INPES_cpl_dcmp=4; JNPES_cpl_dcmp=6
  THRD_cpl_dcmp=1; WPG_cpl_dcmp=6;  MPB_cpl_dcmp="0 143"; APB_cpl_dcmp="0 149"
  OPB_cpl_dcmp="150 169"; IPB_cpl_dcmp="170 179"; WPB_cpl_dcmp="180 199"
  NPROC_ICE_cpl_dcmp=10

  TASKS_cpl_mpi=280; TPN_cpl_mpi=36; INPES_cpl_mpi=4; JNPES_cpl_mpi=8
  THRD_cpl_mpi=1; WPG_cpl_mpi=6;  MPB_cpl_mpi="0 191"; APB_cpl_mpi="0 197"
  OPB_cpl_mpi="198 231"; IPB_cpl_mpi="232 251"; WPB_cpl_mpi="252 279"
  NPROC_ICE_cpl_mpi=20

  TASKS_cpl_c384=480; TPN_cpl_c384=36; INPES_cpl_c384=6; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=24; MPB_cpl_c384="0 287"; APB_cpl_c384="0 311"
  OPB_cpl_c384="312 431"; IPB_cpl_c384="432 479"
  NPROC_ICE_cpl_c384=48

  TASKS_cpl_bmrk=560; TPN_cpl_bmrk=36; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"; WPB_cpl_bmrk="480 559"
  NPROC_ICE_cpl_bmrk=48

  TASKS_cpl_bmrk_mpi=600; TPN_cpl_bmrk_mpi=36; INPES_cpl_bmrk_mpi=6; JNPES_cpl_bmrk_mpi=8
  THRD_cpl_bmrk_mpi=1; WPG_cpl_bmrk_mpi=24; MPB_cpl_bmrk_mpi="0 287"; APB_cpl_bmrk_mpi="0 311"
  OPB_cpl_bmrk_mpi="312 431"; IPB_cpl_bmrk_mpi="432 479"; WPB_cpl_bmrk_mpi="480 599"
  NPROC_ICE_cpl_bmrk_mpi=48

  TASKS_cpl_c192=288; TPN_cpl_c192=36; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"
  NPROC_ICE_cpl_c192=24

  TASKS_cdeps_100=40; TPN_cdeps_100=36
  MPB_cdeps_100="0 11"; APB_cdeps_100="0 11"
  OPB_cdeps_100="12 27"; IPB_cdeps_100="28 39"

  TASKS_cdeps_025=208; TPN_cdeps_025=36
  MPB_cdeps_025="0 39"; APB_cdeps_025="0 39"
  OPB_cdeps_025="40 159"; IPB_cdeps_025="160 207"

elif [[ $MACHINE_ID = stampede.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=48 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=78  ; TPN_thrd=24 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_c384=336 ; TPN_c384=20 ; INPES_c384=8 ; JNPES_c384=6
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4

  TASKS_cpl_atmw=192; TPN_cpl_atmw=48; INPES_cpl_atmw=3; JNPES_cpl_atmw=8
  THRD_cpl_atmw=1; WPG_cpl_atmw=6; APB_cpl_atmw="0 149"; WPB_cpl_atmw="150 191"

  TASKS_cpl_atmw_gdas=560; TPN_cpl_atmw_gdas=12; INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8
  THRD_cpl_atmw_gdas=4; WPG_cpl_atmw_gdas=24; APB_cpl_atmw_gdas="0 311"; WPB_cpl_atmw_gdas="312 559"

  TASKS_cpl_c96=192; TPN_cpl_c96=48; INPES_cpl_c96=3; JNPES_cpl_c96=8
  THRD_cpl_c96=1; WPG_cpl_c96=6;  MPB_cpl_c96="0 143"; APB_cpl_c96="0 149"
  OPB_cpl_c96="150 179"; IPB_cpl_c96="180 191"
  NPROC_ICE_cpl_c96=12

  TASKS_cpl_dflt=200; TPN_cpl_dflt=48; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 169"; IPB_cpl_dflt="170 179"; WPB_cpl_dflt="180 199"
  NPROC_ICE_cpl_dflt=10

  TASKS_cpl_thrd=120; TPN_cpl_thrd=24; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 71"; APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 97"; IPB_cpl_thrd="98 107"; WPB_cpl_thrd="108 119"
  NPROC_ICE_cpl_thrd=10

  TASKS_cpl_dcmp=200; TPN_cpl_dcmp=48; INPES_cpl_dcmp=4; JNPES_cpl_dcmp=6
  THRD_cpl_dcmp=1; WPG_cpl_dcmp=6;  MPB_cpl_dcmp="0 143"; APB_cpl_dcmp="0 149"
  OPB_cpl_dcmp="150 169"; IPB_cpl_dcmp="170 179"; WPB_cpl_dcmp="180 199"
  NPROC_ICE_cpl_dcmp=10

  TASKS_cpl_mpi=280; TPN_cpl_mpi=48; INPES_cpl_mpi=4; JNPES_cpl_mpi=8
  THRD_cpl_mpi=1; WPG_cpl_mpi=6;  MPB_cpl_mpi="0 191"; APB_cpl_mpi="0 197"
  OPB_cpl_mpi="198 231"; IPB_cpl_mpi="232 251"; WPB_cpl_mpi="252 279"
  NPROC_ICE_cpl_mpi=20

  TASKS_cpl_c384=480; TPN_cpl_c384=48; INPES_cpl_c384=6; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=24; MPB_cpl_c384="0 287"; APB_cpl_c384="0 311"
  OPB_cpl_c384="312 431"; IPB_cpl_c384="432 479"
  NPROC_ICE_cpl_c384=48

  TASKS_cpl_bmrk=560; TPN_cpl_bmrk=48; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"; WPB_cpl_bmrk="480 559"
  NPROC_ICE_cpl_bmrk=48

  TASKS_cpl_bmrk_mpi=600; TPN_cpl_bmrk_mpi=40; INPES_cpl_bmrk_mpi=6; JNPES_cpl_bmrk_mpi=8
  THRD_cpl_bmrk_mpi=1; WPG_cpl_bmrk_mpi=24; MPB_cpl_bmrk_mpi="0 287"; APB_cpl_bmrk_mpi="0 311"
  OPB_cpl_bmrk_mpi="312 431"; IPB_cpl_bmrk_mpi="432 479"; WPB_cpl_bmrk_mpi="480 599"
  NPROC_ICE_cpl_bmrk_mpi=48

  TASKS_cpl_c192=288; TPN_cpl_c192=48; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"
  NPROC_ICE_cpl_c192=24

  TASKS_cdeps_100=40; TPN_cdeps_100=48
  MPB_cdeps_100="0 11"; APB_cdeps_100="0 11"
  OPB_cdeps_100="12 27"; IPB_cdeps_100="28 39"

  TASKS_cdeps_025=208; TPN_cdeps_025=48
  MPB_cdeps_025="0 39"; APB_cdeps_025="0 39"
  OPB_cdeps_025="40 159"; IPB_cdeps_025="160 207"

elif [[ $MACHINE_ID = expanse.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=64 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=78  ; TPN_thrd=64 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4

  TASKS_cpl_atmw=200; TPN_cpl_atmw=64; INPES_cpl_atmw=3; JNPES_cpl_atmw=8
  THRD_cpl_atmw=1; WPG_cpl_atmw=6; APB_cpl_atmw="0 149"; WPB_cpl_atmw="150 199"

  TASKS_cpl_atmw_gdas=560; TPN_cpl_atmw_gdas=12; INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8
  THRD_cpl_atmw_gdas=2; WPG_cpl_atmw_gdas=24; APB_cpl_atmw_gdas="0 311"; WPB_cpl_atmw_gdas="312 559"

  TASKS_cpl_c96=192; TPN_cpl_c96=64; INPES_cpl_c96=3; JNPES_cpl_c96=8
  THRD_cpl_c96=1; WPG_cpl_c96=6;  MPB_cpl_c96="0 143"; APB_cpl_c96="0 149"
  OPB_cpl_c96="150 179"; IPB_cpl_c96="180 191"
  NPROC_ICE_cpl_c96=12

  TASKS_cpl_dflt=384; TPN_cpl_dflt=64; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 179"; IPB_cpl_dflt="180 191"; WPB_cpl_dflt="192 383"
  NPROC_ICE_cpl_dflt=10

  TASKS_cpl_thrd=120; TPN_cpl_thrd=32; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 71"; APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 97"; IPB_cpl_thrd="98 107"; WPB_cpl_thrd="108 119"
  NPROC_ICE_cpl_thrd=10

  TASKS_cpl_dcmp=200; TPN_cpl_dcmp=64; INPES_cpl_dcmp=4; JNPES_cpl_dcmp=6
  THRD_cpl_dcmp=1; WPG_cpl_dcmp=6;  MPB_cpl_dcmp="0 143"; APB_cpl_dcmp="0 149"
  OPB_cpl_dcmp="150 169"; IPB_cpl_dcmp="170 179"; WPB_cpl_dcmp="180 199"
  NPROC_ICE_cpl_dcmp=10

  TASKS_cpl_mpi=280; TPN_cpl_mpi=64; INPES_cpl_mpi=4; JNPES_cpl_mpi=8
  THRD_cpl_mpi=1; WPG_cpl_mpi=6;  MPB_cpl_mpi="0 191"; APB_cpl_mpi="0 197"
  OPB_cpl_mpi="198 231"; IPB_cpl_mpi="232 251"; WPB_cpl_mpi="252 279"
  NPROC_ICE_cpl_mpi=20

  TASKS_cpl_bmrk=480; TPN_cpl_bmrk=64; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"
  NPROC_ICE_cpl_bmrk=48

  TASKS_cpl_wwav=640; TPN_cpl_wwav=64; INPES_cpl_wwav=6; JNPES_cpl_wwav=8
  THRD_cpl_wwav=2; WPG_cpl_wwav=24; MPB_cpl_wwav="0 287"; APB_cpl_wwav="0 311"
  OPB_cpl_wwav="312 431"; IPB_cpl_wwav="432 479"; WPB_cpl_wwav="480 639"

  TASKS_cpl_bmrk_mpi=600; TPN_cpl_bmrk_mpi=40; INPES_cpl_bmrk_mpi=6; JNPES_cpl_bmrk_mpi=8
  THRD_cpl_bmrk_mpi=1; WPG_cpl_bmrk_mpi=24; MPB_cpl_bmrk_mpi="0 287"; APB_cpl_bmrk_mpi="0 311"
  OPB_cpl_bmrk_mpi="312 431"; IPB_cpl_bmrk_mpi="432 479"; WPB_cpl_bmrk_mpi="480 599"
  NPROC_ICE_cpl_bmrk_mpi=48

  TASKS_cpl_c192=288; TPN_cpl_c192=64; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"
  NPROC_ICE_cpl_c192=24

  TASKS_cpl_c384=318; TPN_cpl_c384=64; INPES_cpl_c384=3; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=6;  MPB_cpl_c384="0 143"; APB_cpl_c384="0 149"
  OPB_cpl_c384="150 269"; IPB_cpl_c384="270 317"
  NPROC_ICE_cpl_c384=48

else

  echo "Unknown MACHINE_ID ${MACHINE_ID}"
  exit 1

fi

WLCLK_dflt=30
# Longer default walltime on Gaea
if [[ $MACHINE_ID = gaea.* ]]; then
  WLCLK_dflt=180
fi

export WLCLK=$WLCLK_dflt

export_fv3 ()
{
export FV3=true
export S2S=false
export HAFS=false
export DATM_CDEPS=false
export DOCN_CDEPS=false
export THRD=1
export POSTAPP='global'
export NEW_DIAGTABLE='none'
export NEW_FIELDTABLE='none'
export USE_MERRA2=.false.

export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export TASKS=$TASKS_dflt
export TPN=$TPN_dflt
export RESTART_INTERVAL=0
export QUILTING=.true.
export WRITE_GROUP=1
export WRTTASK_PER_GROUP=6
export OUTPUT_HISTORY=.true.
export WRITE_DOPOST=.false.
export NUM_FILES=2
export FILENAME_BASE="'atm' 'sfc'"
export OUTPUT_GRID="'cubed_sphere_grid'"
export OUTPUT_FILE="'netcdf'"
export IDEFLATE=0
export NBITS=0
export IMO=384
export JMO=190

#input file
export FIELD_TABLE=field_table_gfsv16

# Coldstart/warmstart
#rt script for ICs
export MODEL_INITIALIZATION=false
#namelist variable
export WARM_START=.false.
export READ_INCREMENT=.false.
export RES_LATLON_DYNAMICS="''"
export NGGPS_IC=.true.
export EXTERNAL_IC=.true.
export MAKE_NH=.true.
export MOUNTAIN=.false.
export NA_INIT=1

# Radiation
export DO_RRTMGP=.false.
export ICLOUD=0
export IAER=111
export ICLIQ_SW=1
export IOVR=1

# Microphysics
export IMP_PHYSICS=11
export NWAT=6
# GFDL MP
export DNATS=1
export DO_SAT_ADJ=.true.
export LHEATSTRG=.false.
export LSEASPRAY=.false.
export LGFDLMPRAD=.false.
export EFFR_IN=.false.
# Thompson MP
export LRADAR=.true.
export LTAEROSOL=.true.
export EXT_DIAG_THOMPSON=.false.
export sedi_semi=.false.
export sedi_semi_update=.false.
export sedi_semi_decfl=.false.

# GWD
export LDIAG_UGWP=.false.
export DO_UGWP=.false.
export DO_TOFD=.false.
export GWD_OPT=1
export DO_UGWP_V0=.false.
export DO_UGWP_V0_OROG_ONLY=.false.
export DO_GSL_DRAG_LS_BL=.false.
export DO_GSL_DRAG_SS=.false.
export DO_GSL_DRAG_TOFD=.false.
export DO_UGWP_V1=.false.
export DO_UGWP_V1_OROG_ONLY=.false.

# resolution dependent settings
export CDMBWD_c48='0.071,2.1,1.0,1.0'
export CDMBWD_c96='0.14,1.8,1.0,1.0'
export CDMBWD_c192='0.23,1.5,1.0,1.0'
export CDMBWD_c384='1.1,0.72,1.0,1.0'
export CDMBWD_c768='4.0,0.15,1.0,1.0'

# set default
export CDMBWD=${CDMBWD_c96}

# PBL
export SATMEDMF=.false.
export ISATMEDMF=0
export HYBEDMF=.true.
export SHINHONG=.false.
export DO_YSU=.false.
export DO_MYNNEDMF=.false.
export DO_MYJPBL=.false.
export HURR_PBL=.false.
export MONINQ_FAC=1.0

# Shallow/deep convection
export IMFSHALCNV=2
export HWRF_SAMFSHAL=.false.
export IMFDEEPCNV=2
export HWRF_SAMFDEEP=.false.
export RAS=.false.
export RANDOM_CLDS=.false.
export CNVCLD=.true.

# Aerosol convective scavenging
export FSCAV_AERO="'*:0.0'"

# SFC
export DO_MYJSFC=.false.
export DO_MYNNSFCLAY=.false.

# LSM
export LSM=1
export LSOIL_LSM=4
export LANDICE=.true.
export KICE=2
export IALB=1
export IEMS=1

# Ozone / stratospheric H2O
export OZ_PHYS_OLD=.true.
export OZ_PHYS_NEW=.false.
export H2O_PHYS=.false.


export CPL=.false.
export CPLCHM=.false.
export CPLFLX=.false.
export CPLICE=.false.
export CPLWAV=.false.
export CPLWAV2ATM=.false.
export DAYS=1
export NPX=97
export NPY=97
export NPZ=64
export NPZP=65
export NSTF_NAME=2,1,1,0,5
export OUTPUT_FH="12 -1"
export NFHOUT=12
export NFHMAX_HF=12
export NFHOUT_HF=6
export IAU_OFFSET=0
export FHZERO=6
export FNALBC="'global_snowfree_albedo.bosu.t126.384.190.rg.grb'"
export FNVETC="'global_vegtype.igbp.t126.384.190.rg.grb'"
export FNSOTC="'global_soiltype.statsgo.t126.384.190.rg.grb'"
export FNSMCC="'global_soilmgldas.t126.384.190.grb'"
export FNSMCC_control="'global_soilmgldas.statsgo.t1534.3072.1536.grb'"
export FNMSKH_control="'global_slmask.t1534.3072.1536.grb'"
export FNABSC="'global_mxsnoalb.uariz.t126.384.190.rg.grb'"

# Tiled Fix files
export ATMRES=C96
export TILEDFIX=.false.

export ENS_NUM=1
export SYEAR=2016
export SMONTH=10
export SDAY=03
export SHOUR=00
export SECS=`expr $SHOUR \* 3600`
export FHMAX=${FHMAX:-`expr $DAYS \* 24`}
export DT_ATMOS=1800
export FHCYC=24
export FHROT=0
export LDIAG3D=.false.
export QDIAG3D=.false.
export PRINT_DIFF_PGR=.false.
export MAX_OUTPUT_FIELDS=300

# Stochastic physics
export STOCHINI=.false.
export DO_SPPT=.false.
export DO_SHUM=.false.
export DO_SKEB=.false.
export LNDP_TYPE=0
export N_VAR_LNDP=0
export LNDP_EACH_STEP=.false.
export SKEB=-999.
export SPPT=-999.
export SHUM=-999.
export LNDP_VAR_LIST='XXX'
export LNDP_PRT_LIST=-999

#IAU
export IAU_INC_FILES="''"


#Cellular automata
export DO_CA=.false.
export CA_SGS=.false.
export CA_GLOBAL=.false.

export IAU_DRYMASSFIXER=.false.

#waves
export WW3RSTDTHR=12
export DT_2_RST="$(printf "%02d" $(( ${WW3RSTDTHR}*3600 )))"
export DTRST=0
export RSTTYPE=T
export WW3OUTDTHR=1
export DTFLD="$(printf "%02d" $(( ${WW3OUTDTHR}*3600 )))"
export DTPNT="$(printf "%02d" $(( ${WW3OUTDTHR}*3600 )))"
export GOFILETYPE=1
export POFILETYPE=1
export OUTPARS_WAV="WND HS FP DP PHS PTP PDIR"
export CPLILINE='$'
export ICELINE='$'
export WINDLINE='$'
export CURRLINE='$'
export NFGRIDS=0
export NMGRIDS=1
export WW3GRIDLINE="'glo_1deg'  'no' 'no' 'CPL:native' 'no' 'no' 'no' 'no' 'no' 'no'   1  1  0.00 1.00  F"
export FUNIPNT=T
export IOSRV=1
export FPNTPROC=T
export FGRDPROC=T
export UNIPOINTS='points'
export FLAGMASKCOMP=' F'
export FLAGMASKOUT=' F'
export RUN_BEG="${SYEAR}${SMONTH}${SDAY} $(printf "%02d" $(( ${SHOUR}  )))0000"
export RUN_END="2100${SMONTH}${SDAY} $(printf "%02d" $(( ${SHOUR}  )))0000"
export OUT_BEG=$RUN_BEG
export OUT_END=$RUN_END
export RST_BEG=$RUN_BEG
export RST_2_BEG=$RUN_BEG
export RST_END=$RUN_END
export RST_2_END=$RUN_END

# Regional
export WRITE_RESTART_WITH_BCS=.false.

# Diagnostics
export PRINT_DIFF_PGR=.false.

# Coupling
export coupling_interval_fast_sec=0
}

export_cpl ()
{
export FV3=true
export S2S=true
export HAFS=false
export DATM_CDEPS=false
export DOCN_CDEPS=false

export SYEAR=2021
export SMONTH=03
export SDAY=22
export SHOUR=06
export SECS=`expr $SHOUR \* 3600`
export BMIC=.false.

export DAYS=1
export FHMAX=24
export FDIAG=6
export FHZERO=6

# default atm/ocn/ice resolution
export ATMRES=C96
export OCNRES=100
export ICERES=1.00
export NX_GLB=360
export NY_GLB=320
export NPZ=127
export NPZP=128

# default resources
export TASKS=$TASKS_cpl_dflt
export TPN=$TPN_cpl_dflt
export INPES=$INPES_cpl_dflt
export JNPES=$JNPES_cpl_dflt
export THRD=$THRD_cpl_dflt
export WRTTASK_PER_GROUP=$WPG_cpl_dflt

export med_petlist_bounds=$MPB_cpl_dflt
export atm_petlist_bounds=$APB_cpl_dflt
export ocn_petlist_bounds=$OPB_cpl_dflt
export ice_petlist_bounds=$IPB_cpl_dflt
export wav_petlist_bounds=$WPB_cpl_dflt

# component and coupling timesteps
export DT_ATMOS=720
export DT_CICE=${DT_ATMOS}
export DT_DYNAM_MOM6=1800
export DT_THERM_MOM6=3600

# nems.configure defaults
export NEMS_CONFIGURE=nems.configure.cpld_wave.IN
export med_model=cmeps
export atm_model=fv3
export ocn_model=mom6
export ice_model=cice6
export wav_model=ww3

export coupling_interval_slow_sec=${DT_THERM_MOM6}
export coupling_interval_fast_sec=${DT_ATMOS}

export RESTART_N=${FHMAX}
export CPLMODE=nems_frac
export cap_dbug_flag=0
export use_coldstart=false
export use_mommesh=true
export RUNTYPE=startup
export CICERUNTYPE=initial
export eps_imesh=1.0e-1

# FV3 defaults
export FRAC_GRID=.true.
export CCPP_SUITE=FV3_GFS_v16_coupled_nsstNoahmpUGWPv1
export INPUT_NML=cpld_control.nml.IN
export FIELD_TABLE=field_table_gfsv16
export DIAG_TABLE=diag_table_template

export DIAG_TABLE_ADDITIONAL=''

export FHROT=0
export NSOUT=-1
export OUTPUT_FH='6 -1'

#P7 default
export IALB=2
export IEMS=2
export LSM=2
export IOPT_DVEG=4
export IOPT_CRS=2
export IOPT_RAD=3
export IOPT_ALB=1
export IOPT_STC=3

# FV3 P7 settings
export D2_BG_K1=0.20
export D2_BG_K2=0.04
export DZ_MIN=2
export PSM_BC=1
export DDDMP=0.2

# P7 Merra2 Aerosols & NSST
export USE_MERRA2=.true.
export IAER=1011
export NSTF_NAME=2,1,0,0,0

export LHEATSTRG=.true.
export LSEASPRAY=.true.

# P7 UGWP1
export GWD_OPT=2
export DO_UGWP_V1=.true.
export KNOB_UGWP_VERSION=1
export KNOB_UGWP_NSLOPE=1
export DO_UGWP_V0=.false.
export DO_GSL_DRAG_LS_BL=.true.
export DO_GSL_DRAG_SS=.true.
export DO_GSL_DRAG_TOFD=.true.
export DO_UGWP_V1_OROG_ONLY=.false.
export DO_UGWP_V0_NST_ONLY=.false.
export LDIAG_UGWP=.false.

# P7 CA
export DO_CA=.true.
export CA_SGS=.true.
export CA_GLOBAL=.false.
export NCA=1
export NCELLS=5
export NLIVES=12
export NTHRESH=18
export NSEED=1
export NFRACSEED=0.5
export CA_TRIGGER=.true.
export NSPINUP=1
export ISEED_CA=12345

# P7 settings
export FNALBC="'C96.snowfree_albedo.tileX.nc'"
export FNALBC2="'C96.facsf.tileX.nc'"
export FNTG3C="'C96.substrate_temperature.tileX.nc'"
export FNVEGC="'C96.vegetation_greenness.tileX.nc'"
export FNVETC="'C96.vegetation_type.tileX.nc'"
export FNSOTC="'C96.soil_type.tileX.nc'"
export FNSMCC=${FNSMCC_control}
export FNMSKH=${FNMSKH_control}
export FNVMNC="'C96.vegetation_greenness.tileX.nc'"
export FNVMXC="'C96.vegetation_greenness.tileX.nc'"
export FNSLPC="'C96.slope_type.tileX.nc'"
export FNABSC="'C96.maximum_snow_albedo.tileX.nc'"
export LANDICE=".false."
export FSICL=99999
export USE_CICE_ALB=.false.

# P7 default mushy thermo
export KTHERM=2
export TFREEZE_OPTION=mushy

export CPLFLX=.true.
export CPLICE=.true.
export CPL=.true.
export CPLWAV=.true.
export CPLWAV2ATM=.true.
export MIN_SEAICE=1.0e-11

# for FV3: default values will be changed if doing a warm-warm restart
export WARM_START=.false.
export MAKE_NH=.true.
export NA_INIT=1
export EXTERNAL_IC=.true.
export NGGPS_IC=.true.
export MOUNTAIN=.false.

# MOM6 defaults; 1 degree
export MOM_INPUT=MOM_input_template_100
export MOM6_RESTART_SETTING=n
export MOM6_RIVER_RUNOFF=False
export FRUNOFF=''
export CHLCLIM=seawifs_1998-2006_smoothed_2X.nc
export MOM6_USE_LI2016=True
# since CPL_SLOW is set to DT_THERM, this should be always be false
export MOM6_THERMO_SPAN=False
export MOM6_USE_WAVES=True
export MOM6_ALLOW_LANDMASK_CHANGES=False
# MOM6 IAU
export MOM_IAU=False
export MOM_IAU_HRS=6
# MOM6 stochastics
export DO_OCN_SPPT=False
export PERT_EPBL=False
export OCN_SPPT=-999.
export EPBL=-999.

# CICE6 defaults; 1 degree
export CICE_DECOMP=slenderX2
export NPROC_ICE=$NPROC_ICE_cpl_dflt
# SlenderX2
export CICE_DECOMP=slenderX2
export np2=`expr $NPROC_ICE / 2`
export BLCKX=`expr $NX_GLB / $np2`
export BLCKY=`expr $NY_GLB / 2`
export MESHOCN_ICE=mesh.mx${OCNRES}.nc
export CICEGRID=grid_cice_NEMS_mx${OCNRES}.nc
export CICEMASK=kmtu_cice_NEMS_mx${OCNRES}.nc
export RUNID=unknown
# set large; restart frequency now controlled by restart_n in nems.configure
export DUMPFREQ=d
export DUMPFREQ_N=1000
export USE_RESTART_TIME=.false.
export RESTART_EXT=.false.
# setting to true will allow Frazil FW and Salt to be
# included in fluxes sent to ocean
export FRAZIL_FWSALT=.true.
# default to write CICE average history files
export CICE_HIST_AVG=.true.

#wave
export WW3GRIDLINE="'glo_1deg'  'no' 'CPL:native' 'CPL:native' 'CPL:native' 'no' 'no' 'no' 'no' 'no'   1  1  0.00 1.00  F"
export RUN_BEG="${SYEAR}${SMONTH}${SDAY} $(printf "%02d" $(( ${SHOUR}  )))0000"
export RUN_END="2100${SMONTH}${SDAY} $(printf "%02d" $(( ${SHOUR}  )))0000"
export OUT_BEG=$RUN_BEG
export OUT_END=$RUN_END
export RST_BEG=$RUN_BEG
export RST_2_BEG=$RUN_BEG
export RST_END=$RUN_END
export RST_2_END=$RUN_END

# checkpoint restarts
export RESTART_FILE_PREFIX=''
export RESTART_FILE_SUFFIX_HRS=''
export RESTART_FILE_SUFFIX_SECS=''
export RT35D=''
}
export_35d_run ()
{
export CNTL_DIR=""
export LIST_FILES=""
}
export_datm_cdeps ()
{
export FV3=false
export S2S=false
export HAFS=false
export DATM_CDEPS=true
export DOCN_CDEPS=false
export CPLWAV=.false.
export DAYS=1
export FHMAX=24
export THRD=1
export FHROT=0
export WARM_START=.false.

# atm/ocn/ice resolution
export IATM=1760
export JATM=880
export ATM_NX_GLB=$IATM
export ATM_NY_GLB=$JATM
export ATMRES=1760x880
export OCNRES=100
export ICERES=1.00
export NX_GLB=360
export NY_GLB=320

# nems.configure
export NEMS_CONFIGURE=nems.configure.datm_cdeps.IN
export med_model=cmeps
export atm_model=datm
export ocn_model=mom6
export ice_model=cice6
export atm_petlist_bounds=$APB_cdeps_100
export med_petlist_bounds=$MPB_cdeps_100
export ocn_petlist_bounds=$OPB_cdeps_100
export ice_petlist_bounds=$IPB_cdeps_100
export TASKS=$TASKS_cdeps_100
export TPN=$TPN_cdeps_100
# SlenderX2
export CICE_DECOMP=slenderX2
export NPROC_ICE=12
export np2=`expr $NPROC_ICE / 2`
export BLCKX=`expr $NX_GLB / $np2`
export BLCKY=`expr $NY_GLB / 2`

export ENS_NUM=1
export SYEAR=2011
export SMONTH=10
export SDAY=01
export SHOUR=00
export SECS=`expr $SHOUR \* 3600`
export CDATE=${SYEAR}${SMONTH}${SDAY}${SHOUR}

export NFHOUT=6
export DT_ATMOS=900
export DT_CICE=${DT_ATMOS}
export DT_DYNAM_MOM6=1800
export DT_THERM_MOM6=3600
export coupling_interval_slow_sec=${DT_THERM_MOM6}
export coupling_interval_fast_sec=${DT_ATMOS}

export RESTART_N=${FHMAX}
export CPLMODE=nems_orig_data
export cap_dbug_flag=0
export use_coldstart=false
export use_mommesh=true
export RUNTYPE=startup
export CICERUNTYPE=initial
export eps_imesh=1.0e-1
export flux_convergence=0.0
export flux_iteration=2
export flux_scheme=0

export INPUT_NML=input.mom6.nml.IN
export MODEL_CONFIGURE=datm_cdeps_configure.IN
export DIAG_TABLE=diag_table_template

# atm defaults
export DATM_SRC=CFSR
export FILENAME_BASE=cfsr.
export mesh_file=cfsr_mesh.nc
export MESH_ATM=DATM_INPUT/${mesh_file}
export atm_datamode=${DATM_SRC}
export stream_files=DATM_INPUT/${FILENAME_BASE}201110.nc

# MOM6 defaults; 1 degree
export MOM_INPUT=MOM_input_template_100
export MOM6_RESTART_SETTING=n
export MOM6_RIVER_RUNOFF=False
export FRUNOFF=''
export CHLCLIM=seawifs_1998-2006_smoothed_2X.nc
# MOM6 IAU
export MOM_IAU=False
export MOM_IAU_HRS=6
export MOM6_USE_LI2016=False
# MOM6 stochastics
export DO_OCN_SPPT=False
export PERT_EPBL=False
export OCN_SPPT=-999.
export EPBL=-999.
# since coupling_interval_slow is set to DT_THERM, this should be always be false
export MOM6_THERMO_SPAN=False
export MOM6_USE_WAVES=False
export MOM6_ALLOW_LANDMASK_CHANGES=False

# CICE6 defaults; 1 degree
export MESHOCN_ICE=mesh.mx${OCNRES}.nc
export CICEGRID=grid_cice_NEMS_mx${OCNRES}.nc
export CICEMASK=kmtu_cice_NEMS_mx${OCNRES}.nc
export RUNID=unknown
# set large; restart frequency now controlled by restart_n in nems.configure
export DUMPFREQ=d
export DUMPFREQ_N=1000
export USE_RESTART_TIME=.false.
export RESTART_EXT=.false.
# setting to true will allow Frazil FW and Salt to be
# included in fluxes sent to ocean
export FRAZIL_FWSALT=.true.
# default to write CICE average history files
export CICE_HIST_AVG=.true.
# default non-mushy thermo
export KTHERM=1
export TFREEZE_OPTION=linear_salt
export BL_SUFFIX=""
export RT_SUFFIX=""
}
export_hafs_datm_cdeps ()
{
export FV3=false
export S2S=false
export HAFS=true
export DATM_CDEPS=true
export DOCN_CDEPS=false
export THRD=1
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export TASKS=$TASKS_dflt
export TPN=$TPN_dflt

export atm_model=datm

export DATM_IN_CONFIGURE=datm_in
export DATM_STREAM_CONFIGURE=hafs_datm.streams.era5.IN
}
export_hafs_docn_cdeps ()
{
export FV3=true
export S2S=false
export HAFS=true
export DOCN_CDEPS=true
export THRD=1
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export TASKS=$TASKS_dflt
export TPN=$TPN_dflt

export ocn_model=docn
export ocn_datamode=sstdata

export DOCN_IN_CONFIGURE=docn_in
export DOCN_STREAM_CONFIGURE=hafs_docn.streams.IN
}
export_hafs_regional ()
{
export FV3=true
export S2S=false
export HAFS=true
export DATM_CDEPS=false
export DOCN_CDEPS=false
export THRD=1
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export TASKS=$TASKS_dflt
export TPN=$TPN_dflt

# model_configure
export SYEAR=2019
export SMONTH=08
export SDAY=29
export SHOUR=00
export SECS=`expr $SHOUR \* 3600`
export FHMAX=6
export ENS_NUM=1
export DT_ATMOS=900
export CPL=.true.
export RESTART_INTERVAL=0
export FHROT=0
export coupling_interval_fast_sec=0
export QUILTING=.true.
export WRITE_GROUP=1
export WRTTASK_PER_GROUP=6
export OUTPUT_HISTORY=.true.
export WRITE_DOPOST=.false.
export NUM_FILES=2
export FILENAME_BASE="'atm' 'sfc'"
export OUTPUT_GRID="'regional_latlon'"
export OUTPUT_FILE="'netcdf'"
export IDEFLATE=0
export NBITS=0
export NFHOUT=3
export NFHMAX_HF=-1
export NFHOUT_HF=3
export CEN_LON=-62.0
export CEN_LAT=25.0
export LON1=-114.5
export LAT1=-5.0
export LON2=-9.5
export LAT2=55.0
export DLON=0.03
export DLAT=0.03

# input.nml
export CPL_IMP_MRG=.true.

# nems.configure
export med_model=cmeps
export CAP_DBUG_FLAG=0
export RESTART_N=${FHMAX}
export CPLMODE=hafs
export RUNTYPE=startup
export USE_COLDSTART=false
}
