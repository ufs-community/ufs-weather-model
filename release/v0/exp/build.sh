#!/bin/ksh

set -e

# DEBUG switch on (1) or off (0)
CLEAN=1
# DEBUG switch on (1) or off (0)
DEBUG=0
# OPENMP switch on (1) or off (0)
OPENMP=1

# List of valid/tested machines
valid_machines=(hera theia cheyenne macosx linux wcoss_cray wcoss_phase1 wcoss_phase2 jet gaea)
valid_compilers=(intel pgi gnu)

function usage   {
  echo " "
  echo "Usage: "
  echo "build.sh machine compiler homedir"
  echo "    Where: machine  [required] can be : ${valid_machines[@]}"
  echo "           compiler [required] can be : ${valid_compilers[@]}"
  echo "                                        (wcoss, jet, gaea: intel only)"
  echo "           homedir  [optional] can be any valid directory with write permissions"
  echo " "
  echo "Further compile options are set at the top of build.sh: CLEAN, DEBUG, OPENMP"
  echo " "
  exit 1
}

if [[ $1 = "help" ]] ; then usage; fi

# Always specify host and compiler
if [[ $# -lt 2 ]];  then usage; fi
machine=${1}
compiler=${2}
if [[ ${machine} == hera || ${machine} == theia || ${machine} == cheyenne || ${machine} == macosx || ${machine} == linux ]]; then
  arch=${machine}.${compiler}
elif [[ ${machine} == wcoss_cray || ${machine} == wcoss_phase1 || ${machine} == wcoss_phase2 || ${machine} == jet || ${machine} == gaea ]]; then
  if [[ ${compiler} == intel ]]; then
    arch=${machine}
  else
    usage
  fi
fi

set -x

homedir=${3:-`pwd`/../../..}

# Build the various FV3 binaries
cd $homedir/tests
# Set debug flag
if [ "$DEBUG" -eq 1 ]; then
  debug_compile_option="DEBUG=Y"
  mode="debug"
else
  debug_compile_option="DEBUG=N"
  mode="prod"
fi
# Set OpenMP flag
if [ "$OPENMP" -eq 1 ]; then
  openmp_compile_option="OPENMP=Y"
else
  openmp_compile_option="OPENMP=N"
fi

# 32-bit non-hydrostatic
precision_option="32BIT=Y"
precision="32bit"
hydro_option="HYDRO=N"
hydro="nh"
compile_option="$debug_compile_option $openmp_compile_option $hydro_option $precision_option"
./compile.sh $homedir/FV3 $arch "$compile_option" 1
cp $homedir/tests/fv3_1.exe ../NEMS/exe/fv3_gfs_${hydro}.${mode}.${precision}.${compiler}.x
rm $homedir/tests/fv3_1.exe

# 32-bit hydrostatic
#precision_option="32BIT=Y"
#precision="32bit"
#hydro_option="HYDRO=Y"
#hydro="hydro"
#compile_option="$debug_compile_option $openmp_compile_option $hydro_option $precision_option"
#./compile.sh $homedir/FV3 $arch "$compile_option" 1
#cp $homedir/tests/fv3_1.exe ../NEMS/exe/fv3_gfs_${hydro}.${mode}.${precision}.${compiler}.x
#rm $homedir/tests/fv3_1.exe

# 64-bit non-hydrostatic
#precision_option="32BIT=N"
#precision="64bit"
#hydro_option="HYDRO=N"
#hydro="nh"
#compile_option="$debug_compile_option $openmp_compile_option $hydro_option $precision_option"
#./compile.sh $homedir/FV3 $arch "$compile_option" 1
#cp $homedir/tests/fv3_1.exe ../NEMS/exe/fv3_gfs_${hydro}.${mode}.${precision}.${compiler}.x
#rm $homedir/tests/fv3_1.exe

# 64-bit hydrostatic
#precision_option="32BIT=N"
#precision="64bit"
##hydro_option="HYDRO=Y"
#hydro="hydro"
#compile_option="$debug_compile_option $openmp_compile_option $hydro_option $precision_option"
#./compile.sh $homedir/FV3 $arch "$compile_option" 1
#cp $homedir/tests/fv3_1.exe ../NEMS/exe/fv3_gfs_${hydro}.${mode}.${precision}.${compiler}.x
#rm $homedir/tests/fv3_1.exe
