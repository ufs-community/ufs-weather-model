#!/bin/ksh
set -x
set -e

# DEBUG switch on (1) or off (0)
CLEAN=1
# DEBUG switch on (1) or off (0)
DEBUG=0
# OPENMP switch on (1) or off (0)
OPENMP=1

# List of valid/tested machines
valid_machines=( theia cheyenne macosx )
valid_compilers=( intel pgi gnu )

function usage   {
  echo "Usage: "
  echo "build.sh machine compiler homedir"
  echo "    Where: machine  [required] can be : ${valid_machines[@]}"
  echo "           compiler [required] can be : ${valid_compilers[@]}"
  echo "           homedir  [optional] can be any valid directory with write permissions"
  exit 1
}

if [[ $1 = "help" ]] ; then usage; fi

# Always specify host and compiler
if [[ $# -lt 2 ]];  then usage; fi
machine=${1}
compiler=${2}
arch=${machine}.${compiler}
homedir=${3:-`pwd`/../../..}
idir=$(realpath ${homedir})

# Build the various FV3 binaries
cd $homedir/tests
# Set debug flag
if [ "$DEBUG" -eq 1 ]; then
  debug_compile_option="DEBUG=Y"
else
  debug_compile_option="DEBUG=N"
fi
# Set OpenMP flag
if [ "$OPENMP" -eq 1 ]; then
  openmp_compile_option="OPENMP=Y"
else
  openmp_compile_option="OPENMP=N"
fi
# Set precision
if [[ $machine == "theia" ]]; then
  precision_option="32BIT=Y"
  precision="32bit"
elif [[ $machine == "cheyenne" ]]; then
  precision_option="32BIT=N"
  precision="64bit"
elif [[ $machine == "macosx" ]]; then
  precision_option="32BIT=Y"
  precision="32bit"
fi
# Set hydrostatic/non-hydrostatic
hydro_option="HYDRO=N"
#hydro_option="HYDRO=Y"

compile_option="$debug_compile_option $openmp_compile_option $hydro_option $precision_option"
./compile.sh $homedir/FV3 $arch "$compile_option" 1
cp $homedir/tests/fv3_1.exe ../NEMS/exe/fv3_gfs_nh.prod.${precision}.${compiler}.x
rm $homedir/tests/fv3_1.exe
