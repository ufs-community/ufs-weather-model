#!/bin/sh
#SBATCH -e err
#SBATCH -o out
#SBATCH --account=gsienkf
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=20
#SBATCH --job-name="stoch_unit_tests"

RES=96
NPX=`expr $RES + 1`
NPY=`expr $RES + 1`
source ./module-setup.sh
module purge
module use $( pwd -P )
module load modules.stoch_gnu

# compile codes
sh compile_standalone.hera_gnu
if [ ! -f standalone_stochy.x ];then
  echo "compilation errors"
  exit 1
fi

# copy input directory
if [ ! -d INPUT ]; then
   cp -r /scratch2/BMC/gsienkf/Philip.Pegion/stochastic_physics_unit_tests/input_data INPUT
fi
mkdir -p RESTART

# test 3 different domain decompositions and compare to baseline
#layout 1x4
#cp input.nml.template input.nml
sed -i -e "s/LOX/1/g" input.nml
sed -i -e "s/LOY/4/g" input.nml
sed -i -e "s/NPX/$NPX/g" input.nml
sed -i -e "s/NPY/$NPY/g" input.nml
sed -i -e "s/RES/$RES/g" input.nml
sed -i -e "s/_STOCHINI_/.false./g" input.nml
export OMP_NUM_THREADS=2
module list
time srun --label -n 24 standalone_stochy.x 
mkdir stochy_out
mv workg* stochy_out
