#!/bin/sh
#SBATCH -e err
#SBATCH -o out
#SBATCH --account=gsienkf
#SBATCH --qos=debug
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=20
#SBATCH --job-name="stoch_unit_tests"

RES=96
NPX=`expr $RES + 1`
NPY=`expr $RES + 1`

source ./module-setup.sh
module purge
module use $( pwd -P )
module load modules.stoch

# compile codes
sh compile_standalone.hera_intel
if [ ! -f standalone_stochy.x ];then
  echo "compilation errors"
  exit 1
fi
sh compile_compare.sh

# copy input directory
if [ ! -d INPUT ]; then
   cp -r /scratch2/BMC/gsienkf/Philip.Pegion/stochastic_physics_unit_tests/input_data INPUT
fi
mkdir -p RESTART

# test 3 different domain decompositions and compare to baseline
##layout 1x4
cp input.nml.template input.nml
sed -i -e "s/LOX/1/g" input.nml
sed -i -e "s/LOY/4/g" input.nml
sed -i -e "s/NPX/$NPX/g" input.nml
sed -i -e "s/NPY/$NPY/g" input.nml
   sed -i -e "s/RES/$RES/g" input.nml
sed -i -e "s/_STOCHINI_/.false./g" input.nml
export OMP_NUM_THREADS=1
time srun --label -n 24 standalone_stochy.x  >& stdout.1x4_layout
mkdir layout_1x4
mv workg* layout_1x4

#layout 2x2
export OMP_NUM_THREADS=2
cp input.nml.template input.nml
sed -i -e "s/LOX/2/g" input.nml
sed -i -e "s/LOY/2/g" input.nml
sed -i -e "s/NPX/$NPX/g" input.nml
sed -i -e "s/NPY/$NPY/g" input.nml
sed -i -e "s/RES/$RES/g" input.nml
sed -i -e "s/_STOCHINI_/.false./g" input.nml
time srun -n 24 standalone_stochy.x >& stdout.2x2_layout
mkdir layout_2x2
mv workg* layout_2x2

#layout 1x4
export OMP_NUM_THREADS=1
cp input.nml.template input.nml
sed -i -e "s/LOX/4/g" input.nml
sed -i -e "s/LOY/1/g" input.nml
sed -i -e "s/NPX/$NPX/g" input.nml
sed -i -e "s/NPY/$NPY/g" input.nml
sed -i -e "s/RES/$RES/g" input.nml
sed -i -e "s/_STOCHINI_/.false./g" input.nml
time srun -n 24 standalone_stochy.x >& stdout.4x1_layout
mkdir layout_4x1
mv workg* layout_4x1
# restart run
mv stochy_middle.nc INPUT/atm_stoch.res.nc
export OMP_NUM_THREADS=2
cp input.nml.template input.nml
sed -i -e "s/LOX/4/g" input.nml
sed -i -e "s/LOY/1/g" input.nml
sed -i -e "s/NPX/$NPX/g" input.nml
sed -i -e "s/NPY/$NPY/g" input.nml
sed -i -e "s/RES/$RES/g" input.nml
sed -i -e "s/_STOCHINI_/.true./g" input.nml
time srun -n 24 standalone_stochy.x >& stdout.2x2_restart
rm workg* 

compare_output
if [ $? -ne 0 ];then
   echo "unit tests failed"
else
   diff stochy_final.nc stochy_final_2.nc
   if [ $? -eq 0 ];then
      echo "unit tests successful"
#      rm -rf layout_*
      rm logfile*
      rm stochy*nc
      rm ../*.o ../*.mod
      rm ../libstochastic_physics.a
      rm standalone_stochy.x
   else
      echo "restart test failed"
   fi
fi
