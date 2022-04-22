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
DO_CA_SGS=.true.
DO_CA_GLOBAL=.true.

source ./module-setup.sh
module purge
module use $( pwd -P )
module load modules.stoch
EXEC=standalone_ca.x
# compile codes
sh compile_standalone_ca.hera.intel
if [ ! -f $EXEC ];then
  echo "compilation errors"
  exit 1
fi
#sh compile_compare_ca.sh

# copy input directory
if [ ! -d INPUT ]; then
   cp -r /scratch2/BMC/gsienkf/Philip.Pegion/stochastic_physics_unit_tests/input_data INPUT
fi
mkdir -p RESTART

   #layout 1x1
   cp input.nml.ca_template input.nml
   sed -i -e "s/LOX/1/g" input.nml
   sed -i -e "s/LOY/1/g" input.nml
   sed -i -e "s/NPX/$NPX/g" input.nml
   sed -i -e "s/NPY/$NPY/g" input.nml
   sed -i -e "s/RES/$RES/g" input.nml
   sed -i -e "s/CA_SGS/${DO_CA_SGS}/g" input.nml
   sed -i -e "s/CA_GLOBAL/${DO_CA_GLOBAL}/g" input.nml
   sed -i -e "s/WARM_START/.false./g" input.nml
   export OMP_NUM_THREADS=1
   time srun --label -n 6 $EXEC >& stdout.1x1
   mkdir ca_layout_1x1
   mv ca_out* ca_layout_1x1
   ct=1
   while [ $ct -le 6 ];do
      mv RESTART/mid_run.ca_data.tile${ct}.nc INPUT/ca_data.tile${ct}.nc
      mv RESTART/ca_data.tile${ct}.nc RESTART/run1_end_ca_data.tile${ct}.nc
      ct=`expr $ct + 1`
   done
      
   cp input.nml.ca_template input.nml
   sed -i -e "s/LOX/1/g" input.nml
   sed -i -e "s/LOY/1/g" input.nml
   sed -i -e "s/NPX/$NPX/g" input.nml
   sed -i -e "s/NPY/$NPY/g" input.nml
   sed -i -e "s/RES/$RES/g" input.nml
   sed -i -e "s/CA_SGS/${DO_CA_SGS}/g" input.nml
   sed -i -e "s/CA_GLOBAL/${DO_CA_GLOBAL}/g" input.nml
   sed -i -e "s/WARM_START/.true./g" input.nml
   time srun --label -n 6 $EXEC >& stdout.1x1_restart
   mkdir ca_layout_1x1_restart
   mv ca_out* ca_layout_1x1_restart
   ct=1
   while [ $ct -le 6 ];do
      diff RESTART/ca_data.tile${ct}.nc RESTART/run1_end_ca_data.tile${ct}.nc
      if [ $? -ne 0 ];then
         echo "restart test failed"
         exit 1
      fi   
      ct=`expr $ct + 1`
   done
   if [ $? -eq 0 ];then
      echo "unit test 1 successful"
   fi
   set OMP_NUM_THREADS=2
      
   cp input.nml.ca_template input.nml
   ct=1
   while [ $ct -le 6 ];do
      rm INPUT/ca_data.tile${ct}.nc
      ct=`expr $ct + 1`
   done
   sed -i -e "s/LOX/1/g" input.nml
   sed -i -e "s/LOY/1/g" input.nml
   sed -i -e "s/NPX/$NPX/g" input.nml
   sed -i -e "s/NPY/$NPY/g" input.nml
   sed -i -e "s/RES/$RES/g" input.nml
   sed -i -e "s/CA_SGS/${DO_CA_SGS}/g" input.nml
   sed -i -e "s/CA_GLOBAL/${DO_CA_GLOBAL}/g" input.nml
   sed -i -e "s/WARM_START/.true./g" input.nml
   time srun --label -n 6 $EXEC >& stdout.1x1_restart2
   mkdir ca_layout_1x1_restart2
   mv ca_out* ca_layout_1x1_restart2
   ct=1
   while [ $ct -le 6 ];do
      diff RESTART/ca_data.tile${ct}.nc RESTART/run1_end_ca_data.tile${ct}.nc
      if [ $? -ne 0 ];then
         echo "restart test failed"
         exit 1
      fi   
      ct=`expr $ct + 1`
   done
   if [ $? -eq 0 ];then
      echo "unit test 1 successful"
   fi
   set OMP_NUM_THREADS=2
   cp input.nml.ca_template input.nml
   sed -i -e "s/LOX/1/g" input.nml
   sed -i -e "s/LOY/4/g" input.nml
   sed -i -e "s/NPX/$NPX/g" input.nml
   sed -i -e "s/NPY/$NPY/g" input.nml
   sed -i -e "s/RES/$RES/g" input.nml
   sed -i -e "s/CA_SGS/${DO_CA_SGS}/g" input.nml
   sed -i -e "s/CA_GLOBAL/${DO_CA_GLOBAL}/g" input.nml
   sed -i -e "s/WARM_START/.true./g" input.nml
   time srun --label -n 24 $EXEC >& stdout.1x4_restart
   mkdir ca_layout_1x4_restart
   mv ca_out* ca_layout_1x4_restart
   ct=1
   while [ $ct -le 6 ];do
      diff RESTART/ca_data.tile${ct}.nc RESTART/run1_end_ca_data.tile${ct}.nc
      if [ $? -ne 0 ];then
         echo "restart test failed"
         exit 1
      fi   
      ct=`expr $ct + 1`
   done
   if [ $? -eq 0 ];then
      echo "unit test 2 successful"
   fi
   cp input.nml.ca_template input.nml
   sed -i -e "s/LOX/1/g" input.nml
   sed -i -e "s/LOY/4/g" input.nml
   sed -i -e "s/NPX/$NPX/g" input.nml
   sed -i -e "s/NPY/$NPY/g" input.nml
   sed -i -e "s/RES/$RES/g" input.nml
   sed -i -e "s/CA_SGS/${DO_CA_SGS}/g" input.nml
   sed -i -e "s/CA_GLOBAL/${DO_CA_GLOBAL}/g" input.nml
   sed -i -e "s/WARM_START/.false./g" input.nml
   time srun --label -n 24 $EXEC >& stdout.1x4
   mkdir ca_layout_1x4
   mv ca_out* ca_layout_1x4
   ct=1
   while [ $ct -le 6 ];do
      diff RESTART/ca_data.tile${ct}.nc RESTART/run1_end_ca_data.tile${ct}.nc
      if [ $? -ne 0 ];then
         echo "restart test failed"
         exit 1
      fi   
      ct=`expr $ct + 1`
   done
   if [ $? -eq 0 ];then
      echo "unit test 3 successful"
   fi
   diff ca_layout_1x4/ca_out.tile12.nc intel/ca_layout_1x4
   if [ $? -ne 0 ];then
      echo "unit test 4 successful"
   else
      echo "unit test 4 failed"
   fi
