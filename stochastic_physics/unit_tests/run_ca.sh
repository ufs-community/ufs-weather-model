#!/bin/sh
#SBATCH -e err
#SBATCH -o out
#SBATCH --account=gsienkf
#SBATCH --qos=debug
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=40
#SBATCH --time=20
#SBATCH --job-name="stoch_unit_tests"

source ./module-setup.sh
module purge
module use $( pwd -P )
module load modules.hera.intel
EXEC=standalone_ca.x

ulimit -s unlimited
export OMP_STACKSIZE=512M
export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=1

cp input.nml.noise input.nml
sed -i -e "s/NOISE/0/g" input.nml
echo "option 0 run 1"
sleep 5
time srun --label -n 384 $EXEC  >& stdout_option_0
mkdir option_0
mv ca_out* option_0

cp input.nml.noise input.nml
sed -i -e "s/NOISE/1/g" input.nml
echo "option 1 run 1"
sleep 5
time srun --label -n 384 $EXEC  >& stdout_option_1
mkdir option_1
mv ca_out* option_1

cp input.nml.noise input.nml
sed -i -e "s/NOISE/2/g" input.nml
echo "option 2 run 1"
sleep 5
time srun --label -n 384 $EXEC  >& stdout_option_2
mkdir option_2
mv ca_out* option_2
exit

cp input.nml.noise input.nml
sed -i -e "s/NOISE/2/g" input.nml
echo "option 2 run 2"
sleep 5
time srun --label -n 384 $EXEC  >& stdout_option_2b
mkdir option_2b
mv ca_out* option_2b

cp input.nml.noise input.nml
sed -i -e "s/NOISE/1/g" input.nml
echo "option 1 run 2"
sleep 5
time srun --label -n 384 $EXEC  >& stdout_option_1b
mkdir option_1b
mv ca_out* option_1b

cp input.nml.noise input.nml
sed -i -e "s/NOISE/0/g" input.nml
echo "option 0 run 2"
sleep 5
time srun --label -n 384 $EXEC  >& stdout_option_0b
mkdir option_0b
mv ca_out* option_0b
