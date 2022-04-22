#!/bin/sh
compile_all=1
DEBUG=YES
#DEBUG=NO
source ./module-setup.sh
module purge
module use $( pwd -P )
if [ $DEBUG == 'YES' ]; then
   module load modules.stoch_gnu_dbg
else
   module load modules.stoch_gnu
fi
#module list
rm standalone_ca.x
FC=mpif90
FMS_INC=${FMS_ROOT}/include_r4
FMS_LIB=${FMS_ROOT}/lib
INCS="-I. -I${FMS_INC} -I${NETCDF}/include"
if [ $DEBUG == 'YES' ]; then
   FLAGS="-DDEBUG -ggdb -fbacktrace -cpp -fcray-pointer -ffree-line-length-none -fno-range-check -fdefault-real-8 -g -O0 -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check -fopenmp -c "$INCS
   FLAGS2=$FLAGS
else
   FLAGS="-fbacktrace -cpp -fcray-pointer -ffree-line-length-none -fno-range-check -fdefault-real-8 -fdefault-double-8 -g -O2 -fopenmp -c "$INCS
   FLAGS2=$FLAGS
fi
cd ..
if [ $compile_all -eq 1 ];then
   rm -f *.i90 *.i *.o *.mod lib*a
   $FC ${FLAGS} kinddef.F90
   $FC ${FLAGS} mpi_wrapper.F90
   $FC ${FLAGS2} unit_tests/fv_arrays_stub.F90
   $FC ${FLAGS2} unit_tests/fv_mp_stub_mod.F90
   $FC ${FLAGS2} unit_tests/fv_control_stub.F90
   $FC ${FLAGS2} unit_tests/atmosphere_stub.F90
   $FC ${FLAGS2} random_numbers.F90
   $FC ${FLAGS} halo_exchange.fv3.F90
   $FC ${FLAGS} mersenne_twister.F90
   $FC ${FLAGS} plumes.F90 
   $FC ${FLAGS} update_ca.F90
   $FC ${FLAGS} cellular_automata_sgs.F90
   $FC ${FLAGS} cellular_automata_global.F90
   ar rv libcellular_automata.a *.o
fi
   $FC ${FLAGS} update_ca.F90
exit
if [ $DEBUG == 'YES' ]; then
   $FC -fdec -ggdb -fbacktrace -cpp -fcray-pointer -ffree-line-length-none -fno-range-check -fdefault-real-8 -fdefault-double-8 -g -O0 -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check -I. -fopenmp -o unit_tests/standalone_ca.x unit_tests/standalone_ca.F90 ${INCS} -I${NETCDF}/include -L. -lcellular_automata -L${FMS_LIB} -lfms_r4 -L${ESMF_LIB} -Wl,-rpath,${ESMF_LIB} -lesmf -L${NETCDF}/lib -lnetcdff -lnetcdf -L${HDF5_LIBRARIES} -lhdf5_hl -lhdf5 \
-L${ZLIB_LIBRARIES} -lz -ldl
else
   $FC -fdec -fbacktrace -cpp -fcray-pointer -ffree-line-length-none -fno-range-check -fdefault-real-8 -fdefault-double-8 -g -O2 -I. -fopenmp -o unit_tests/standalone_ca.x unit_tests/standalone_ca.F90 ${INCS} -I${NETCDF}/include -L. -lcellular_automata -L${FMS_LIB} -lfms_r4 -L${ESMF_LIB} -Wl,-rpath,${ESMF_LIB} -lesmf -L${NETCDF}/lib -lnetcdff -lnetcdf -L${HDF5_LIBRARIES} -lhdf5_hl -lhdf5 \
-L${ZLIB_LIBRARIES} -lz -ldl
fi

