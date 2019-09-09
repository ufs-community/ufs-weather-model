#!/bin/bash

# The list of compilers here must cover all system listed in build_ccpp.sh -> VALID_MACHINES
case "$MACHINE_ID" in
        wcoss_cray)
            export LD=ftn
            export CC=cc
            export CXX=CC
            export FC=ftn
            export F77=ftn
            export F90=ftn
            ;;
        wcoss_dell_p3)
            export CC=mpiicc
            export CXX=mpiicpc
            export FC=mpiifort
            export F77=mpiifort
            export F90=mpiifort
            ;;
        gaea.intel)
            export LD=ftn
            export CC=cc
            export CXX=CC
            export FC=ftn
            export F77=ftn
            export F90=ftn
            ;;
        jet.intel)
            export CC=mpiicc
            export CXX=mpiicpc
            export FC=mpiifort
            export F77=mpiifort
            export F90=mpiifort
            ;;
        hera.intel)
            export CC=mpiicc
            export CXX=mpiicpc
            export FC=mpiifort
            export F77=mpiifort
            export F90=mpiifort
            ;;
        theia.intel)
            export CC=mpiicc
            export CXX=mpiicpc
            export FC=mpiifort
            export F77=mpiifort
            export F90=mpiifort
            ;;
        theia.gnu)
            export CC=mpicc
            export CXX=mpicxx
            export FC=mpif90
            export F77=mpif77
            export F90=mpif90
            ;;
        theia.pgi)
            export CPP="mpicc -E"
            export CC=mpicc
            export CXX=mpic++
            export FC=mpif90
            export F77=mpif77
            export F90=mpif90
            ;;
        cheyenne.intel)
            export CC=mpicc
            export CXX=mpicxx
            export FC=mpif90
            export F77=mpif77
            export F90=mpif90
            ;;
        cheyenne.intel-impi)
            export CC=mpiicc
            export CXX=mpiicpc
            export FC=mpiifort
            export F77=mpiifort
            export F90=mpiifort
            ;;
        cheyenne.gnu)
            export CC=mpicc
            export CXX=mpicxx
            export FC=mpif90
            export F77=mpif77
            export F90=mpif90
            ;;
        cheyenne.pgi)
            export CPP="mpicc -E"
            export CC=mpicc
            export CXX=mpicxx
            export FC=mpif90
            export F77=mpif77
            export F90=mpif90
            ;;
        endeavor.intel)
            export CC=mpiicc
            export CXX=mpiicpc
            export FC=mpiifort
            export F77=mpiifort
            export F90=mpiifort
            ;;
        stampede.intel)
            export CC=mpicc
            export CXX=mpicxx
            export FC=mpif90
            export F77=mpif77
            export F90=mpif90
            ;;
        supermuc_phase2.intel)
            export CC=mpiicc
            export CXX=mpiicpc
            export FC=mpiifort
            export F77=mpiifort
            export F90=mpiifort
            ;;
        macosx.gnu)
            # set in generic modulefile
            ;;
        linux.gnu)
            # set in generic modulefile
            ;;
        *)
            echo "ERROR: MACHINE_ID ${MACHINE_ID} not configured in set_compilers.sh"
            exit 1 
esac

echo "Compilers set for ${MACHINE_ID}."