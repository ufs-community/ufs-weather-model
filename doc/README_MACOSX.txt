# Dom Heinzeller (dom.heinzeller@noaa.gov), 05/18/2018

In order to build and run the FV3 trunk (May 2018) on Mac OS X, the following installation steps are recommended.
The version numbers for the "brew" correspond to the default versions in April/May 2018 and will change to newer versions
in the future. Unless problems occur during the manual builds in steps 10-12, these differences can be ignored.

1. install homebrew (enter sudo password when requested)
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

2. sudo-create /usr/local/src, /usr/local/esmf-7.1.0r and /usr/local/NCEPlibs-20180401. Change permissions to your user name / user group
    # /usr/local/src
    sudo mkdir /usr/local/src
    sudo chown YOUR_USERNAME /usr/local/src
    sudo chgrp YOUR_GROUPNAME /usr/local/src
    # /usr/local/esmf-7.1.0r
    sudo mkdir /usr/local/esmf-7.1.0r
    sudo chown YOUR_USERNAME /usr/local/esmf-7.1.0r
    sudo chgrp YOUR_GROUPNAME /usr/local/esmf-7.1.0r
    # /usr/local/NCEPlibs-20180401
    sudo mkdir /usr/local/NCEPlibs-20180401
    sudo chown YOUR_USERNAME /usr/local/NCEPlibs-20180401
    sudo chgrp YOUR_GROUPNAME /usr/local/NCEPlibs-20180401

3. Install gcc-7.2.0, gfortran-7.2.0
    brew install gcc --verbose --without-multilib

4. Install clang-5.0.0 with openmp support
    brew install llvm

5. Install mpich-3.2.1
    brew install mpich

6. Install netCDF library
    brew install -v netcdf

7. Install libpng-1.6.34
    brew install libpng

8. Install udunits-2.2.25
    brew install udunits

9. Install ncview-2.1.7
    brew install ncview

10. Install ESMF 7.1.0r

    # Download esmf_7_1_0r_src.tar.gz from https://www.earthsystemcog.org/projects/esmf/download/ to /usr/local/src

    cd /usr/local/src
    tar -xvf esmf_7_1_0r_src.tar.gz
    cd esmf
    export NETCDF=/usr/local
    export ESMF_DIR=`pwd`
    export ESMF_INSTALL_PREFIX=/usr/local/esmf-7.1.0r
    export ESMF_COMPILER=gfortranclang
    export ESMF_CXXCOMPILER=mpicxx
    export ESMF_CXXLINKER=mpicxx
    export ESMF_F90COMPILER=mpif90
    export ESMF_F90LINKER=mpif90
    export ESMF_BOPT=O
    export ESMF_OPTLEVEL=2
    export ESMF_COMM=mpich
    export ESMF_MPIRUN=mpirun
    export ESMF_NETCDF=1
    export ESMF_NETCDF_INCLUDE=$NETCDF/include
    export ESMF_NETCDF_LIBPATH=$NETCDF/lib
    export ESMF_NETCDF=split
    export LDFLAGS="-Wl,-no_compact_unwind"
    make info 2>&1 | tee log.info
    # Ignore warnings of type "ld: warning: could not create compact unwind for ... stack subq instruction is too different from dwarf stack size"
    make 2>&1 | tee log.make
    make check 2>&1 | tee log.check
    # Found 40 multi-processor system tests, 40 passed and 0 failed.
    # Found 3190 non-exhaustive multi-processor unit tests, 2790 passed and 400 failed
    # --> ignore and proceed
    make install 2>&1 | tee log.install
    make installcheck 2>&1 | tee log.installcheck

    cd $ESMF_INSTALL_PREFIX
    #
    # Fix wrong path to libesmf.dylib in ESMF binaries - this will hopefully be addressed in future ESMF releases
    #
    install_name_tool -change $ESMF_DIR/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib $ESMF_INSTALL_PREFIX/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib bin/binO/Darwin.gfortranclang.64.mpich.default/ESMF_Info
    install_name_tool -change $ESMF_DIR/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib $ESMF_INSTALL_PREFIX/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib bin/binO/Darwin.gfortranclang.64.mpich.default/ESMF_InfoC
    install_name_tool -change $ESMF_DIR/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib $ESMF_INSTALL_PREFIX/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib bin/binO/Darwin.gfortranclang.64.mpich.default/ESMF_WebServController
    install_name_tool -change $ESMF_DIR/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib $ESMF_INSTALL_PREFIX/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib bin/binO/Darwin.gfortranclang.64.mpich.default/ESMF_Regrid
    install_name_tool -change $ESMF_DIR/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib $ESMF_INSTALL_PREFIX/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib bin/binO/Darwin.gfortranclang.64.mpich.default/ESMF_RegridWeightGen
    install_name_tool -change $ESMF_DIR/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib $ESMF_INSTALL_PREFIX/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib bin/binO/Darwin.gfortranclang.64.mpich.default/ESMF_Scrip2Unstruct
    #
    # Fix wrong ID in libesmf.dylib - this will hopefully be addressed in future ESMF releases
    #
    install_name_tool -id $ESMF_INSTALL_PREFIX/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib

    # Clean up
    cd /usr/local/src
    rm -fr esmf
    export -n ESMF_DIR
    export -n ESMF_INSTALL_PREFIX
    export -n ESMF_COMPILER
    export -n ESMF_CXXCOMPILER
    export -n ESMF_CXXLINKER
    export -n ESMF_F90COMPILER
    export -n ESMF_F90LINKER
    export -n ESMF_BOPT
    export -n ESMF_OPTLEVEL
    export -n ESMF_COMM
    export -n ESMF_MPIRUN
    export -n ESMF_NETCDF
    export -n ESMF_NETCDF_INCLUDE
    export -n ESMF_NETCDF_LIBPATH
    export -n ESMF_NETCDF
    export -n NETCDF
    export -n LDFLAGS

11. Build external NCEP libraries (use date tag 20180401 to allow for different versions in the future)

    # Obtain source code from gitub and build in /usr/local/src
    cd /usr/local/src
    git clone https://github.com/climbfuji/NCEPlibs.git
    cd NCEPlibs
    ./make_ncep_libs.sh -s macosx -c gnu -d /usr/local/NCEPlibs-20180401 -o 1 2>&1 | tee log.make

12. Build model. Change to top-level directory of NEMSfv3gfs.

    # Set environment variables for MACOSX with clang/gfortran
    . modulefiles/macosx.gnu/fv3

    # Build model
    cd tests
    ./compile.sh $PWD/../FV3 macosx.gnu 2>&1 | tee log.compile

13. Set up the run directory using the template on Theia or Cheyenne at some location on your machine:

    a) copy the contents of the run directory templates to where you want to run the model, change to this directory
       (these folders are read-only, i.e. users might have to add the write-flag after copying/rsyncing them)

        theia:    /scratch4/BMC/gmtb/Dom.Heinzeller/macosx_rundirs/C96_trunk_20180427/gnu/
        cheyenne: /glade/p/work/heinzell/fv3/macosx_rundirs/C96_trunk_20180427/gnu/

    b) edit run_macosx.sh in change the variable FV3_BUILD_DIR to the top-level directory of your FV3-build

    c) source the setenv_develop.sh script and execute the model run using the wrapper run_macosx.sh
        ./run_macosx.sh 2>&1 | tee run_macosx.log
        # or, with X OpenMP threads
        OMP_NUM_THREADS=X ./run_macosx.sh 2>&1 | tee run_macosx.log

    d) go and get yourself a cup of coffee ...
