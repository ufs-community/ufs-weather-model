# Dom Heinzeller (dom.heinzeller@noaa.gov), 08/12/2019

Target systems: macOS High Sierra / macOS Mojave with LLVM clang + GNU gfortran compilers and mpich MPI library.

In order to build and run the FV3 trunk (August 2019) with possible CCPP extensions by GMTB on macOS,
the following installation steps are recommended. The version numbers for the "brew" correspond to the default versions
in August 2019 and will change to newer versions in the future. Unless problems occur during the manual builds in
steps 6-11, these differences can be ignored. It is also assumed that the bash shell is used in the following.

1. Install homebrew (enter sudo password when requested)
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

2. Create /usr/local/src, /usr/local/esmf-8.0.0_bs40 and /usr/local/NCEPlibs-20190811. Change permissions to your user name / user group

    # change to root mode
    sudo su
    # /usr/local/src
    mkdir /usr/local/src
    chown YOUR_USERNAME /usr/local/src
    chgrp YOUR_GROUPNAME /usr/local/src
    # /usr/local/esmf-8.0.0_bs40
    mkdir /usr/local/esmf-8.0.0_bs40
    chown YOUR_USERNAME /usr/local/esmf-8.0.0_bs40
    chgrp YOUR_GROUPNAME /usr/local/esmf-8.0.0_bs40
    # /usr/local/NCEPlibs-20190811
    mkdir /usr/local/NCEPlibs-20190811
    chown YOUR_USERNAME /usr/local/NCEPlibs-20190811
    chgrp YOUR_GROUPNAME /usr/local/NCEPlibs-20190811
    # leave root mode
    exit

3. Use homebrew to install the following packages

    a) Install gcc-9.1.0, gfortran-9.1.0
    brew install gcc@9

    b) Install clang-8.0.1 with openmp support
    brew install llvm

    c) Install mpich-3.3.1
    brew install mpich

    d) Install netcdf-4.6.3
    brew install netcdf

    e) Install libpng-1.6.37
    brew install libpng

    f) Install udunits-2.2.27
    brew install udunits

    g) Install cmake-3.15.2
    brew install cmake

    h) Install xquartz-2.7.11 (optional, only needed for ncview)
    brew cask install xquartz

    i) Install ncview-2.1.7 (optional)
    brew install ncview

4. Create a shell setup script ~/setenv_develop_nemsfv3gfs.sh to set the required paths for compiling NEMSfv3gfs. Contents:

####################################################################################
#!/bin/bash

echo "Setting environment variables for develop-NEMSfv3gfs"

export PATH="/usr/local/opt/llvm/bin:$PATH"
export CPPFLAGS="-I/usr/local/opt/llvm/include"
export LDFLAGS="-L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib"
export LIBS_OPENMP="-L/usr/local/opt/llvm/lib -lomp"

export CC=/usr/local/opt/llvm/bin/clang
export CXX=/usr/local/opt/llvm/bin/clang++
export FC=/usr/local/bin/gfortran
export F90=/usr/local/bin/gfortran
export F77=/usr/local/bin/gfortran

export MPICC="mpicc -cc=${CC}"
export MPICXX="mpicxx -cxx=${CXX}"
export MPIFORT=/usr/local/bin/mpifort
export MPIF77=/usr/local/bin/mpif77
export MPIF90=/usr/local/bin/mpif90

export HDF5=/usr/local
export NETCDF=/usr/local
export ESMFMKFILE=/usr/local/esmf-8.0.0_bs40/lib/esmf.mk
export NCEPLIBS_DIR=/usr/local/NCEPlibs-20190811
export MKL_DIR=/opt/intel/compilers_and_libraries_2019.4.233/mac/mkl

####################################################################################

5. Source this shell script for compiling the libraries and the model below
    . ~/setenv_develop_nemsfv3gfs.sh

6. Install ESMF 8.0.0_bs40

    # Download ESMF_8_0_0_beta_snapshot_40 from https://sourceforge.net/p/esmf/esmf/ref/master/tags/
    # to /usr/local/src (creates a directory esmf-esmf-...), rename it to esmf-8.0.0_bs40 and tar it
    # up for later use
    cd /usr/local/src
    mv esmf-esmf-... esmf-8.0.0_bs40
    tar -cvzf esmf-8.0.0_bs40.tar.gz esmf-8.0.0_bs40

    # Compile and install ESMF
    cd esmf-8.0.0_bs40
    export ESMF_DIR=`pwd`
    export ESMF_COMPILER=gfortranclang
    export ESMF_CXXCOMPILER=$MPICXX
    export ESMF_CXXLINKER=$MPICXX
    export ESMF_F90COMPILER=$MPIF90
    export ESMF_F90LINKER=$MPIF90
    export ESMF_BOPT=O
    export ESMF_OPTLEVEL=2
    export ESMF_COMM=mpich
    export ESMF_MPIRUN=mpirun
    export ESMF_NETCDF=1
    export ESMF_NETCDF_INCLUDE=$NETCDF/include
    export ESMF_NETCDF_LIBPATH=$NETCDF/lib
    export ESMF_NETCDF_LIBS="-lnetcdff -lnetcdf -lmpichf90"
    export ESMF_NETCDF=split
    export ESMF_INSTALL_PREFIX=/usr/local/esmf-8.0.0_bs40
    export ESMF_INSTALL_BINDIR=bin
    export ESMF_INSTALL_LIBDIR=lib
    export ESMF_INSTALL_MODDIR=mod
    #
    make info 2>&1 | tee log.info
    make 2>&1 | tee log.make
    # "make check" is optional and can take very long time
    make check 2>&1 | tee log.check
    # SYSTEM TESTS SUMMARY
    # Found 45 multi-processor system tests, 38 passed and 7 failed.
    # UNIT TESTS SUMMARY
    # Found 3466 non-exhaustive multi-processor unit tests, 3394 passed and 72 failed.
    # --> ignore and proceed
    make install 2>&1 | tee log.install
    make installcheck 2>&1 | tee log.installcheck
    #
    cd $ESMF_INSTALL_PREFIX
    #
    # Fix wrong path to libesmf.dylib in ESMF binaries - this will hopefully be addressed in future ESMF releases
    #
    install_name_tool -change $ESMF_DIR/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib $ESMF_INSTALL_PREFIX/lib/libesmf.dylib bin/ESMF_Info
    install_name_tool -change $ESMF_DIR/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib $ESMF_INSTALL_PREFIX/lib/libesmf.dylib bin/ESMF_InfoC
    install_name_tool -change $ESMF_DIR/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib $ESMF_INSTALL_PREFIX/lib/libesmf.dylib bin/ESMF_WebServController
    install_name_tool -change $ESMF_DIR/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib $ESMF_INSTALL_PREFIX/lib/libesmf.dylib bin/ESMF_Regrid
    install_name_tool -change $ESMF_DIR/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib $ESMF_INSTALL_PREFIX/lib/libesmf.dylib bin/ESMF_RegridWeightGen
    install_name_tool -change $ESMF_DIR/lib/libO/Darwin.gfortranclang.64.mpich.default/libesmf.dylib $ESMF_INSTALL_PREFIX/lib/libesmf.dylib bin/ESMF_Scrip2Unstruct
    #
    # Fix wrong ID in libesmf.dylib - this will hopefully be addressed in future ESMF releases
    #
    install_name_tool -id $ESMF_INSTALL_PREFIX/lib/libesmf.dylib lib/libesmf.dylib
    install_name_tool -id $ESMF_INSTALL_PREFIX/lib/libesmf_fullylinked.dylib lib/libesmf_fullylinked.dylib
    #
    # Clean up
    cd /usr/local/src
    rm -fr esmf-8.0.0_bs40
    export -n ESMF_DIR
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
    export -n ESMF_NETCDF_LIBS
    export -n ESMF_NETCDF
    export -n ESMF_INSTALL_PREFIX
    export -n ESMF_INSTALL_BINDIR
    export -n ESMF_INSTALL_LIBDIR
    export -n ESMF_INSTALL_MODDIR

7. Build external NCEP libraries (use date tag 20190811 to allow for different versions in the future)

    # Obtain source code from gitub and build in /usr/local/src
    cd /usr/local/src
    git clone https://github.com/NCAR/NCEPlibs.git NCEPlibs-20190811
    cd NCEPlibs-20190811
    # Requires exporting CC, F90, MPIF90 (done by setenv_develop_nemsfv3gfs.sh)
    ./make_ncep_libs.sh -s macosx -c gnu -d /usr/local/NCEPlibs-20190811 -o 1 2>&1 | tee log.make

8. Download and install Intel Math Kernel Library MKL (full package) from https://software.intel.com/en-us/mkl
   to /opt/intel (default location) using the installer script (requires sudo/root access)

9. Pro-tip (optional): have your computer remember your GitHub username and password.
   See https://help.github.com/en/articles/caching-your-github-password-in-git.

    git config --global credential.helper osxkeychain

10. Download the model from GitHub:
    cd $HOME
    mkdir NEMSfv3gfs
    cd NEMSfv3gfs
    git clone --branch=gmtb/develop https://github.com/NCAR/NEMSfv3gfs NEMSfv3gfs-gmtb-develop-20190811
    cd NEMSfv3gfs-gmtb-develop-20190811
    git submodule init
    git submodule update

11. Build model. Change to top-level directory of NEMSfv3gfs-gmtb-develop-20190811

    . ~/setenv_develop_nemsfv3gfs.sh
    cd tests
    # Note: omit '32BIT=Y' to compile dynamics in double precision (slower to run)
    ./compile.sh $PWD/../FV3 macosx.gnu '32BIT=Y CCPP=N' 2>&1 | tee log.compile                                     # without CCPP
    ./compile.sh $PWD/../FV3 macosx.gnu '32BIT=Y CCPP=Y' 2>&1 | tee log.compile                                     # with CCPP, dynamic mode
    ./compile.sh $PWD/../FV3 macosx.gnu '32BIT=Y CCPP=Y STATIC=Y SUITES=FV3_GFS_2017_gfdlmp' 2>&1 | tee log.compile # with CCPP, static mode, GFS suite

12. Set up the run directory using the template on Theia or Cheyenne at some location on your machine:

    a) copy the contents of the run directory templates to where you want to run the model, change to this directory
       (these folders are read-only, i.e. users might have to add the write-flag after copying/rsyncing them)

        theia:    /scratch4/BMC/gmtb/Dom.Heinzeller/rundirs/20190811/macosx/fv3_gfdlmp/
        cheyenne: /glade/p/ral/jntp/GMTB/NEMSfv3gfs/rundirs/20190811/macosx/fv3_gfdlmp/

    b) edit run_macosx.sh, set variables and change the variable FV3_BUILD_DIR to the top-level directory of your FV3-build

    c) source ~/setenv_develop_nemsfv3gfs.sh and execute the model run using the wrapper run_macosx.sh
        . ~/setenv_develop_nemsfv3gfs.sh
        ./run_macosx.sh 2>&1 | tee run_macosx.log
        # or, with N OpenMP threads (use N=1 for the dynamic CCPP build)
        OMP_NUM_THREADS=N ./run_macosx.sh 2>&1 | tee run_macosx.log

    d) go and get yourself a cup of coffee ...
