# Dom Heinzeller (dom.heinzeller@noaa.gov), 08/21/2019

In order to build and run the FV3 trunk (August 2019) with possible CCPP extensions by GMTB on CentOS Linux,
the following installation steps are recommended. The version numbers correspond to the default versions in
August 2019 and will change to newer versions in the future. Unless problems occur during the manual builds in
step 4, these differences can be ignored. It is also assumed that the bash shell is used in the following.

1. Install CentOS 7 (minimal install or more) or start up cloud instance (e.g. Amazon AWS) with Centos 7

2. Install standard packages as root (su)

    yum update
    yum install -y net-tools
    yum install -y gcc-gfortran
    yum install -y gcc-c++
    yum install -y screen
    yum install -y wget
    yum install -y cmake
    yum install -y bzip2
    yum install -y texinfo
    yum install -y autogen
    yum install -y dejagnu
    yum install -y byacc
    yum install -y curl-devel
    yum install -y m4
    yum install -y git
    yum install -y libxml2-devel

3. Install thirdparty libraries

    mkdir -p /usr/local/src && cd /usr/local/src

    # gcc-8.3.0
    wget http://ftp.mirrorservice.org/sites/sourceware.org/pub/gcc/releases/gcc-8.3.0/gcc-8.3.0.tar.gz
    tar -xvf gcc-8.3.0.tar.gz
    cd gcc-8.3.0
    ./contrib/download_prerequisites # if this hangs, need to change the ftp URLs to https URLs in the script (behind firewall?)
    ./configure \
    --disable-multilib \
    --enable-languages=c,c++,fortran \
    --prefix=/usr/local 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr gcc-8.3.0

    export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH"

    # mpich-3.3.1
    wget http://www.mpich.org/static/downloads/3.3.1/mpich-3.3.1.tar.gz
    tar -xvzf mpich-3.3.1.tar.gz
    cd mpich-3.3.1
    ./configure \
    --prefix=/usr/local 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr mpich-3.3.1

    # netcdf-4.7.0
    wget https://www.zlib.net/zlib-1.2.11.tar.gz
    wget https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
    # go to https://www.hdfgroup.org/downloads/hdf5/source-code/ and download hdf5-1.10.5.tar.bz2
    wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-c-4.7.0.tar.gz
    wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-fortran-4.4.5.tar.gz
    #
    tar -xvzf zlib-1.2.11.tar.gz
    cd zlib-1.2.11
    ./configure \
    --prefix=/usr/local \
    --static 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    make distclean 2>&1 | tee log.distclean
    ./configure \
    --prefix=/usr/local 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr zlib-1.2.11
    #
    tar -xvzf szip-2.1.1.tar.gz
    cd szip-2.1.1
    ./configure \
    --prefix=/usr/local 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr szip-2.1.1
    #
    tar -xvjf hdf5-1.10.5.tar.bz2
    cd hdf5-1.10.5
    CC=mpicc \
    CXX=mpicxx \
    FC=mpif90 \
    ./configure \
    --prefix=/usr/local \
    --enable-parallel \
    --with-zlib=/usr/local \
    --with-szlib=/usr/local 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr hdf5-1.10.5
    #
    tar -xvzf netcdf-c-4.7.0.tar.gz
    cd netcdf-c-4.7.0
    CC=mpicc \
    CXX=mpicxx \
    FC=mpif90 \
    ./configure \
    --prefix=/usr/local \
    --enable-parallel-tests 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr netcdf-c-4.7.0
    #
    tar -xvzf netcdf-fortran-4.4.5.tar.gz
    cd netcdf-fortran-4.4.5
    CC=mpicc \
    CXX=mpicxx \
    FC=mpif90 \
    ./configure \
    --prefix=/usr/local \
    --enable-parallel-tests 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    #
    export NETCDF=/usr/local

    # NCEP libraries
    git clone https://github.com/NCAR/NCEPlibs.git
    mv NCEPlibs NCEPlibs-20190820
    cd NCEPlibs-20190820
    mkdir /usr/local/NCEPlibs-20190820
    ./make_ncep_libs.sh -s linux -c gnu -d /usr/local/NCEPlibs-20190820 -o 1
    cd ..
    rm -fr NCEPlibs-20190820

    export NCEPLIBS_DIR=/usr/local/NCEPlibs-20190820

    # Download ESMF_8_0_0_beta_snapshot_50 from https://sourceforge.net/p/esmf/esmf/ref/master/tags/
    # to /home/ubuntu/src (creates a directory esmf-esmf-...), rename it to esmf-8.0.0_bs50 and tar it
    # up for later use
    cd /usr/local/src
    mv esmf-esmf-... esmf-8.0.0_bs50
    tar -cvzf esmf-8.0.0_bs50.tar.gz esmf-8.0.0_bs50
    #
    # Install esmf-8.0.0_bs50
    tar -xvzf esmf-8.0.0_bs50.tar.gz
    cd esmf-8.0.0_bs50
    export ESMF_DIR=`pwd`
    export ESMF_INSTALL_PREFIX=/usr/local/esmf-8.0.0_bs50
    export ESMF_CXXCOMPILER=mpicxx
    export ESMF_CXXLINKER=mpicxx
    export ESMF_F90COMPILER=mpif90
    export ESMF_F90LINKER=mpif90
    export ESMF_COMM=mpich3
    export ESMF_MPIRUN=mpiexec
    export ESMF_NETCDF=nc-config
    export ESMF_INSTALL_BINDIR=bin
    export ESMF_INSTALL_LIBDIR=lib
    export ESMF_INSTALL_MODDIR=mod
    make info 2>&1 | tee log.info
    make 2>&1 | tee log.make
    # "make check" is optional and can take very long time
    make check 2>&1 | tee log.check
    # SYSTEM TESTS SUMMARY
    # Found 45 multi-processor system tests, 44 passed and 1 failed.
    # UNIT TESTS SUMMARY
    # Found 3466 non-exhaustive multi-processor unit tests, 3466 passed and 0 failed.
    # --> ignore and proceed
    make install 2>&1 | tee log.install
    make installcheck 2>&1 | tee log.installcheck
    cd ..
    rm -fr esmf-8.0.0_bs50
    export -n ESMF_DIR
    export -n ESMF_INSTALL_PREFIX
    export -n ESMF_CXXCOMPILER
    export -n ESMF_CXXLINKER
    export -n ESMF_F90COMPILER
    export -n ESMF_F90LINKER
    export -n ESMF_COMM
    export -n ESMF_MPIRUN
    export -n ESMF_NETCDF
    export -n ESMF_INSTALL_BINDIR
    export -n ESMF_INSTALL_LIBDIR
    export -n ESMF_INSTALL_MODDIR

    export ESMFMKFILE=/usr/local/esmf-8.0.0_bs50/lib/esmf.mk

4. Compile NEMSfv3gfs as normal user

    mkdir ~/scratch

    # clone this branch of NEMSfv3gfs to ~/scratch/NEMSfv3gfs

    cd ~/scratch/NEMSfv3gfs/tests

    export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH"
    export NCEPLIBS_DIR=/usr/local/NCEPlibs-20190820
    export NETCDF=/usr/local
    export ESMFMKFILE=/usr/local/esmf-8.0.0_bs50/lib/esmf.mk
    export CC=mpicc
    export CXX=mpicxx
    export F77=mpif77
    export F90=mpif90
    export FC=mpif90

    ./compile.sh $PWD/../FV3 linux.gnu 'CCPP=N' 2>&1 | tee log.compile # without CCPP
    ./compile.sh $PWD/../FV3 linux.gnu 'CCPP=Y' 2>&1 | tee log.compile # with CCPP, dynamic mode
    ./compile.sh $PWD/../FV3 linux.gnu 'CCPP=Y STATIC=Y SUITES=FV3_GFS_2017_gfdlmp' 2>&1 | tee log.compile # with CCPP, static mode, GFS suite

5. Set up the run directory using the template on Theia or Cheyenne at some location on your machine:

    a) copy the contents of the run directory templates to where you want to run the model, change to this directory
       (these folders are read-only, i.e. users might have to add the write-flag after copying/rsyncing them)

        theia:    /scratch4/BMC/gmtb/Dom.Heinzeller/rundirs/20190811/linux/fv3_gfdlmp/
        cheyenne: /glade/p/ral/jntp/GMTB/NEMSfv3gfs/rundirs/20190811/linux/fv3_gfdlmp/

    b) edit run_linux.sh, set variables and change the variable FV3_BUILD_DIR to the top-level directory of your FV3-build

    c) set the environment variables as above (consider creating a shell script to source?) and run the model

       export NCEPLIBS_DIR=/usr/local/NCEPlibs-20190820
       export NETCDF=/usr/local
       export ESMFMKFILE=/usr/local/esmf-8.0.0_bs50/lib/esmf.mk

       ./run_linux.sh 2>&1 | tee run_linux.log
       # or, with N OpenMP threads (use N=1 for the dynamic CCPP build)
       OMP_NUM_THREADS=N ./run_linux.sh 2>&1 | tee run_linux.log

    d) go and get yourself a cup of coffee ...
