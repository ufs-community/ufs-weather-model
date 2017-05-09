#!/bin/tcsh -f
#
# $Id: fre-nctools-chaco-test,v 1.2 2015/01/16 21:47:07 fms Exp $
# ------------------------------------------------------------------------------
# FMS/FRE Project: Program to Create the fre-nctools Package
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2012
# Designed and written by V Balaji, Amy Langenhorst and Aleksey Yakovlev
#
set pkgName = fre-nctools
set home_dir = `pwd`/../..

#set pkgVersion = $1
#echo $pkgVersion

#Loop though these directories to build the tools there.
set freNCToolsSrc = 'tools/{make_hgrid,make_solo_mosaic,fregrid}'

# This first run assumes you are in the git directory
set srcDir = `pwd`

#Build in a temporary directory.
set tmpDir = `pwd`/build    #$HOME/tmp.$$
mkdir -p $tmpDir
if ( $status ) then
  echo "Error during mkdir '$tmpDir'..."
  exit 1
endif
pushd $tmpDir

echo "////////////////////////////////////////////////////////////////////////////////"
echo "//////////////////////////////////////////////////////// Environment Settings //"
echo "////////////////////////////////////////////////////////////////////////////////"

#Original setup is for cray so for now require input only on a different platform.
set system_site = ${1}
if ( "$system_site" == "" ) then
  echo "Usage: BUILD_TOOLS.csh wcoss_cray or theia"
  exit 1
else

source $MODULESHOME/init/tcsh
#source ${PWD}/../../../IC_scripts/ENV.GAEA
#setenv FRE_SYSTEM_SITE gaea
source ${PWD}/../../../modulefiles/fv3gfs/fre-nctools.${system_site}
setenv FRE_SYSTEM_SITE ${system_site}

setenv MPICH_UNEX_BUFFER_SIZE 256m
setenv MPICH_MAX_SHORT_MSG_SIZE 64000
setenv MPICH_PTL_UNEX_EVENTS 160k
setenv KMP_STACKSIZE 2g
setenv F_UFMTENDIAN big
if ( ${system_site} == theia ) then
  setenv HDF5_DIR $HDF5
  setenv NETCDF_DIR $NETCDF
endif
alias make make HDF5_HOME=${HDF5_DIR}  NETCDF_HOME=${NETCDF_DIR} NC_BLKSZ=64K SITE=${FRE_SYSTEM_SITE} -f fre-nctools.mk

module list

echo "////////////////////////////////////////////////////////////////////////////////"
echo "//////////////////////////////////////////////////////////// Directory Layout //"
echo "////////////////////////////////////////////////////////////////////////////////"

mkdir -p share/src
cp -r ${srcDir}/{shared,tools} share/src
mkdir -p $FRE_SYSTEM_SITE/bin
mkdir -p share/bin

echo "Done..."

foreach freNCToolsDir ( $freNCToolsSrc )
  echo "////////////////////////////////////////////////////////////////////////////////"
  echo "////////////////////////////////////////////////////////////////// $freNCToolsDir:t"
  echo "////////////////////////////////////////////////////////////////////////////////"

  pushd share/src/$freNCToolsDir
  cp fre-nctools.mk_${system_site} fre-nctools.mk
  set targets=` grep "TARGETS  :=" fre-nctools.mk | cut -f2 -d'=' `
  if ( $?NOPARALLEL ) then
     set targets=` grep "TARGETS  :=" fre-nctools.mk | cut -f2 -d'=' | sed 's/ \S*_parallel/ /'`
  endif
  echo "Making $targets"

  make clean
  make
  if ( $status ) then
   echo "Error: make failed for $targets"
   exit 1
  endif

  foreach target ( $targets )
    if ( -f $target ) then
      #mv $target $tmpDir/$FRE_SYSTEM_SITE/bin
      mv $target $home_dir/exec
    else
      echo "Error during '$target' build"
      exit 1
    endif
  end
  make clean
  popd
end

echo "////////////////////////////////////////////////////////////////////////////////"
echo "///////////////////////////////////////////////////////////////// filter_topo //"
echo "////////////////////////////////////////////////////////////////////////////////"

cd ../tools/filter_topo
./make.csh_${system_site}
mv filter_topo $home_dir/exec/.

echo "\n////////// CLEANING UP TEMPORARY BUILD AREA //////////\n"
rm -fr $tmpDir

