#
# $Id: fre-nctools.mk,v 20.0 2013/12/14 00:33:28 fms Exp $
# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build Regridding Executables
# ------------------------------------------------------------------------------
# afy    Ver   1.00  Initial version (Makefile, ver 17.0.4.2)       June 10
# afy    Ver   1.01  Add rules to build MPI-based executable        June 10
# afy    Ver   1.02  Simplified according to fre-nctools standards  June 10
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2011
# This program is distributed under the terms of the GNU General Public
# License. See the file COPYING contained in this directory
#
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
#
include env.$(SITE)

#MPICC    := mpicc
#CC       := icc
CFLAGS   := -O3 -traceback
CFLAGS_O2:= -O2 -traceback
INCLUDES := -I${NETCDF_HOME}/include -I./ -I../shared -I../../shared/mosaic
CLIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz -limf $(CLIBS2) $(STATIC)

TARGETS  := make_hgrid make_hgrid_parallel

SOURCES  := make_hgrid.c create_conformal_cubic_grid.c create_gnomonic_cubic_grid.c create_grid_from_file.c create_lonlat_grid.c
SOURCES  += mpp_domain.c mpp_io.c tool_util.c
SOURCES  += read_mosaic.c create_xgrid.c interp.c

OBJECTS  := $(SOURCES:c=o)

HEADERS = fre-nctools.mk ../shared/mpp.h  ../shared/mpp_domain.h  ../shared/mpp_io.h ../shared/tool_util.h   \
          ../../shared/mosaic/constant.h ../../shared/mosaic/create_xgrid.h  \
          ../../shared/mosaic/interp.h  ../../shared/mosaic/mosaic_util.h  \
          ../../shared/mosaic/read_mosaic.h ./create_hgrid.h

all: $(TARGETS)

make_hgrid: $(OBJECTS) mosaic_util.o mpp.o
	$(CC) -o $@ $^ $(CLIBS)

make_hgrid_parallel: $(OBJECTS) mosaic_util_parallel.o mpp_parallel.o
	$(MPICC) -o $@ $^ $(CLIBS)

make_hgrid.o: make_hgrid.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

create_gnomonic_cubic_grid.o: create_gnomonic_cubic_grid.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

create_conformal_cubic_grid.o: create_conformal_cubic_grid.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

create_lonlat_grid.o: create_lonlat_grid.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

create_grid_from_file.o: create_grid_from_file.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

mosaic_util.o: ../../shared/mosaic/mosaic_util.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

mosaic_util_parallel.o: ../../shared/mosaic/mosaic_util.c $(HEADERS)
	$(MPICC) -Duse_libMPI $(CFLAGS) $(INCLUDES) -o $@ -c $< 

read_mosaic.o: ../../shared/mosaic/read_mosaic.c $(HEADERS)
	$(CC) -Duse_netCDF $(CFLAGS) $(INCLUDES) -c $< 

read_mosaic_parallel.o: ../../shared/mosaic/read_mosaic.c $(HEADERS)
	$(MPICC) -Duse_libMPI -Duse_netCDF $(CFLAGS) $(INCLUDES) -o $@ -c $< 

interp.o: ../../shared/mosaic/interp.c $(HEADERS)
	$(CC) -Duse_netCDF $(CFLAGS) $(INCLUDES) -c $< 

mpp_io.o: ../shared/mpp_io.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

mpp_domain.o: ../shared/mpp_domain.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

mpp.o: ../shared/mpp.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

mpp_parallel.o: ../shared/mpp.c $(HEADERS)
	$(MPICC) -Duse_libMPI $(CFLAGS) $(INCLUDES) -o $@ -c $< 

tool_util.o: ../shared/tool_util.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

create_xgrid.o: ../../shared/mosaic/create_xgrid.c $(HEADERS)
	$(CC) $(CFLAGS_O2) $(INCLUDES) -c $< 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	-rm -f *.o $(TARGETS)
