#
# $Id: fre-nctools.mk,v 20.0 2013/12/14 00:32:02 fms Exp $
# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build Regridding Executables
# ------------------------------------------------------------------------------
# afy    Ver   1.00  Initial version (Makefile, ver 17.0.4.2)       June 10
# afy    Ver   1.01  Add rules to build MPI-based executable        June 10
# afy    Ver   1.02  Simplified according to fre-nctools standards  June 10
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2010
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
#
# MPICC and CC are defined in env.$(SITE)
include ./env.$(SITE)

#MPICC    := mpicc
#CC       := icc
CFLAGS   := -O3 -g -traceback
CFLAGS_O2:= -O2 -g -traceback
INCLUDES := -I${NETCDF_HOME}/include -I./ -I../shared -I../../shared/mosaic
CLIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz -limf $(CLIBS2) $(STATIC)

TARGETS  := fregrid fregrid_parallel

SOURCES  := fregrid.c bilinear_interp.c conserve_interp.c fregrid_util.c
SOURCES  += create_xgrid.c gradient_c2l.c interp.c read_mosaic.c
SOURCES  += mpp_domain.c mpp_io.c tool_util.c affinity.c

OBJECTS  := $(SOURCES:c=o)

HEADERS = fre-nctools.mk ../shared/mpp.h  ../shared/mpp_domain.h  ../shared/mpp_io.h ../shared/tool_util.h   \
          ../../shared/mosaic/constant.h ../../shared/mosaic/create_xgrid.h  \
          ../../shared/mosaic/gradient_c2l.h ../../shared/mosaic/interp.h  \
          ../../shared/mosaic/mosaic_util.h  ../../shared/mosaic/read_mosaic.h 

all: $(TARGETS)

fregrid: $(OBJECTS) mosaic_util.o mpp.o
	$(CC) -o $@ $^ $(CLIBS)

fregrid_parallel: $(OBJECTS) mosaic_util_parallel.o mpp_parallel.o
	$(MPICC) -o $@ $^ $(CLIBS)

mosaic_util.o: ../../shared/mosaic/mosaic_util.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

mosaic_util_parallel.o: ../../shared/mosaic/mosaic_util.c $(HEADERS)
	$(MPICC) -Duse_libMPI $(CFLAGS) $(INCLUDES) -o $@ -c $< 

read_mosaic.o: ../../shared/mosaic/read_mosaic.c $(HEADERS)
	$(CC) -Duse_netCDF $(CFLAGS) $(INCLUDES) -c $< 

gradient_c2l.o: ../../shared/mosaic/gradient_c2l.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

interp.o: ../../shared/mosaic/interp.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

mpp_io.o: ../shared/mpp_io.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

mpp_domain.o: ../shared/mpp_domain.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

affinity.o: ../shared/affinity.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

mpp.o: ../shared/mpp.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

mpp_parallel.o: ../shared/mpp.c $(HEADERS)
	$(MPICC) -Duse_libMPI $(CFLAGS) $(INCLUDES) -o $@ -c $< 

tool_util.o: ../shared/tool_util.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

conserve_interp.o: conserve_interp.c $(HEADERS)
	$(CC) $(CFLAGS_O2) $(INCLUDES) -c $< 

bilinear_interp.o: bilinear_interp.c $(HEADERS)
	$(CC) $(CFLAGS_O2) $(INCLUDES) -c $< 

create_xgrid.o: ../../shared/mosaic/create_xgrid.c $(HEADERS)
	$(CC) $(CFLAGS_O2) $(INCLUDES) -c $< 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	-rm -f *.o $(TARGETS)
