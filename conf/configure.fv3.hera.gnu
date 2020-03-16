## NEMS configuration file
##
## Platform: Hera
## Compiler: GNU with OpenMPI

SHELL=/bin/sh

################################################################################
## Include the common configuration parts

ifdef InNemsMakefile
include $(TOP)/conf/configure.nems.NUOPC
endif

############
# commands #
############
FC  = mpif90
CC  = mpicc
CXX = mpicxx
LD  = $(FC)

#########
# flags #
#########
# default is 64-bit OpenMP non-hydrostatic build using AVX2
DEBUG =
REPRO =
VERBOSE =
OPENMP = Y
AVX2 = Y
HYDRO = N
CCPP = N
STATIC = N

include       $(ESMFMKFILE)
ESMF_INC    = $(ESMF_F90COMPILEPATHS)

NEMSIOINC = -I$(NEMSIO_INC)
NCEPLIBS = $(NEMSIO_LIB) $(BACIO_LIB4) $(SP_LIBd) $(W3EMC_LIBd) $(W3NCO_LIBd)

##############################################
# Need to use at least GNU Make version 3.81 #
##############################################
need := 3.81
ok := $(filter $(need),$(firstword $(sort $(MAKE_VERSION) $(need))))
ifneq ($(need),$(ok))
$(error Need at least make version $(need).  Load module gmake/3.81)
endif

NETCDF_ROOT = $(NETCDF)
INCLUDE = -I$(NETCDF_ROOT)/include
NETCDF_INC = -I$(NETCDF_ROOT)/include
ifneq ($(findstring netcdf/4,$(LOADEDMODULES)),)
  NETCDF_LIB += -L$(NETCDF)/lib -lnetcdff -lnetcdf
else
  NETCDF_LIB += -L$(NETCDF)/lib -lnetcdff -lnetcdf
endif

FPPFLAGS := -cpp -Wp,-w $(INCLUDE)
CFLAGS := $(INCLUDE)

FFLAGS := $(INCLUDE) -fcray-pointer -ffree-line-length-none -fno-range-check -fbacktrace

CPPDEFS += -Duse_libMPI -Duse_netCDF -DSPMD -DUSE_LOG_DIAG_FIELD_INFO  -DUSE_GFSL63 -DGFS_PHYS -Duse_WRTCOMP
CPPDEFS += -DNEW_TAUCTMAX -DINTERNAL_FILE_NML -DNO_INLINE_POST

ifeq ($(HYDRO),Y)
CPPDEFS +=
else
CPPDEFS += -DMOIST_CAPPA -DUSE_COND
endif

ifeq ($(NAM_phys),Y)
CPPDEFS += -DNAM_phys
endif

ifeq ($(32BIT),Y)
CPPDEFS += -DOVERLOAD_R4 -DOVERLOAD_R8
FFLAGS +=
else
FFLAGS += -fdefault-real-8 -fdefault-double-8
endif

ifeq ($(AVX2),Y)
FFLAGS +=
CFLAGS +=
else
FFLAGS +=
CFLAGS +=
endif

ifeq ($(MULTI_GASES),Y)
CPPDEFS += -DMULTI_GASES
endif

FFLAGS_OPT = -O2 -fno-range-check
FFLAGS_REPRO = -O2 -g -fbacktrace -fno-range-check
FFLAGS_DEBUG = -g -O0 -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check -fbacktrace -fno-range-check -Wall

TRANSCENDENTALS :=
FFLAGS_OPENMP = -fopenmp
FFLAGS_VERBOSE = -v -V

CFLAGS += -D__IFC

CFLAGS_OPT = -O2
CFLAGS_REPRO = -O2
CFLAGS_OPENMP = -fopenmp
CFLAGS_DEBUG = -O0 -g

# Optional Testing compile flags.  Mutually exclusive from DEBUG, REPRO, and OPT
# *_TEST will match the production if no new option(s) is(are) to be tested.
FFLAGS_TEST = -O3
CFLAGS_TEST = -O2

LDFLAGS :=
LDFLAGS_OPENMP := -fopenmp
LDFLAGS_VERBOSE := -Wl,-V,--verbose,-cref,-M

# start with blank LIBS
LIBS :=

ifeq ($(REPRO),Y)
CPPDEFS += -DREPRO
CFLAGS += $(CFLAGS_REPRO)
FFLAGS += $(FFLAGS_REPRO)
FAST :=
else ifeq ($(DEBUG),Y)
CPPDEFS += -DDEBUG
CFLAGS += $(CFLAGS_DEBUG)
FFLAGS += $(FFLAGS_DEBUG)
FAST :=
else ifeq ($(TEST),Y)
CFLAGS += $(CFLAGS_TEST)
FFLAGS += $(FFLAGS_TEST)
FAST :=
else
CFLAGS += $(CFLAGS_OPT)
FFLAGS += $(FFLAGS_OPT)
FAST := $(TRANSCENDENTALS)
endif

ifeq ($(OPENMP),Y)
CPPDEFS += -DOPENMP
CFLAGS += $(CFLAGS_OPENMP)
FFLAGS += $(FFLAGS_OPENMP)
LDFLAGS += $(LDFLAGS_OPENMP)
endif

ifeq ($(VERBOSE),Y)
CFLAGS += $(CFLAGS_VERBOSE)
FFLAGS += $(FFLAGS_VERBOSE)
LDFLAGS += $(LDFLAGS_VERBOSE)
endif

ifeq ($(CCPP),Y)
CPPDEFS += -DCCPP
CFLAGS += -I$(PATH_CCPP)/include
FFLAGS += -I$(PATH_CCPP)/include
ifeq ($(STATIC),Y)
CPPDEFS += -DSTATIC
LDFLAGS += -L$(PATH_CCPP)/lib -lccppphys -lccpp $(NCEPLIBS) -lxml2
else
LDFLAGS += -L$(PATH_CCPP)/lib -lccpp
endif
endif

LDFLAGS += $(LIBS)

ifdef InNemsMakefile
FFLAGS += $(ESMF_INC)
CPPFLAGS += -cpp -traditional
EXTLIBS = $(NCEPLIBS) $(ESMF_LIB) $(LDFLAGS) $(NETCDF_LIB)
endif
