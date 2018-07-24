## NEMS configuration file
##
## Platform: Odin
## Compiler: Intel with MVAPICH2

SHELL=/bin/sh

################################################################################
## Include the common configuration parts

ifdef InNemsMakefile
include $(TOP)/conf/configure.nems.NUOPC
endif

###################### PHYS_MODE ##### CHEM_MODE ###############################
#
#
#

PHYS_MODE       =compile
CHEM_MODE       =compile
ifeq ($(PHYS_MODE),compile)
       PHYS_LIB = $(TOP)/atmos/gsm/gsmphys
       PHYS_INC = $(TOP)/atmos/gsm/gsmphys
       PHYS_DIR = $(TOP)/atmos/gsm/gsmphys
endif
ifeq ($(CHEM_MODE),compile)
        CHEM_LIB = $(TOP)/chem
        CHEM_INC = $(TOP)/chem/gocart/src/Config/
        CHEM_DIR = $(TOP)/chem
       CHEM_MOD = $(TOP)/chem/gocart/${ARCH}/include
      ESMADIR = chem/gocart
endif

############
# commands #
############
FC = ftn
CC = cc
CXX = CC
LD = ftn -mkl=sequential

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

NETCDF_ROOT = $(NETCDF4)
INCLUDE = -I$(NETCDF_ROOT)/include
NETCDF_INC = -I$(NETCDF_ROOT)/include
ifneq ($(findstring netcdf/4,$(LOADEDMODULES)),)
  NETCDF_LIB += -L$(NETCDF4)/lib -lnetcdff -lnetcdf
else
  NETCDF_LIB += -L$(NETCDF4)/lib -lnetcdff -lnetcdf
endif

FPPFLAGS := -fpp -Wp,-w $(INCLUDE)
CFLAGS := $(INCLUDE)

FFLAGS := $(INCLUDE) -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -nowarn -sox -align array64byte

ifeq ($(HYDRO),Y)
CPPDEFS += -Duse_libMPI -Duse_netCDF -DSPMD -DUSE_LOG_DIAG_FIELD_INFO -Duse_LARGEFILE -DUSE_GFSL63 -DGFS_PHYS -Duse_WRTCOMP
else
CPPDEFS += -Duse_libMPI -Duse_netCDF -DSPMD -DUSE_LOG_DIAG_FIELD_INFO -Duse_LARGEFILE -DUSE_GFSL63 -DGFS_PHYS -DMOIST_CAPPA -DUSE_COND -Duse_WRTCOMP
endif

CPPDEFS += -DNEW_TAUCTMAX -DINTERNAL_FILE_NML

ifeq ($(32BIT),Y)
CPPDEFS += -DOVERLOAD_R4 -DOVERLOAD_R8
FFLAGS += -i4 -real-size 32
else
FFLAGS += -i4 -real-size 64 -no-prec-div -no-prec-sqrt
endif

FFLAGS += -xCORE-AVX-I #-axavx 
CFLAGS += -xCORE-AVX-I #-axavx 

FFLAGS_OPT = -O2 -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3
FFLAGS_REPRO = -O2 -debug minimal -fp-model source -qoverride-limits -g -traceback
#FFLAGS_DEBUG = -g -O0 -check bounds -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -traceback -ftrapuv
FFLAGS_DEBUG = -g -O0 -check bounds -traceback 

TRANSCENDENTALS := -fast-transcendentals
FFLAGS_OPENMP = -qopenmp
FFLAGS_VERBOSE = -v -V -what

CFLAGS += -D__IFC -sox -fp-model source

CFLAGS_OPT = -O2 -debug minimal
CFLAGS_REPRO = -O2 -debug minimal
CFLAGS_OPENMP = -qopenmp
CFLAGS_DEBUG = -O0 -g -ftrapuv -traceback

# Optional Testing compile flags.  Mutually exclusive from DEBUG, REPRO, and OPT
# *_TEST will match the production if no new option(s) is(are) to be tested.
FFLAGS_TEST = -O3 -debug minimal -fp-model source -qoverride-limits
CFLAGS_TEST = -O2

LDFLAGS :=
LDFLAGS_OPENMP := -qopenmp
LDFLAGS_VERBOSE := -Wl,-V,--verbose,-cref,-M

# start with blank LIBS
LIBS :=

#ifneq ($(REPRO),)
#CFLAGS += $(CFLAGS_REPRO)
#FFLAGS += $(FFLAGS_REPRO)
#FAST :=
#else ifneq ($(DEBUG),)
CFLAGS += $(CFLAGS_DEBUG)
FFLAGS += $(FFLAGS_DEBUG)
FAST :=
#else ifneq ($(TEST),)
#CFLAGS += $(CFLAGS_TEST)
#FFLAGS += $(FFLAGS_TEST)
#FAST :=
#else
#CFLAGS += $(CFLAGS_OPT)
#FFLAGS += $(FFLAGS_OPT)
#FAST := $(TRANSCENDENTALS)
#endif

ifneq ($(OPENMP),)
CFLAGS += $(CFLAGS_OPENMP)
FFLAGS += $(FFLAGS_OPENMP)
LDFLAGS += $(LDFLAGS_OPENMP)
endif

ifneq ($(VERBOSE),)
CFLAGS += $(CFLAGS_VERBOSE)
FFLAGS += $(FFLAGS_VERBOSE)
LDFLAGS += $(LDFLAGS_VERBOSE)
endif

LDFLAGS += $(LIBS)

ifdef InNemsMakefile
FFLAGS += $(ESMF_INC)
CPPFLAGS += -traditional
EXTLIBS = $(NCEPLIBS) $(ESMF_LIB) $(LDFLAGS) $(NETCDF_LIB)
endif
