SHELL = /bin/sh

inside_nems := $(wildcard ../../../conf/configure.nems)
ifneq ($(strip $(inside_nems)),)
    include ../../../conf/configure.nems
else
    exist_configure_fv3 := $(wildcard ../FV3/conf/configure.fv3)
    ifneq ($(strip $(exist_configure_fv3)),)
        include ../FV3/conf/configure.fv3
    else
        $(error "../FV3/conf/configure.fv3 file is missing. Run ./configure")
    endif
    $(info )
    $(info Build standalone FV3 stochastic_physics ...)
    $(info )
endif

LIBRARY  = libstochastic_physics.a

FFLAGS   += -I../FV3/gfsphysics/ -I../FV3/atmos_cubed_sphere -I$(FMS_DIR)

SRCS_F   =

SRCS_f90 =  \
                ./plumes.f90

SRCS_f   =  \
		./stochy_gg_def.f                           \
		./stochy_layout_lag.f                       \
		./four_to_grid_stochy.f                     \
		./glats_stochy.f                            \
		./sumfln_stochy.f                           \
		./gozrineo_stochy.f                         \
		./get_ls_node_stochy.f                      \
		./get_lats_node_a_stochy.f                  \
		./setlats_a_stochy.f                        \
		./setlats_lag_stochy.f                      \
		./epslon_stochy.f                           \
		./getcon_lag_stochy.f                       \
		./pln2eo_stochy.f                           \
		./dozeuv_stochy.f                           \
		./dezouv_stochy.f

SRCS_F90 = \
		./kinddef.F90                               \
		./mpi_wrapper.F90                           \
		./halo_exchange.fv3.F90                     \
		./spectral_layout.F90                       \
		./getcon_spectral.F90                       \
		./stochy_namelist_def.F90                   \
		./compns_stochy.F90                         \
		./stochy_internal_state_mod.F90             \
		./stochastic_physics.F90                    \
		./stochy_patterngenerator.F90               \
		./stochy_data_mod.F90                       \
		./get_stochy_pattern.F90                    \
		./initialize_spectral_mod.F90               \
		./cellular_automata_global.F90              \
		./cellular_automata_sgs.F90                 \
		./update_ca.F90                             \
                ./lndp_apply_perts.F90

SRCS_c   =

DEPEND_FILES = $(SRCS_f) $(SRCS_f90) $(SRCS_F) $(SRCS_F90)

OBJS_f   = $(SRCS_f:.f=.o)
OBJS_f90 = $(SRCS_f90:.f90=.o)
OBJS_F   = $(SRCS_F:.F=.o)
OBJS_F90 = $(SRCS_F90:.F90=.o)
OBJS_c   = $(SRCS_c:.c=.o)

OBJS = $(OBJS_f) $(OBJS_f90) $(OBJS_F) $(OBJS_F90) $(OBJS_c)

all default: depend $(LIBRARY)

$(LIBRARY): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

.PHONY: clean
clean:
	@echo "Cleaning stochastic_physics ... "
	@echo
	$(RM) -f $(LIBRARY) *__genmod.f90 *.o */*.o *.mod *.i90 *.lst *.i depend

MKDEPENDS = ../FV3/mkDepends.pl
include ../FV3/conf/make.rules

# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
    -include depend
endif

