SHELL = /bin/sh

inside_nems := $(wildcard ../../../conf/configure.nems)
ifneq ($(strip $(inside_nems)),)
    include ../../../conf/configure.nems
else
    exist_configure_fv3 := $(wildcard ../conf/configure.fv3)
    ifneq ($(strip $(exist_configure_fv3)),)
        include ../conf/configure.fv3
    else
        $(error "../conf/configure.fv3 file is missing. Run ./configure")
    endif
    $(info )
    $(info Build standalone FV3 gfsphysics ...)
    $(info )
endif

ifneq (,$(findstring MACOSX,$(CPPDEFS)))
LIBRARY = libgfsphys.dylib
else
LIBRARY = libgfsphys.so
endif
VER_MAJOR = 1
VER_MINOR = 0
VER_PATCH = 0

FFLAGS   += -I../fms -I../fms/include -fPIC

CPPDEFS += -DNEW_TAUCTMAX -DSMALL_PE -DNEMS_GSM

ifneq (,$(findstring CCPP,$(CPPDEFS)))
include ./CCPP_SCHEMES.mk
SRCS_f   = $(SCHEMES_f)
SRCS_F   = $(SCHEMES_F)
SRCS_f90 = $(SCHEMES_f90)
SRCS_F90 = $(SCHEMES_F90)
else
SRCS_f   =  \
	   ./physics/cnvc90.f                                                        \
	   ./physics/date_def.f                                                      \
	   ./physics/dcyc2.f                                                         \
	   ./physics/gfs_phy_tracer_config.f                                         \
	   ./physics/gocart_tracer_config_stub.f                                     \
	   ./physics/gscond.f                                                        \
	   ./physics/gwdc.f                                                          \
	   ./physics/gwdps.f                                                         \
	   ./physics/h2o_def.f                                                       \
	   ./physics/iounitdef.f                                                     \
	   ./physics/mersenne_twister.f                                              \
	   ./physics/mfdeepcnv.f                                                     \
	   ./physics/mfpbl.f                                                         \
	   ./physics/mfshalcnv.f                                                     \
	   ./physics/module_bfmicrophysics.f                                         \
	   ./physics/moninedmf.f                                                     \
	   ./physics/namelist_soilveg.f                                              \
	   ./physics/ozne_def.f                                                      \
	   ./physics/ozphys.f                                                        \
	   ./physics/physparam.f                                                     \
	   ./physics/precpd.f                                                        \
	   ./physics/rad_initialize.f                                                \
	   ./physics/radiation_aerosols.f                                            \
	   ./physics/radiation_astronomy.f                                           \
	   ./physics/radiation_clouds.f                                              \
	   ./physics/radiation_gases.f                                               \
	   ./physics/radiation_surface.f                                             \
	   ./physics/radlw_datatb.f                                                  \
	   ./physics/radlw_main.f                                                    \
	   ./physics/radlw_param.f                                                   \
	   ./physics/radsw_datatb.f                                                  \
	   ./physics/radsw_main.f                                                    \
	   ./physics/radsw_param.f                                                   \
	   ./physics/rascnvv2.f                                                      \
	   ./physics/rayleigh_damp.f                                                 \
	   ./physics/set_soilveg.f                                                   \
	   ./physics/GFS_surface_loop_control.f                                      \
	   ./physics/sfc_diag.f                                                      \
	   ./physics/sfc_diff.f                                                      \
	   ./physics/sfc_drv.f                                                       \
	   ./physics/sfc_nst.f                                                       \
	   ./physics/sfc_sice.f                                                      \
	   ./physics/sflx.f                                                          \
	   ./physics/tridi.f
SRCS_f90 = \
	   ./physics/GFS_calpreciptype.f90                                           \
	   ./physics/GFS_MP_generic_post.f90                                         \
	   ./physics/GFS_MP_generic_pre.f90                                          \
	   ./physics/GFS_zhao_carr_pre.f90                                           \
	   ./physics/GFS_rad_time_vary.fv3.f90                                       \
	   ./physics/GFS_radupdate.f90                                               \
	   ./physics/funcphys.f90                                                    \
	   ./physics/gcycle.f90                                                      \
	   ./physics/get_prs_fv3.f90                                                 \
	   ./physics/GFS_DCNV_generic.f90                                            \
	   ./physics/GFS_SCNV_generic.f90                                            \
	   ./physics/GFS_PBL_generic.f90                                             \
	   ./physics/GFS_suite_interstitial.ipd.F90                                  \
	   ./physics/GFS_phys_time_vary.fv3.f90                                      \
	   ./physics/GFS_stochastics.f90                                             \
	   ./physics/GFS_surface_generic.f90                                         \
	   ./physics/h2ointerp.f90                                                   \
	   ./physics/module_nst_model.f90                                            \
	   ./physics/module_nst_parameters.f90                                       \
	   ./physics/module_nst_water_prop.f90                                       \
	   ./physics/ozinterp.f90                                                    \
	   ./physics/physcons.f90                                                    \
	   ./physics/radcons.f90                                                     \
	   ./physics/wam_f107_kp_mod.f90
SRCS_F   = \
	   ./physics/aer_cloud.F                                                     \
	   ./physics/cldwat2m_micro.F                                                \
	   ./physics/machine.F                                                       \
	   ./physics/sfcsub.F                                                        \
	   ./physics/num_parthds.F                                                   \
	   ./physics/wv_saturation.F
SRCS_F90 = \
	   ./physics/GFDL_parse_tracers.F90                                          \
	   ./GFS_layer/GFS_abstraction_layer.F90                                     \
	   ./GFS_layer/GFS_diagnostics.F90                                           \
	   ./GFS_layer/GFS_driver.F90                                                \
	   $(GFS_SUITE_INTERSTITIAL)                                                 \
	   ./physics/GFS_rrtmg_pre.F90                                               \
	   ./physics/GFS_rrtmg_post.F90                                              \
	   ./physics/memcheck.F90                                                    \
	   ./physics/rrtmg_sw_pre.F90                                                \
	   ./physics/rrtmg_sw_post.F90                                               \
	   ./physics/rrtmg_lw_pre.F90                                                \
	   ./physics/rrtmg_lw_post.F90                                               \
	   ./physics/GFS_debug.F90                                                   \
	   ./GFS_layer/GFS_physics_driver.F90                                        \
	   ./GFS_layer/GFS_radiation_driver.F90                                      \
	   ./GFS_layer/GFS_restart.F90                                               \
	   ./GFS_layer/GFS_typedefs.F90                                              \
	   ./IPD_layer/IPD_driver.F90                                                \
	   ./IPD_layer/IPD_typedefs.F90
endif

SRCS_c   =

ifneq (,$(findstring CCPP,$(CPPDEFS)))
include ./CCPP_CAPS.mk
else
CAPS_F90 = 
endif

DEPEND_FILES = $(SRCS_f) $(SRCS_f90) $(SRCS_F) $(SRCS_F90) $(CAPS_F90)

OBJS_f   = $(SRCS_f:.f=.o)
OBJS_f90 = $(SRCS_f90:.f90=.o)
OBJS_F   = $(SRCS_F:.F=.o)
OBJS_F90 = $(SRCS_F90:.F90=.o)
OBJS_c   = $(SRCS_c:.c=.o)

OBJS = $(OBJS_f) $(OBJS_f90) $(OBJS_F) $(OBJS_F90) $(OBJS_c)

CAPS = $(CAPS_F90:.F90=.o)

all default: depend $(LIBRARY)

ifneq (,$(findstring MACOSX,$(CPPDEFS)))
LIBRARY_FULL_NAME = $(subst .dylib,.$(VER_MAJOR).$(VER_MINOR).$(VER_PATCH).dylib,$(LIBRARY))
$(LIBRARY): $(OBJS) $(CAPS)
	$(FC) -shared $(OBJS) $(CAPS) $(LDFLAGS) $(NCEPLIBS) -o $(LIBRARY_FULL_NAME)
	ln -sf $(LIBRARY_FULL_NAME) $(LIBRARY)
	ln -sf $(LIBRARY_FULL_NAME) $(subst .dylib,.$(VER_MAJOR).dylib,$(LIBRARY))
else
$(LIBRARY): $(OBJS) $(CAPS)
	$(FC) -shared -Wl,-soname,$(LIBRARY).$(VER_MAJOR) $(OBJS) $(CAPS) $(LDFLAGS) $(NCEPLIBS) -o $(LIBRARY).$(VER_MAJOR).$(VER_MINOR).$(VER_PATCH)
	ln -sf $(LIBRARY).$(VER_MAJOR).$(VER_MINOR).$(VER_PATCH) $(LIBRARY)
	ln -sf $(LIBRARY).$(VER_MAJOR).$(VER_MINOR).$(VER_PATCH) $(LIBRARY).$(VER_MAJOR)
endif

# this is the place to override default (implicit) compilation rules
# and create specific (explicit) rules

# this has no effect, because radiation_aerosols.o is in physics, not in gfsphys
./radiation_aerosols.o : ./gfsphys/radiation_aerosols.f
	$(FC) $(CPPDEFS) $(FFLAGS) $(OTHER_FFLAGS) -xCORE-AVX-I -c $< -o $@

# Reduce optimization for sfc_sice for bit-for-bit reproducibility
FFLAGS_REDUCED_OPT=$(subst -O2,-O1,$(subst -xCORE-AVX2,-xCORE-AVX-I,$(FFLAGS)))
./physics/sfc_sice.o : ./physics/sfc_sice.f
	$(FC) $(CPPDEFS) $(FFLAGS_REDUCED_OPT) $(OTHER_FFLAGS) -c $< -o $@

./GFS_layer/GFS_diagnostics.o : ./GFS_layer/GFS_diagnostics.F90
	$(FC) $(CPPDEFS) $(FFLAGS) $(OTHER_FFLAGS) -O0 -c $< -o $@

ifneq (,$(findstring PGIFIX,$(CPPDEFS)))
$(CAPS):
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHER_FFLAGS) -c $< -o $@
	# Apply a fix specific to the PGI compiler (rename objects in cap object files)
	./pgifix.py $@
else
$(CAPS):
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHER_FFLAGS) -c $< -o $@
endif

# Do preprocessing of the IPD-CCPP driver in two steps to be
# able to look at the actual .f90 file that gets compiled
./IPD_layer/IPD_CCPP_driver.o: ./IPD_layer/IPD_CCPP_driver.F90
	$(CPP) $(CPPDEFS) $(CPPFLAGS) $< > $*.tmp.f90
	$(FC) $(FFLAGS) $(OTHER_FFLAGS) -c $*.tmp.f90 -o $@

.PHONY: clean
clean:
	@echo "Cleaning gfsphysics  ... "
	@echo
	$(RM) -f $(LIBRARY) *__genmod.f90 *.o */*.o *.tmp.f90 */*.tmp.f90 *.mod *.lst *.i depend
	$(RM) -f $(LIBRARY).$(VER_MAJOR)
	$(RM) -f $(LIBRARY).$(VER_MAJOR).$(VER_MINOR).$(VER_PATCH)

MKDEPENDS = ../mkDepends.pl
include ../conf/make.rules

include ./depend

# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
   -include depend
endif
