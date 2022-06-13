.. ---------------------------------------------------
.. This file is included in the InputsOutputs.rst file
.. ---------------------------------------------------

.. _InputNML:

---------------------------
Namelist file *input.nml*
---------------------------

The atmosphere model reads many parameters from a Fortran namelist file, named *input.nml*.  This file contains
several Fortran namelist records, some of which are always required, others of which are only used when selected
physics options are chosen.

The following link provides an in depth description of the namelist settings:

https://dtcenter.ucar.edu/GMTB/v5.0.0/sci_doc/

The following link describes the various physics-related namelist records:

https://dtcenter.ucar.edu/GMTB/v5.0.0/sci_doc/CCPPsuite_nml_desp.html

The following link describes the stochastic physics namelist records:

https://stochastic-physics.readthedocs.io/en/ufs-v1.1.0/namelist_options.html

The following link describes some of the other namelist records (dynamics, grid, etc):

https://noaa-emc.github.io/FV3_Dycore_ufs-v1.1.0/html/index.html

The namelist section ``&interpolator_nml`` is not used in this release, and any modifications to
it will have no effect on the model results.

----------
fms_io_nml
----------

The namelist section ``&fms_io_nml`` of ``input.nml`` contains variables that control
reading and writing of restart data in netCDF format.  There is a global switch to turn on/off
the netCDF restart options in all of the modules that read or write these files. The two namelist
variables that control the netCDF restart options are ``fms_netcdf_override`` and ``fms_netcdf_restart``.
The default values of both flags are .true., so by default, the behavior of the entire model is
to use netCDF IO mode. To turn off netCDF restart, simply set ``fms_netcdf_restart`` to .false..
The namelist variables used in ``&fms_io_nml`` are described in :numref:`Table %s <fms_io_nml>`.

.. _fms_io_nml:

.. list-table:: *Description of the &fms_io_nml namelist section.*
   :widths: 25 40 15 10
   :header-rows: 1

   * - Variable Name
     - Description
     - Data Type
     - Default Value
   * - fms_netcdf_override
     - If true, ``fms_netcdf_restart`` overrides the individual ``do_netcdf_restart`` value.  If false, individual module settings has a precedence over the global setting, therefore ``fms_netcdf_restart`` is ignored.
     - logical
     - .true.
   * - fms_netcdf_restart
     - If true, all modules using restart files will operate under netCDF mode.  If false, all modules using restart files will operate under binary mode.  This flag is effective only when ``fms_netcdf_override`` is .true. When ``fms_netcdf_override`` is .false., individual module setting takes over.
     - logical
     - .true.
   * - threading_read
     - Can be 'single' or 'multi'
     - character(len=32)
     - 'multi'
   * - format
     - Format of restart data.  Only netCDF format is supported in fms_io.
     - character(len=32)
     - 'netcdf'
   * - read_all_pe
     - Reading can be done either by all PEs (default) or by only the root PE.
     - logical
     - .true.
   * - iospec_ieee32
     - If set, call mpp_open single 32-bit ieee file for reading.
     - character(len=64)
     - '-N ieee_32'
   * - max_files_w
     - Maximum number of write files
     - integer
     - 40
   * - max_files_r
     - Maximum number of read files
     - integer
     - 40
   * - time_stamp_restart
     -  If true, ``time_stamp`` will be added to the restart file name as a prefix.
     - logical
     - .true.
   * - print_chksum
     - If true, print out chksum of fields that are read and written through save_restart/restore_state.
     - logical
     - .false.
   * - show_open_namelist_file_warning
     - Flag to warn that open_namelist_file should not be called when INTERNAL_FILE_NML is defined.
     - logical
     - .false.
   * - debug_mask_list
     - Set ``debug_mask_list`` to true to print out mask_list reading from mask_table.
     - logical
     - .false.
   * - checksum_required
     -  If true, compare checksums stored in the attribute of a field against the checksum after reading in the data.
     - logical
     - .true.

This release of the UFS Weather Model sets the following variables in the ``&fms_io_nml`` namelist:

.. code-block:: console

   &fms_io_nml
     checksum_required = .false.
     max_files_r = 100
     max_files_w = 100
   /

------
namsfc
------

The namelist section ``&namsfc`` contains the filenames of the static datasets (i.e., *fix files*).
:numref:`Table %s <FixFiles>` contains a brief description of the climatological information in these files.
The variables used in ``&namsfc`` to set the filenames are described in :numref:`Table %s <namsfc_nml>`.

.. _namsfc_nml:

.. list-table:: *List of common variables in the *namsfc* namelist section used to set the filenames of
                static datasets.*
   :widths: 15 40 15 20
   :header-rows: 1

   * - Variable Name
     - File contains
     - Data Type
     - Default Value
   * - fnglac
     - Climatological glacier data
     - character*500
     - 'global_glacier.2x2.grb'
   * - fnmxic
     - Climatological maximum ice extent
     - character*500
     - 'global_maxice.2x2.grb'
   * - fntsfc
     - Climatological surface temperature
     - character*500
     - 'global_sstclim.2x2.grb'
   * - fnsnoc
     - Climatological snow depth
     - character*500
     - 'global_snoclim.1.875.grb'
   * - fnzorc
     - Climatological surface roughness
     - character*500
     - 'global_zorclim.1x1.grb'
   * - fnalbc
     - Climatological snowfree albedo
     - character*500
     - 'global_albedo4.1x1.grb'
   * - fnalbc2
     - Four albedo fields for seasonal mean climatology
     - character*500
     - 'global_albedo4.1x1.grb'
   * - fnaisc
     - Climatological sea ice
     - character*500
     - 'global_iceclim.2x2.grb'
   * - fntg3c
     - Climatological deep soil temperature
     - character*500
     - 'global_tg3clim.2.6x1.5.grb'
   * - fnvegc
     - Climatological vegetation cover
     - character*500
     - 'global_vegfrac.1x1.grb'
   * - fnvetc
     - Climatological vegetation type
     - character*500
     - 'global_vegtype.1x1.grb'
   * - fnsotc
     - Climatological soil type
     - character*500
     - 'global_soiltype.1x1.grb'
   * - fnsmcc
     - Climatological soil moisture
     - character*500
     - 'global_soilmcpc.1x1.grb'
   * - fnmskh
     - High resolution land mask field
     - character*500
     - 'global_slmask.t126.grb'
   * - fnvmnc
     - Climatological minimum vegetation cover
     - character*500
     - 'global_shdmin.0.144x0.144.grb'
   * - fnvmxc
     - Climatological maximum vegetation cover
     - character*500
     - 'global_shdmax.0.144x0.144.grb'
   * - fnslpc
     - Climatological slope type
     - character*500
     - 'global_slope.1x1.grb'
   * - fnabsc
     - Climatological maximum snow albedo
     - character*500
     - 'global_snoalb.1x1.grb'

A sample subset of this namelist is shown below:

.. code-block:: console

   &namsfc
     FNGLAC   = 'global_glacier.2x2.grb'
     FNMXIC   = 'global_maxice.2x2.grb'
     FNTSFC   = 'RTGSST.1982.2012.monthly.clim.grb'
     FNSNOC   = 'global_snoclim.1.875.grb'
     FNZORC   = 'igbp'
     FNALBC   = 'global_snowfree_albedo.bosu.t126.384.190.rg.grb'
     FNALBC2  = 'global_albedo4.1x1.grb'
     FNAISC   = 'CFSR.SEAICE.1982.2012.monthly.clim.grb'
     FNTG3C   = 'global_tg3clim.2.6x1.5.grb'
     FNVEGC   = 'global_vegfrac.0.144.decpercent.grb'
     FNVETC   = 'global_vegtype.igbp.t126.384.190.rg.grb'
     FNSOTC   = 'global_soiltype.statsgo.t126.384.190.rg.grb'
     FNSMCC   = 'global_soilmgldas.t126.384.190.grb'
     FNMSKH   = 'seaice_newland.grb'
     FNVMNC   = 'global_shdmin.0.144x0.144.grb'
     FNVMXC   = 'global_shdmax.0.144x0.144.grb'
     FNSLPC   = 'global_slope.1x1.grb'
     FNABSC   = 'global_mxsnoalb.uariz.t126.384.190.rg.grb'
   /

Additional variables for the ``&namsfc`` namelist can be found in the ``FV3/ccpp/physics/physics/sfcsub.F``
file.

---------------
atmos_model_nml
---------------

The namelist section ``&atmos_model_nml`` contains information used by the atmosphere model.
The variables used in ``&atmos_model_nml`` are described in :numref:`Table %s <atmos_model_nml>`.

.. _atmos_model_nml:

.. list-table:: *List of common variables in the *atmos_model_nml* namelist section.
   :widths: 10 40 15 15
   :header-rows: 1

   * - Variable Name
     - Description
     - Data Type
     - Default Value
   * - blocksize
     - Number of columns in each ``block`` sent to the physics. OpenMP threading is done over the number of blocks.
       For best performance this number should divide the number of grid cells per processor:
       ``((npx-1)*(npy-1)/(layout\_x)*(layout\_y))``. A description of these variables is provided
       `here <https://noaa-emc.github.io/FV3_Dycore_ufs-v1.1.0/html/group__Parameters__List.html>`_.
     - integer
     - 1
   * - chksum_debug
     - If true, compute checksums for all variables passed into the GFS physics, before and after each physics timestep.       This is very useful for reproducibility checking.
     - logical
     - .false.
   * - dycore_only
     - If true, only the dynamical core (and not the GFS physics) is executed when running the model,
       essentially running the model as a solo dynamical core.
     - logical
     - .false.
   * - debug
     - If true, turn on additional diagnostics for the atmospheric model.
     - logical
     - .false.
   * - sync
     - If true, initialize timing identifiers.
     - logical
     - .false.
   * - ccpp_suite
     - Name of the CCPP physics suite
     - character(len=256)
     - FV3_GFS_v15p2, set in ``build.sh``
   * - avg_max_length
     - Forecast interval (in seconds) determining when the maximum values of diagnostic fields in FV3
       dynamics are computed.
     - real
     - 3600.

A sample of this namelist is shown below:

.. code-block:: console

   &atmos_model_nml
     blocksize = 32
     chksum_debug = .false.
     dycore_only = .false.
     ccpp_suite = 'FV3_GFS_v16beta'
   /

The namelist section relating to the FMS diagnostic manager ``&diag_manager_nml`` is described in :numref:`Section %s <DiagManagerNML>`.

---------------
gfs_physics_nml
---------------

The namelist section ``&gfs_physics_nml`` contains physics-related information used by the atmosphere model and
some of the variables are only relevant for specific parameterizations and/or configurations.
The small set of variables used in ``&gfs_physics_nml`` are described in :numref:`Table %s <gfs_physics_nml>`.

.. _gfs_physics_nml:

.. list-table:: *List of common variables in the *gfs_physics_nml* namelist section.
   :widths: 10 40 15 15
   :header-rows: 1

   * - Variable Name
     - Description
     - Data Type
     - Default Value
   * - cplflx
     - Flag to activate atmosphere-ocean coupling. If true, turn on receiving exchange fields from other components such as ocean.
     - logical
     - .false.
   * - use_med_flux
     - Flag to receive atmosphere-ocean fluxes from mediator. If true, atmosphere-ocean fluxes will be received into the CCPP physics and used there, instead of calculating them.
     - logical
     - .false.

A sample subset of this namelist is shown below:

.. code-block:: console

   &gfs_physics_nml
     use_med_flux = .true.
     cplflx       = .true.
   /

Additional variables for the ``&gfs_physics_nml`` namelist can be found in the ``FV3/ccpp/data/GFS_typedefs.F90``
file.
