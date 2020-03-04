.. _InputsOutputs:
  
*************************
Inputs and Outputs
*************************

This chapter describes the input and output files needed for executing the model in the various supported configurations. 

=============
Input files
=============

There are three types of files needed to execute a run: static datasets (*fix* files containing climatological
information), files that depend on grid resolution and initial conditions, and model configuration files (such as namelists).

------------------------------------
Static datasets (i.e., *fix files*)
------------------------------------
The static input files are listed and described in :numref:`Table %s <FixFiles>`.

.. _FixFiles:

.. list-table:: *Fix files containing climatological information*
   :widths: 40 50
   :header-rows: 1

   * - Filename
     - Description
   * - aerosol.dat
     - External aerosols data file
   * - CFSR.SEAICE.1982.2012.monthly.clim.grb
     - CFS reanalysis of monthly sea ice climatology
   * - co2historicaldata_YYYY.txt
     - Monthly CO2 in PPMV data for year YYYY
   * - global_albedo4.1x1.grb
     - Four albedo fields for seasonal mean climatology: 2 for strong zenith angle dependent (visible and near IR)
       and 2 for weak zenith angle dependent
   * - global_glacier.2x2.grb
     - Glacier points, permanent/extreme features
   * - global_h2oprdlos.f77
     - Coefficients for the parameterization of photochemical production and loss of water (H2O)
   * - global_maxice.2x2.grb
     - Maximum ice extent, permanent/extreme features
   * - global_mxsnoalb.uariz.t126.384.190.rg.grb
     - Climatological maximum snow albedo
   * - global_o3prdlos.f77
     - Monthly mean ozone coefficients 
   * - global_shdmax.0.144x0.144.grb
     - Climatological maximum vegetation cover
   * - global_shdmin.0.144x0.144.grb
     - Climatological minimum vegetation cover
   * - global_slope.1x1.grb
     - Climatological slope type
   * - global_snoclim.1.875.grb
     - Climatological snow depth
   * - global_snowfree_albedo.bosu.t126.384.190.rg.grb
     - Climatological snowfree albedo
   * - global_soilmgldas.t126.384.190.grb
     - Climatological soil moisture
   * - global_soiltype.statsgo.t126.384.190.rg.grb
     - Soil type from the STATSGO dataset
   * - global_tg3clim.2.6x1.5.grb
     - Climatological deep soil temperature
   * - global_vegfrac.0.144.decpercent.grb
     - Climatological vegetation fraction
   * - global_vegtype.igbp.t126.384.190.rg.grb
     - Climatological vegetation type
   * - global_zorclim.1x1.grb
     - Climatological surface roughness
   * - RTGSST.1982.2012.monthly.clim.grb
     - Monthly, climatological, real-time global sea surface temperature
   * - seaice_newland.grb
     - High resolution land mask
   * - sfc_emissivity_idx.txt
     - External surface emissivity data table
   * - solarconstant_noaa_an.txt
     - External solar constant data table

---------------------------------------------
Grid description and initial condition files
---------------------------------------------
The input files containing grid information and the initial conditions are listed and described in :numref:`Table %s <GridICFiles>`.

.. _GridICFiles:

.. list-table:: *Input files containing grid information and initial conditions*
   :widths: 35 50 15
   :header-rows: 1

   * - Filename
     - Description
     - Date-dependent
   * - C96_grid.tile[1-6].nc
     - C96 grid information for tiles 1-6
     - 
   * - gfs_ctrl.nc
     - NCEP NGGPS tracers, ak, and bk
     - ✔
   * - gfs_data.tile[1-6].nc
     - Initial condition fields (ps, u, v, u, z, t, q, O3). May include spfo3, spfo, spf02 if multiple gases are used
     - ✔
   * - oro_data.tile[1-6].nc
     - Model terrain (topographic/orographic information) for grid tiles 1-6
     - 
   * - sfc_ctrl.nc
     - Control parameters for surface input: forecast hour, date, number of soil levels
     -
   * - sfc_data.tile[1-6].nc
     - Surface properties for grid tiles 1-6
     - ✔

------------------------------------
Model configuration files
------------------------------------
The configuration files used by the UFS Weather Model are listed here and described below:

- *diag_table*
- *field_table*
- *input.nml*
- *model_configure*
- *nems.configure*
- *suite_[suite_name].xml* (used only at build time)


*diag_table* file
------------------------------------
There are three sections in file *diag_table*: Header (Global), File, and Field. These are described below.

**Header Description**

The Header section must reside in the first two lines of the *diag_table* file and contain the title and date
of the experiment (see example below).  The title must be a Fortran character string. The base date is the
reference time used for the time units, and must be greater than or equal to the model start time. The base date
consists of six space-separated integers in the following format:  ``year month day hour minute second``.  Here is an example:

.. code-block:: console

   20161003.00Z.C96.64bit.non-mono
   2016 10 03 00 0 0

**File Description**

The File Description lines are used to specify the name of the file(s) to which the output will be written. They
contain one or more sets of six required and five optional fields (optional fields are denoted by square brackets
``[ ]``).  The lines containing File Descriptions can be intermixed with the lines containing Field Descriptions as
long as files are defined before fields that are to be written them.  File entries have the following format:
 
.. code-block:: console

   "file_name", output_freq, "output_freq_units", file_format, "time_axis_units", "time_axis_name"
   [, new_file_freq, "new_file_freq_units"[, "start_time"[, file_duration, "file_duration_units"]]]

These file line entries are described in :numref:`Table %s <FileDescription>`.

.. _FileDescription:

.. list-table:: *Description of the six required and five optional fields used to define output file sampling rates.*
   :widths: 20 25 55
   :header-rows: 1

   * - File Entry
     - Variable Type
     - Description
   * - file_name
     - CHARACTER(len=128)
     - Output file name without the trailing ".nc"
   * - output_freq
     - INTEGER
     - | The period between records in the file_name:
       |  > 0  output frequency in output_freq_units.
       |  = 0  output frequency every time step (output_freq_units is ignored)
       |  =-1  output at end of run only (output_freq_units is ignored)
   * - output_freq_units
     - CHARACTER(len=10)
     - The units in which output_freq is given.  Valid values are “years”, “months”, “days”, “minutes”, “hours”, or “seconds”.
   * - file_format
     - INTEGER
     - Currently only the netCDF file format is supported.  = 1  netCDF
   * - time_axis_units
     - CHARACTER(len=10)
     - The units to use for the time-axis in the file.  Valid values are “years”, “months”, “days”, “minutes”, “hours”,
       or “seconds”.
   * - time_axis_name
     - CHARACTER(len=128)
     - Axis name for the output file time axis.  The character string must contain the string 'time'.
       (mixed upper and lowercase allowed.)
   * - new_file_freq
     - INTEGER, OPTIONAL
     - Frequency for closing the existing file, and creating a new file in new_file_freq_units.
   * - new_file_freq_units
     - CHARACTER(len=10), OPTIONAL
     - Time units for creating a new file:  either years, months, days, minutes, hours, or seconds.
       NOTE: If the new_file_freq field is present, then this field must also be present.
   * - start_time
     - CHARACTER(len=25), OPTIONAL
     - Time to start the file for the first time.  The format of this string is the same as the global date.    
       NOTE: The new_file_freq and the new_file_freq_units fields must be present to use this field.
   * - file_duration
     - INTEGER, OPTIONAL
     - How long file should receive data after start time in file_duration_units.  This optional field can only be
       used if the start_time field is present.  If this field is absent, then the file duration will be equal to the
       frequency for creating new files.  NOTE: The file_duration_units field must also be present if this field is present.
   * - file_duration_units
     - CHARACTER(len=10), OPTIONAL
     - File duration units. Can be either years, months, days, minutes, hours, or seconds.  NOTE: If the file_duration field
       is present, then this field must also be present.

**Field Description**

The field section of the diag_table specifies the fields to be output at run time.  Only fields registered
with ``register_diag_field()``, which is an API in the FMS ``diag_manager`` routine, can be used in the *diag_table*.
 
Registration of diagnostic fields is done using the following syntax 

.. code-block:: console

   diag_id = register_diag_field(module_name, diag_name, axes, ...)
 
in file ``FV3/atmos_cubed_sphere/tools/fv_diagnostics.F90``.  As an example, the sea level pressure is registered as:
 
.. code-block:: console

   id_slp = register_diag_field (trim(field), 'slp', axes(1:2), &   Time, 'sea-level pressure', 'mb', missing_value=missing_value, range=slprange )
 
All data written out by ``diag_manager`` is controlled via the *diag_table*.  A line in the field section of the
*diag_table* file contains eight variables with the following format:
 
.. code-block:: console

   "module_name", "field_name", "output_name", "file_name", "time_sampling", "reduction_method", "regional_section", packing
 
These field section entries are described in :numref:`Table %s <FieldDescription>`.

.. _FieldDescription:

.. list-table:: *Description of the eight variables used to define the fields written to the output files.*
   :widths: 16 24 55
   :header-rows: 1

   * - Field Entry
     - Variable Type
     - Description
   * - module_name
     - CHARACTER(len=128)
     - Module that contains the field_name variable.  (e.g. dynamic, gfs_phys, gfs_sfc)
   * - field_name
     - CHARACTER(len=128)
     - The name of the variable as registered in the model.
   * - output_name
     - CHARACTER(len=128)
     - Name of the field as written in file_name.
   * - file_name
     - CHARACTER(len=128)
     - Name of the file where the field is to be written.
   * - time_sampling
     - CHARACTER(len=50)
     - Currently not used.  Please use the string "all".
   * - reduction_method
     - CHARACTER(len=50)
     - The data reduction method to perform prior to writing data to disk.  Current supported option is .false..  See ``FMS/diag_manager/diag_table.F90`` for more information.
   * - regional_section
     - CHARACTER(len=50)
     - Bounds of the regional section to capture. Current supported option is “none”. See ``FMS/diag_manager/diag_table.F90`` for more information. 
   * - packing
     - INTEGER
     - Fortran number KIND of the data written.  Valid values:  1=double precision, 2=float, 4=packed 16-bit integers, 8=packed 1-byte (not tested).

Comments can be added to the diag_table using the hash symbol (``#``).
 
A brief example of the diag_table is shown below.  ``“...”`` denote where lines have been removed.

.. code-block:: console

   20161003.00Z.C96.64bit.non-mono
   2016 10 03 00 0 0
 
   "grid_spec",     -1,  "months",   1, "days",  "time"
   "atmos_4xdaily",  6,  "hours",    1, "days",  "time"
   "atmos_static"   -1,  "hours",    1, "hours", "time"
   "fv3_history",    0,  "hours",    1, "hours", "time"
   "fv3_history2d",  0,  "hours",    1, "hours", "time"
 
   #
   #=======================
   # ATMOSPHERE DIAGNOSTICS
   #=======================
   ###
   # grid_spec
   ###
    "dynamics", "grid_lon",  "grid_lon",  "grid_spec", "all", .false.,  "none", 2,
    "dynamics", "grid_lat",  "grid_lat",  "grid_spec", "all", .false.,  "none", 2,
    "dynamics", "grid_lont", "grid_lont", "grid_spec", "all", .false.,  "none", 2,
    "dynamics", "grid_latt", "grid_latt", "grid_spec", "all", .false.,  "none", 2,
    "dynamics", "area",      "area",      "grid_spec", "all", .false.,  "none", 2,
   ###
   # 4x daily output
   ###
    "dynamics",  "slp",       "slp",      "atmos_4xdaily", "all", .false.,  "none", 2
    "dynamics",  "vort850",   "vort850",  "atmos_4xdaily", "all", .false.,  "none", 2
    "dynamics",  "vort200",   "vort200",  "atmos_4xdaily", "all", .false.,  "none", 2
    "dynamics",  "us",        "us",       "atmos_4xdaily", "all", .false.,  "none", 2
    "dynamics",  "u1000",     "u1000",    "atmos_4xdaily", "all", .false.,  "none", 2
    "dynamics",  "u850",      "u850",     "atmos_4xdaily", "all", .false.,  "none", 2
    "dynamics",  "u700",      "u700",     "atmos_4xdaily", "all", .false.,  "none", 2
    "dynamics",  "u500",      "u500",     "atmos_4xdaily", "all", .false.,  "none", 2
    "dynamics",  "u200",      "u200",     "atmos_4xdaily", "all", .false.,  "none", 2
    "dynamics",  "u100",      "u100",     "atmos_4xdaily", "all", .false.,  "none", 2
    "dynamics",  "u50",       "u50",      "atmos_4xdaily", "all", .false.,  "none", 2
    "dynamics",  "u10",       "u10",      "atmos_4xdaily", "all", .false.,  "none", 2

   ...
   ###
   # gfs static data
   ###
    "dynamics",  "pk",        "pk",       "atmos_static",  "all", .false.,  "none", 2
    "dynamics",  "bk",        "bk",       "atmos_static",  "all", .false.,  "none", 2
    "dynamics",  "hyam",     "hyam",      "atmos_static",  "all", .false.,  "none", 2
    "dynamics",  "hybm",     "hybm",       "atmos_static",  "all", .false.,  "none", 2
    "dynamics",  "zsurf",    "zsurf",      "atmos_static",  "all", .false.,  "none", 2
   ###
   # FV3 variables needed for NGGPS evaluation
   ###
   "gfs_dyn",    "ucomp",      "ugrd",     "fv3_history",    "all",  .false.,  "none",  2
   "gfs_dyn",    "vcomp",      "vgrd",     "fv3_history",    "all",  .false.,  "none",  2
   "gfs_dyn",    "sphum",      "spfh",     "fv3_history",    "all",  .false.,  "none",  2
   "gfs_dyn",    "temp",       "tmp",      "fv3_history",    "all",  .false.,  "none",  2
   ...
   "gfs_phys",  "ALBDO_ave",    "albdo_ave", "fv3_history2d", "all", .false., "none",  2
   "gfs_phys",  "cnvprcp_ave",  "cprat_ave", "fv3_history2d", "all", .false., "none",  2
   "gfs_phys",  "cnvprcpb_ave", "cpratb_ave","fv3_history2d", "all", .false., "none",  2
   "gfs_phys",  "totprcp_ave",  "prate_ave", "fv3_history2d", "all", .false., "none",  2
   ...
   "gfs_sfc",   "crain",   "crain",    "fv3_history2d",  "all",  .false.,  "none",  2
   "gfs_sfc",   "tprcp",   "tprcp",    "fv3_history2d",  "all",  .false.,  "none",  2
   "gfs_sfc",   "hgtsfc",  "orog",     "fv3_history2d",  "all",  .false.,  "none",  2
   "gfs_sfc",   "weasd",   "weasd",    "fv3_history2d",  "all",  .false.,  "none",  2
   "gfs_sfc",   "f10m",    "f10m",     "fv3_history2d",  "all",  .false.,  "none",  2
  ...

More information on the content of this file can be found in ``FMS/diag_manager/diag_table.F90``.  

.. note:: None of the lines in the *diag_table* can span multiple lines.

*field_table* file
------------------------------------
The FMS field and tracer managers are used to manage tracers and specify tracer options.  All tracers
advected by the model must be registered in an ASCII table called *field_table*.  The field table consists
of entries in the following format:
 
The first line of an entry should consist of three quoted strings:
 - The first quoted string will tell the field manager what type of field it is. The string ``“TRACER”`` is used to
   declare a field entry. 
 - The second quoted string will tell the field manager which model the field is being applied to.  The supported
   type at present is ``“atmos_mod”`` for the atmosphere model. 
 - The third quoted string should be a unique tracer name that the model will recognize.
 
The second and following lines are called ``methods``.  These lines can consist of two or three quoted strings.
The first string will be an identifier that the querying module will ask for. The second string will be a name
that the querying module can use to set up values for the module. The third string, if present, can supply
parameters to the calling module that can be parsed and used to further modify values. 
 
An entry is ended with a  forward slash (/) as the final character in a row.  Comments can be inserted in the field table by having a hash symbol (#) as the first character in the line.
 
Below is an example of a field table entry for the tracer called ``“sphum”``:

.. code-block:: console

   # added by FRE: sphum must be present in atmos
   # specific humidity for moist runs
    "TRACER", "atmos_mod", "sphum"
              "longname",     "specific humidity"
              "units",        "kg/kg"
              "profile_type", "fixed", "surface_value=3.e-6" /

In this case, methods applied to this tracer include setting the long name to "specific humidity", the units
to "kg/kg". Finally a field named "profile_type" will be given a child field called "fixed", and that field
will be given a field called "surface_value" with a real value of 3.E-6.  The “profile_type” options are listed
in :numref:`Table %s <TracerTable>`.  If the profile type is “fixed” then the tracer field values are set equal
to the surface value.  If the profile type is “profile” then the top/bottom of model and surface values are read
and an exponential profile is calculated, with the profile being dependent on the number of levels in the component model. 

.. _TracerTable:

.. list-table:: *Tracer profile setup from FMS/tracer_manager/tracer_manager.F90.*
   :widths: 20 25 55
   :header-rows: 1

   * - Method Type
     - Method Name
     - Method Control
   * - profile_type
     - fixed
     - surface_value = X
   * - profile_type
     - profile
     - surface_value = X, top_value = Y (atmosphere)

For the case of
 
.. code-block:: console

   "profile_type","profile","surface_value = 1e-12, top_value = 1e-15"
  
in a 15 layer model this would return values of surf_value = 1e-12 and multiplier = 0.6309573,  i.e 1e-15 = 1e-12*(0.6309573^15).

A ``method`` is a way to allow a component module to alter the parameters it needs for various tracers. In essence,
this is a way to modify a default value. A namelist can supply default parameters for all tracers and a method, as
supplied through the field table, will allow the user to modify the default parameters on an individual tracer basis.
The lines in this file can be coded quite flexibly. Due to this flexibility, a number of restrictions are required.
See ``FMS/field_manager/field_manager.F90`` for more information.

*input.nml* file
------------------------------------

The atmosphere model reads many parameters from a Fortran namelist file, named *input.nml*.  This file contains 
several Fortran namelist records, some of which are always required, others of which are only used when selected
physics options are chosen.

The following link describes the various physics-related namelist records:

https://dtcenter.org/GMTB/UFS/sci_doc/CCPPsuite_nml_desp.html

The following link describes the stochastic physics namelist records

https://stochastic-physics.readthedocs.io/en/ufs-v1.0.0/namelist_options.html

The following link describes some of the other namelist records (dynamics, grid, etc):

https://www.gfdl.noaa.gov/wp-content/uploads/2017/09/fv3_namelist_Feb2017.pdf

The namelist section relating to the FMS diagnostic manager is described in the last section of this chapter.

*model_configure* file
------------------------------------

This file contains settings and configurations for the NUOPC/ESMF main component, including the simulation
start time, the processor layout/configuration, and the I/O selections.  :numref:`Table %s <ModelConfigParams>`
shows the following parameters that can be set in *model_configure* at run-time.

.. _ModelConfigParams:

.. list-table:: *Parameters that can be set in model_configure at run-time.*
   :widths: 20 30 15 20
   :header-rows: 1

   * - Parameter
     - Meaning
     - Type
     - Default Value
   * - print_esmf
     - flag for ESMF PET files
     - logical
     - .true.
   * - PE_MEMBER01
     - total number of tasks for ensemble number 1
     - integer
     - 150 (for c96 with quilt)
   * - start_year
     - start year of model integration
     - integer
     - 2019
   * - start_month
     - start month of model integration
     - integer
     - 09
   * - start_day
     - start day of model integration
     - integer
     - 12
   * - start_hour  
     - start hour of model integration 
     - integer 
     - 00
   * - start_minute
     - start minute of model integration
     - integer
     - 0
   * - start_second
     - start second of model integration
     - integer 
     - 0
   * - nhours_fcst
     - total forecast length
     - integer
     - 48
   * - dt_atmos
     - atmosphere time step in second
     - integer
     - 1800 (for C96)
   * - output_1st_tstep_rst
     - output first time step history file after restart
     - logical
     - .false.
   * - memuse_verbose
     - flag for printing out memory usage
     - logical
     - .false.
   * - atmos_nthreads
     - number of threads for atmosphere
     - integer
     - 4
   * - restart_interval
     - frequency to output restart file 
     - integer
     - 0 (write restart file at the end of integration)
   * - quilting 
     - flag to turn on quilt
     - logical
     - .true.
   * - write_groups
     - total number of groups
     - integer
     - 2
   * - write_tasks_per_group
     - total number of write tasks in each write group
     - integer
     - 6
   * - output_history
     - flag to output history files
     - logical
     - .true.
   * - num_files 
     - number of output files 
     - integer
     - 2
   * - filename_base
     - file name base for the output files 
     - character(255)  
     - 'atm' 'sfc'
   * - output_grid 
     - output grid 
     - character(255)
     - gaussian_grid
   * - output_file
     - output file format
     - character(255)
     - nemsio
   * - imo
     - i-dimension for output grid 
     - integer
     - 384
   * - jmo
     - j-dimension for output grid
     - integer
     - 190
   * - nfhout
     - history file output frequency 
     - integer
     - 3
   * - nfhmax_hf
     - forecast length of high history file
     - integer
     - 0 (0:no high frequency output)
   * - nfhout_hf
     - high history file output frequency
     - integer
     - 1
   * - nsout
     - output frequency of number of time step
     - integer
     - -1 (negative: turn off the option, 1: output history file at every time step)

:numref:`Table %s <ModelConfigParamsNotChanged>` shows the following parameters in *model_configure* that
are not usually changed.

.. _ModelConfigParamsNotChanged:

.. list-table:: *Parameters that are not usually changed in model_configure at run-time.*
   :widths: 20 30 15 20
   :header-rows: 1

   * - Parameter
     - Meaning
     - Type
     - Default Value
   * - total_member
     - total number of ensemble member
     - integer
     - 1
   * - RUN_CONTINUE
     - Flag for more than one NEMS run 
     - logical 
     - .false.
   * - ENS_SPS
     - flag for the ensemble stochastic coupling flag 
     - logical
     - .false.
   * - calendar
     - type of calendar year
     - character(*)
     - 'gregorian'
   * - fhrot
     - forecast hour at restart for nems/earth grid component clock in coupled model
     - integer
     - 0
   * - cpl
     - flag for coupling with MOM6/CICE5
     - logical
     - .false.
   * - write_dopost
     - flag to do post on write grid component 
     - logical
     - .false.
   * - ideflate
     - lossless compression level 
     - integer
     - 1 (0:no compression, range 1-9)
   * - nbits
     - lossy compression level
     - integer
     - 14 (0: lossless, range 1-32)
   * - write_nemsioflip
     - flag to flip the vertical level for nemsio file
     - logical
     - .true.
   * - write_fsyncflag
     - flag to check if a file is synced to disk
     - logical
     - .true.
   * - iau_offset
     - IAU offset lengdth 
     - integer
     - 0

*nems.configure* file
------------------------------------
This file contains information about the various NEMS components and their run sequence. In the current release,
this is an atmosphere-only model, so this file is simple and does not need to be changed.  A sample of the file contents is below:

.. code-block:: console

  EARTH_component_list: ATM
  ATM_model:            fv3
  runSeq::
    ATM
  ::

*The SDF (Suite Definition File) file*
---------------------------------------
There are two SDFs currently supported: *suite_FV3_GFS_v15p2.xml* and *suite_FV3_GFS_v16beta.xml*.

=============
Output files
=============

The following files are output when running *fv3.exe* in the default configuration (six files of each kind,
corresponding to the six tiles of the model grid):

- *atmos_4xdaily.tile[1-6].nc*
- *atmos_static.tile[1-6].nc*
- *sfcfHHH.nc*
- *atmfHHH.nc*
- *grid_spec.tile[1-6].nc*

Note that the sfcf* and atmf* files are not output on the 6 tiles, but instead as a single global gaussian grid file.  The specifications of the output files (type, projection, etc) may be overridden in the *model_configure* input file.

Standard output files are *logf???*, and out and err as specified by the job submission. ESMF may also produce log
files (controlled by variable print_esmf in the *model_configure* file), called *PET???.ESMF_LogFile*.

==============================================================
Additional Information about the FMS Diagnostic Manager
==============================================================

The UFS Weather Model output is managed through the FMS (Flexible Modeling System) diagnostic manager (``FMS/diag_manager``)
and is configured using the *diag_table* file. Data can be written at any number of sampling and/or averaging intervals
specified at run-time.  More information about the FMS diagnostic manager can be found at: 
https://data1.gfdl.noaa.gov/summer-school/Lectures/July16/03_Seth1_DiagManager.pdf

------------------------------
Diagnostic Manager namelist
------------------------------
The ``diag_manager_nml`` namelist contains values to control the behavior of the diagnostic manager.   Some
of the more common namelist options are described in :numref:`Table %s <DiagManager>`.  See 
``FMS/diag_manager/diag_manager.F90`` for the complete list.

.. _DiagManager:

.. list-table:: *Namelist variables used to control the behavior of the diagnostic manager.*
   :widths: 15 10 30 10
   :header-rows: 1

   * - Namelist variable
     - Type
     - Description
     - Default value
   * - max_files
     - INTEGER
     - Maximum number of files allowed in diag_table
     - 31
   * - max_output_fields
     - INTEGER
     - Maximum number of output fields allowed in diag_table
     - 300
   * - max_input_fields
     - INTEGER
     - Maximum number of registered fields allowed
     - 300
   * - prepend_date
     - LOGICAL
     - Prepend the file start date to the output file.  .TRUE. is only supported if the diag_manager_init routine is called with the optional time_init parameter.
     - .TRUE.
   * - do_diag_field_log
     - LOGICAL
     - Write out all registered fields to a log file
     - .FALSE.
   * - use_cmor
     - LOGICAL
     - Override the missing_value to the CMOR value of -1.0e20
     - .FALSE.
   * - issue_oor_warnings
     - LOGICAL
     - Issue a warning if a value passed to diag_manager is outside the given range
     - .TRUE.
   * - oor_warnings_fatal
     - LOGICAL
     - Treat out-of-range errors as FATAL
     - .FALSE.
   * - debug_diag_manager
     - LOGICAL
     - Check if the diag table is set up correctly
     - .FALSE.

This release of the UFS Weather Model uses the following namelist:

.. code-block:: console

   &diag_manager_nml
     prepend_date = .false.
   /

