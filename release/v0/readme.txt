---------------------------------------------------------------------
Instruction for Building and Running FV3GFS Forecast Experiments  
---------------------------------------------------------------------
Revision History
04/30/2017 -- Fanglin Yang, First version.


Note: The workflow has only been tested on WCOSS Cray.  Options need to be added 
      and tested for running the model on Theia and other machines

1.  check out  https://svnemc.ncep.noaa.gov/projects/nems/apps/NEMSfv3gfs/trunk, 
    and use ?? to compile the model.  
    Four executable files will be created under ./trunk/NEMS/exe
      fv3_gfs_hydro.prod.32bit.x*  
      fv3_gfs_hydro.prod.64bit.x*  
      fv3_gfs_nh.prod.32bit.x*  
      fv3_gfs_nh.prod.64bit.x

2. check out https://svnemc.ncep.noaa.gov/projects/fv3gfs/tags/FV3GFS_V0_RELEASE
   use ./FV3GFS_V0_RELEASE/sorc/fre-nctools.fd/BUILD_TOOLS.csh $system_site to 
   compile regrid utilities, where system_site can be cray, theia or gaea 

3. Copy the directory of fixed fields (binary files) from WCOSS Surge
   cp -rp /gpfs/hps/emc/global/noscrub/Fanglin.Yang/svn/FV3GFS_V0_RELEASE/fix to ./FV3GFS_V0_RELEASE/fix 
   
4. Copy the directory of initial conditions from WCOSS Surge
   cp -rp /gpfs/hps/emc/global/noscrub/Fanglin.Yang/svn/FV3GFS_V0_RELEASE/ICs to ./FV3GFS_V0_RELEASE/ICs 
   It contains cold start initial conditions for the following three cases. They are converted from 
   gfs analysese of Q3FY17 NEMS GFS retrospective parallels.
     2016092900  Hurricane Matthew
     2016011812  Winter East Coast Blizzard 
     2016081200  Louisiana Flooding 

 
5.  Check ./FV3GFS_V0_RELEASE/exp/runjob_cray.sh to see how to submit forecast batch jobs on WCOSS Cray.
    Modify the following two lines in runjob_cray.sh to select model resolution and case to run
      export CASE=C384                  ;#resolution, C96 (~100km), C382 (~25km) or C768 (~13km)
      export CDATE=2016092900           ;#initial condition dates  2016092900 2016011812 2016081200

    Users need to make sure the following directories are correctly set up
     export BASE_GSM=/gpfs/hps/emc/global/noscrub/$LOGNAME/svn/FV3GFS_V0_RELEASE   ;# source directory
     export FIX_FV3=$BASE_GSM/fix/fix_fv3                  ;#model fixed fields
     export IC_DIR=$BASE_GSM/ICs                           ;#forecast initial conditions

     # temporary running directory
       export DATA=/gpfs/hps/stmp/$LOGNAME/${CASE}${PSLOT}${CDATE}

     # directory to save output
       export ROTDIR=/gpfs/hps/ptmp/$LOGNAME/$PSLOT/$CASE

     # NEMS FV3GFS forecast executable directory
       export FCSTEXECDIR=/gpfs/hps/emc/global/noscrub/Fanglin.Yang/svn/fv3gfs/NEMSfv3gfs/trunk/NEMS/exe

     Use "bsub<runjob_cray.sh" to submit batch job on WCOSS Cray.

     The default option is to run the model with non-hydrostatic 32-bit fv3 dycore. 
     Please check runjob_cray.sh to see other options.  Inline documentations are included.

6.  Forecast output saved in $ROTDIR are on global lat-lon 0.25-degree grid in netCDF format. 

   

