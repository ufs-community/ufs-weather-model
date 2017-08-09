---------------------------------------------------------------------
Instructions for Building and Running FV3GFS v0 release Forecast Experiments
---------------------------------------------------------------------

Section 1: How to get code, compile and run experiment
======================================================

Note: The workflow has only been tested on WCOSS Cray, theia and jet.

1.  Check out release version at:
--------------------------------
      https://svnemc.ncep.noaa.gov/projects/nems/apps/NEMSfv3gfs/tags/fv3gfs.v0release

     % svn co https://svnemc.ncep.noaa.gov/projects/nems/apps/NEMSfv3gfs/tags/fv3gfs.v0release
     % cd fv3gfs.v0release

     Code is also accessible via git at NCEP's VLab.

     Initiate VLab account creation using your noaa.gov credentials via 
     https://vlab.ncep.noaa.gov/redmine/

     To request access to the fv3gfs git repo, fill in the web form located here
     https://vlab.ncep.noaa.gov/group/fv3gfs/home

     Clone the fv3gfs repo
     % git clone https://user.name@vlab.ncep.noaa.gov/git/comfv3
     % cd comfv3
     % git submodule init
     % git submodule update


2.  Compile source code:
------------------------
    a) compile the nems fv3 forecast model.

     % cd release/v0/exp
     % ./build.sh machine_name(wcoss_cray, theia, or jet)
      Four executable files will be created under fv3gfs.v0release/NEMS/exe:
        fv3_gfs_hydro.prod.32bit.x*
        fv3_gfs_hydro.prod.64bit.x*
        fv3_gfs_nh.prod.32bit.x*
        fv3_gfs_nh.prod.64bit.x
      It will take about one hour to finish compiling.

   b) compile remap source code

     % cd ../sorc/fre-nctools.fd
      now the current directory is: fv3gfs.v0release/release/v0/sorc/fre-nctools.fd
     % ./BUILD_TOOLS.csh machine_name (wcoss_cray, theia or jet)
      Six executables will be located under:  fv3gfs.v0release/release/v0/exec
          filter_topo
          fregrid
          fregrid_parallel
          make_hgrid
          make_hgrid_parallel
          make_solo_mosaic

3. Run a default experiment:
----------------------------
   The default setting is to run an experiment with non-hydrostatic 32-bit fv3 dycore
   (executable fv3_gfs_nh.prod.32bit.x with C96 resolution, output are on 6 tile cubed-sphere grid,
    the outputs are then remapped to 1 degree global lat-lon grid in netcdf format)

     % cd ../../exp
       now the current directory is: fv3gfs.v0release/release/v0/exp

        on wcoss cray:
     % bsub < runjob_cray.sh
       job error file: err_cray, output print file: out_cray
       run directory is at: /gpfs/hps/stmp/$LOGNAME/C96fv3gfs2016092900
       results are saved at: /gpfs/hps/ptmp/$LOGNAME/fv3gfs/C96/gfs.20160929/00

        on theia:
       If you do not belong to nems group, please change the job account from nems to your account,
        change runjob_theia.sh line 5:

          #PBS -A nems

        to find your account on theia, please issue following command in any of your local directory:
          % account_params
       If you are Non-NCEPDEV people, please change directories DATA and ROTDIR to an area they can write to
        change runjob_theia.sh line 41 and 45:
          export DATA=/scratch4/NCEPDEV/stmp3/$LOGNAME/${CASE}${PSLOT}${CDATE}
          export ROTDIR=/scratch4/NCEPDEV/stmp3/$LOGNAME/$PSLOT/$CASE


       to submit a job on theia:
     % qsub runjob_theia.sh
       job error file: err_theia, output print file: out_theia
       run directory is at: /scratch4/NCEPDEV/stmp3/$LOGNAME/C96fv3gfs2016092900
       results are saved at: /scratch4/NCEPDEV/stmp3/$LOGNAME/fv3gfs/C96/gfs.20160929/00

4. Compare with baseline:
-------------------------
   A baseline for default experiment is set up on supported platforms. To compare with baseline:

     % cd .
       current directory is now: fv3gfs.v0release/release/v0/exp
     % ./diff_baseline.sh machine_name(wcoss_cray, theia, or jet) case(C96, C384, C768)
       default: ./diff_baseline.sh will compare results on wcoss_cray for C96

     on theia:
        If you are Non-NCEPDEV people, please change directory dir1 to location of their output
        change diff_baseline.sh line 13
          dir1=/scratch4/NCEPDEV/stmp3/${LOGNAME}/fv3gfs/${CASE}/gfs.20160929/00

   A message will be given at the end of the script.

   on wcoss cray, the baseline is at:
   /gpfs/hps/nems/noscrub/emc.nemspara/FV3GFS_V0_RELEASE/baseline/fv3gfs_nh_32bit/${CASE}/gfs.20160929/00

   on theia, the baseline is at:
   /scratch4/NCEPDEV/nems/noscrub/emc.nemspara/FV3GFS_V0_RELEASE/baseline/fv3gfs_nh_32bit/${CASE}/gfs.20160929/00


Section 2: Questions and Answers:
=================================

Q1: Where are the fix files?
A1: The fix files (binary files) are located on WCOSS Surge/Luna, Theia and Jet
   a) wcoss_cray:
   /gpfs/hps/emc/nems/noscrub/emc.nemspara/FV3GFS_V0_RELEASE/fix
   b) theia:
   /scratch4/NCEPDEV/nems/noscrub/emc.nemspara/FV3GFS_V0_RELEASE/fix

Q2: Where are the initial conditions:
A2: Initial conditions from WCOSS Surge are in a fixed location on WCOSS Surge/Luna, Theia and Jet
   a) wcoss_cray:
    /gpfs/hps/emc/nems/noscrub/emc.nemspara/FV3GFS_V0_RELEASE/ICs
   b) theia:
    /scratch4/NCEPDEV/nems/noscrub/emc.nemspara/FV3GFS_V0_RELEASE/ICs

   The initial conditions contain cold start initial conditions for the following three cases.
     2016092900  Hurricane Matthew
     2016011812  Winter East Coast Blizzard
     2016081200  Louisiana Flooding
   They are converted from gfs analysese of Q3FY17 NEMS GFS retrospective parallels.

Q3: How to run a non-default experiment?
A4: First view the job card to see how to submit forecast batch jobs on WCOSS Cray.
      https://svnemc.ncep.noaa.gov/projects/nems/apps/NEMSfv3gfs/tags/fv3gfs.v0release/release/v0/exp
      runjob_cray.sh    runjob_theia.sh   runjob_jet.sh

    a) run different resolution:
      export CASE=C96                   # resolution: C96 (~100km), C384 (~25km) or C768 (~13km)
      change CASE to  C384 or  C768
      on wcoss_cray, only CASE needs to be changed.
      on Theia, please add change line 8 in runjob_theia.sh:
        for C384, change line 8 to:
        #PBS -l nodes=96:ppn=12
        for C768, change line 8 to:
        #PBS -l nodes=192:ppn=12

    b) run different case:
      export CDATE=2016092900            # initial condition dates:  2016092900 2016011812 2016081200
      change default CDATE from 2016092900 to 2016011812 or 2016081200


    c) Run different executable:
       export MODE=32bit                 # dycore precision:   32bit, 64bit
       export TYPE=nh                    # hydrostatic option: nh, hydro
       change MODE to 64bit or TYPE to hydro, you will get different forecast executable, please see section 2a.


    d) Remap to different global lat-lon resolution:
       export master_grid=1deg           # 1deg 0p5deg 0p25deg 0p125deg

    e) run with different thread number
       export nth_f=2                    # number of threads: 1, 2, 4
       On theia, if number of thread is changed, please also change line 8,
       the value of ppn needs to be equal to task_per_node

    f) run with different nodes
       export NODES=16                   # number of nodes for forecast job
       export REMAP_TASKS=32             # number of mpi tasks for remap job
       If change forecast job node, please follow: NODES*task_per_node = layout_x*layout_y*6
       If change remap tasks, REMAP_TASKS must be less than the number of latitude

     After changes are made to runjobs_${machine}.sh, please following example in section 3 to submit a job.

Q4. Is there Any other variable that can be changed in the run script?

     Following variables in the job script runjobs_${machine}.sh can be modified:

     a) temporary running directory
        export DATA=/gpfs/hps/stmp/$LOGNAME/${CASE}${PSLOT}${CDATE}

     b) directory to save output
        export ROTDIR=/gpfs/hps/ptmp/$LOGNAME/$PSLOT/$CASE

     c) NEMS FV3GFS forecast executable directory(if forecast executable is moved to a different location,
        the FCSTEXECDIR needs to point to that location)
        export FCSTEXECDIR=/gpfs/hps/emc/global/noscrub/Fanglin.Yang/svn/fv3gfs/NEMSfv3gfs/trunk/NEMS/exe

     Please check runjob_${machine}.sh to see other options. Inline documentations are included.

Q5. What nceplibs libraries are used and how to get these libraries?
A5: Following nceplibs libraries are used in this release:

      bacio-intel/2.0.1
      ip-intel/2.0.0
      sp-intel/2.0.2
      w3nco-intel/2.0.6
      w3emc-intel/2.2.0

    nceplibs source code can be found at: https://svnemc.ncep.noaa.gov/trac/nceplibs


Q6. What other external libraries are used?
A6: External libraries used in this release include:

     a) intel/mpi
     Intel 14-15

     b) netcdf4.3.2
     download and installation: http://www.unidata.ucar.edu/blogs/news/entry/netcdf_4_3_2
     questions: support@unidata.ucar.edu

     c) ESMF7.0.0:
     download and installation: https://www.earthsystemcog.org/news/detail/126/
     questions: esmf_support@list.woc.noaa.gov

Q7: Where to get help for fv3 dycore questions?
A7: GFDL fv3 team at: oar.gfdl.fvgfs_support@noaa.gov

Q8: Where to get help for gfs physics  questions?
A8: GFS physics team

