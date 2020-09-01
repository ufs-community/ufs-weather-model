.. _SDFandNamelistExamplePractices:

********************************************
SDF and Namelist Samples and Best Practices
********************************************

The public release of the UFS MR Weather App includes four supported physics suites:
GFS_v15p2, GFS_v15p2_no_nsst, GFS_v16beta, and GFS_v16beta_no_nsst. You will
find the Suite Definition Files (SDFs) for these suites in

https://github.com/NOAA-EMC/fv3atm/tree/ufs-v1.1.0/ccpp/suites

(no other SDFs are available with this release). You will find the namelists for the C96 configuration here:

https://github.com/ufs-community/ufs-weather-model/tree/ufs-v1.1.0/parm/ccpp_v15p2_c96.nml.IN

and

https://github.com/ufs-community/ufs-weather-model/tree/ufs-v1.1.0/parm/ccpp_v16beta_c96.nml.IN

As noted in the file names, these namelists are for the operational (v15p2) and developmental (v16beta)
GFS suites. Each of these namelists are relevant to the suites with and without the SST prediction scheme, that is,
they are relevant for the suite that employs NSST and for the suite that employs the simple ocean
model (`no_nsst`). The only difference in the namelist regarding how SST prediction is
addressed is variable `nstf_name`. For more information about this variable and for information about
namelist options for higher resolution configurations, please consult the
`CCPP v4.1.0 Scientific Documentation <https://dtcenter.org/GMTB/v4.1.0/sci_doc/>`_.

The four CCPP suites for the UFS MR Weather App release are supported in four grid resolutions:
C96, C192, C384, and C768, with 64 vertical levels.

An in depth description of the namelist settings, SDFs, and parameterizations used
in all supported suites can be found in the `CCPP v4.1.0 Scientific Documentation <https://dtcenter.org/GMTB/v4.1.0/sci_doc/>`_.
Note both suites do not
use stochastic physics by default, but the stochastic physics can be activated following the
instructions described in the `stochastic physics v1.1 user's guide <https://stochastic-physics.readthedocs.io/en/release-ufs-v1.1.0/>`_.

Both the SDF and the *input.nml* contain information about how to specify the physics suite.
Some of this information is redundant, and the user must make sure they are compatible. The
safest practice is to use the SDF and namelist provided for each suite, since those are
supported configurations.

Changes to the SDF must be accompanied by corresponding changes to the namelist. While there
is not a one-to-one correspondence between the namelist and the SDF, :numref:`Table %s <PBLVarOptions>`
shows some variables in the namelist that must match the SDF.

.. _PBLVarOptions:

.. list-table:: *Variables related to PBL options*
   :widths: 15 30 10 10 20 30
   :header-rows: 1

   * - Namelist option
     - Meaning
     - Possible Values
     - Default
     - Used with CCPP scheme
     - Recommentation
   * -
     - **PBL-related variables**
     -
     -
     -
     -
   * - do_myjpbl
     - Flag to activate the MYJ PBL scheme
     - T
     - F
     - mypbl_wrapper
     - Set to F for GFSv15p2* and GFSv16beta*
   * - do_myjsfc
     - Flag to activate the MYJ PBL surface layer scheme
     - T, F
     - F
     - myjsfc_wrapper
     - Set to F for GFSv15p2* and GFSv16beta*
   * - do_mynnedmf
     - Flag to activate the MYNN-EDMF scheme
     - T, F
     - F
     - mynnedmf_wrapper
     - Set to F for GFSv15p2* and GFSv16beta*
   * - do_ysu
     - Flag to activate the YSU PBL scheme
     - T, F
     - F
     - ysudif
     - Set to F for GFSv15p2* and GFSv16beta*
   * - hybedmf
     - Flag to activate the K-based PBL scheme
     - T, F
     - F
     - hedmf
     - Set to T for GFSv15p2* and GFSv16beta*
   * - isatedmf
     - Flag for version of scale-aware TKE-based EDMF scheme
     - 0, 1
     - 0
     - 0=satmedmfvdif, 1=satmedmfvdifq
     - Set to 0 for GFSv15p2* and 1 for GFSv16beta*
   * - ism
     - Flag to choose a land surface model to use
     - 0, 1, 2
     - 1
     - 1=lsm_noah, 2=lsm_ruc
     - Set to 1 for GFSv15p2* and GFSv16beta*
   * - satedmf
     - Flag to activate the scale-aware TKE-based EDMF scheme
     - T, F
     - F
     - satmedmfvdif or satmedmfvdifq
     - Set to T for GFSv15p2* and GFSv16beta*
   * - shinhong
     - Flag to activate the Shin-Hong PBL parameterization
     - T, F
     - F
     - shinhongdif
     - Set to F for GFSv15p2* and GFSv16beta*
   * -
     - **Convection-releated flags**
     -
     -
     -
     -
   * - cscnv
     - Flag to activate the Chikira-Sugiyama deep convection scheme
     - T, F
     - F
     - cs_conv
     - Set to F for GFSv15p2* and GFSv16beta*
   * - do_aw
     - Flag to activate the Arakawa-Wu extension to the Chikira-Sugiyama deep convection scheme
     - T, F
     - F
     - cs_conv_aw_adj
     - Set to F for GFSv15p2* and GFSv16beta*
   * - imfdeepcnv
     - Flag to choose a mass flux deep convective scheme
     - -1, 2, 3, 4
     - -1
     - -1=no deep convection*, 2=samfshalcnv, 3=cu_gf_driver, 4=cu_ntiedtke
     - Set to 2 for GFSv15p2* and GFSv16beta*
   * - imfshalcvn
     - Flag to choose a mass flux shallow convective scheme
     - -1, 2, 3, 4
     - -1
     - -1=no deep convection*, 2=samfshalcnv, 3=cu_gf_driver, 4=cu_ntiedtke
     - Set to 2 for GFSv15p2* and GFSv16beta*

\*Even when imfdeepcvn=-1, the Chikira-Sugiyama deep convection scheme may be specified using cscnv=T.

Other miscellaneous changes to the SDF that must be accompanied by corresponding changes in
the namelist are listed in :numref:`Table %s <MiscVarOptions>`.

.. _MiscVarOptions:

.. list-table:: *Miscellaneous namelist variables and their relation to the SDF*
   :widths: 15 30 10 10 20 30
   :header-rows: 1

   * - Namelist option
     - Meaning
     - Possible Values
     - Default
     - Used with CCPP scheme
     - Recommendation
   * -
     - **Miscellaneous variables**
     -
     -
     -
     -
   * - do_myjsfc
     - Flag to activate the MYJ PBL surface scheme
     - T, F
     - F
     - mynnsfc_wrapper
     - Set to F for GFSv15p2* and GFSv16beta*
   * - do_shoc
     - Flag to activate the Simplified Higher-Order Closure (SHOC) parameterization
     - T, F
     - F
     - shoc
     - Set to F for GFSv15p2* and GFSv16beta*
   * - do_ugwp**
     - Flag to activate the unified Gravity Wave Physics parameterization
     - T, F
     - F
     - cires_ugwp
     - Set to F for GFSv15p2* and GFSv16beta*
   * - imp_physics
     - Flag to choose a microphysics scheme
     - 8, 10, 11
     - 99
     - 8=mp_thompson, 10=m_micro, 11=gfdl_cloud_microphysics
     - Set to 11 for GFSv15p2* and GFSv16beta*
   * - lsm
     - Flag to choose a land surface model to use
     - 0, 1, 2
     - 1
     - 1=lsm_noah, 2=lsm_ruc
     - Set to 1 for GFSv15p2* and GFSv16beta*
   * - lsoil
     - Number of soil layers
     - 4, 9
     - 4
     - 4 for lsm_noah, 9 for lsm_ruc
     - Set to 4 for GFSv15p2* and GFSv16beta*
   * - h2o_phys
     - Flag for stratosphere h2o scheme
     - T, F
     -
     - h2ophys
     - Set to T for GFSv15p2* and GFSv16beta*
   * - oz_phys_2015
     - Flag for new (2015) ozone physics
     - T, F
     -
     - ozphys_2015
     - Set to T for GFSv15p2* and GFSv16beta*

\*\*The CIRES Unified Gravity Wave Physics (cires_ugwp) scheme is used in GFSv15p2* and GFSv16beta* SDFs with do_ugwp=F in the namelist. In this setting, the cires_ugwp calls the operational GFS v15.2 orographic gravity wave drag (gwdps) scheme. When do_ugwp=T, the  cires_ugwp scheme calls an experimental orographic gravity wave (gwdps_v0).

**Note that some namelist variables are not available for use with CCPP.**

   * **do_deep**. In order to disable deep convection, it is necessary to remove the deep convection scheme from the SDF.
   * **shal_cnv**. In order to disable shallow convection, it is necessary to remove the deep convection scheme from the SDF.
   * **ldiag3d** and **ldiag_ugwp**. Must be F for CCPP runs.
   * **gwd_opt**. Ignored in CCPP-supported suites.

**When certain parameterizations are turned on, additional namelist options can be used (they are ignored otherwise).
Some examples are shown in** :numref:`Table %s <EnabledNMLOptions>`.

.. _EnabledNMLOptions:

.. list-table:: *Enabled namelist variables*
   :widths: 10 50
   :header-rows: 1

   * - Namelist setting
     - Enabled namelist variables
   * - do_ugwp=T
     - All variables in namelist record &cires_ugwp_nml plus do_tofd. Additionally, if namelist variable cnvgwd=T and
       the third and fourth position of namelist array cdmbgwd are both 1, then the convective gravity wave drag that
       is part of cires_ugwp will be called. (Not supported with the UFS)
   * - do_mynnedmf=T
     - bl_mynn_tkeadvect, bl_mynn_edmf, bl_mynn_edmf_mom (Not supported with the UFS)
   * - imp_physics=99
     - psautco and prautco (Not supported with the UFS)
   * - imp_physics=10
     - mg_* (Not supported with UFS)
   * - imp_physics=11
     - All variables in namelist record gfdl_cloud_microphysics_nml and lgfdlmprad
   * - satedmf=T
     - isatedmf
