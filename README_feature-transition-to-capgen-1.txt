Branches feature/transition-to-capgen-1 where created for the following repositories
on 09/19/2020 with these original hashes:

ufs-weather-model
-----------------

commit 2a6528d0442befda2914e7aaaa67ce7794681e09 (HEAD -> develop, noaa-emc/develop, noaa-emc/HEAD)
Author: Jessica Meixner <jessica.meixner@noaa.gov>
Date:   Thu Sep 17 12:34:23 2020 -0400

    update to latest WW3 (#199)

    * update to latest WW3,
    * update to latest NEMS and stochastic physics branches
    * point to new baseline in rt.sh
    * wcoss_dell_p3 log, wcoss_cray log, hera.intel log,* orion.intel log

submodules
----------

f61416fef691d9ba39a40df1ce72aa574f54c390 FMS (2019.01.03)
6bc61df3c363f9134a46439ff4a5a4a803daafb1 FV3 (heads/feature/transition-to-capgen-1)
   8b59ebc039dafe1c20ed6dd21cb38ca564852b98 atmos_cubed_sphere (heads/dev/emc)
   f06e053db04eaea602d43d6221081ba54fb6cc95 ccpp/framework (ccpp_transition_to_vlab_master_20190705-213-g1dda224)
   0808cc2e8938ba66003b46746858143a9d75addb ccpp/physics (ccpp_transition_to_vlab_master_20190705-771-g0808cc2) 
      566bee9cd6f9977e82d75d9b4964b20b1ff6163d physics/rte-rrtmgp (1.2.1-22-g566bee9)
9d05172b711f4ab5d6f978dbe575bd67a681b55a NEMS (heads/develop)
   ca171b95095db4fcd0fc7b01c23d073d90becd99 tests/produtil/NCEPLIBS-pyprodutil (ufs-v1.1.0-18-gca171b9)
96e3f3a8fa0389a4b110b0fa23e7a414f6d92038 WW3 (6.07.1-50-g96e3f3a8)
ffdd19bc6c1df747394b7e9958a76238fcd44242 stochastic_physics (ufs-v1.0.0-70-gffdd19b)

All feature/transition-to-capgen-1 are under the NCAR GitHub account and are protected.


Change log
----------

09/19/2020: remove .gitmodules from ccpp-framework, change .gitmodules in fv3atm and ufs-weather-model,
            add CODEOWNERS to fv3atm and ufs-weather-model

10/05/2020: merge https://github.com/NCAR/ccpp-physics/pull/498, remove other empty CCPP arg table lines
            in ccpp-physics that do not have entries in the metadata files (empty subroutines), replace
            horizontal_dimension with horizontal_loop_extent in GFS_typedefs.F90

