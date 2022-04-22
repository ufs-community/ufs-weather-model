!>@brief The module 'stochy_data_mod' contains the initilization routine that read the stochastic phyiscs
!! namelist and determins the number of random patterns.
module stochy_data_mod

! set up and initialize stochastic random patterns.

 use spectral_transforms, only: len_trie_ls,len_trio_ls,ls_dim,ls_max_node,&
                              skeblevs,levs,jcap,lonf,latg,initialize_spectral
 use stochy_namelist_def
 use constants_mod, only : radius
 use mpi_wrapper, only: mp_bcst, is_rootpe, mype
 use stochy_patterngenerator_mod, only: random_pattern, patterngenerator_init,&
 getnoise, patterngenerator_advance,ndimspec,chgres_pattern,computevarspec_r
 use stochy_internal_state_mod
! use mersenne_twister_stochy, only : random_seed
 use mersenne_twister, only : random_seed
 use compns_stochy_mod, only : compns_stochy

 implicit none
 private
 public :: init_stochdata,init_stochdata_ocn

 type(random_pattern), public, save, allocatable, dimension(:) :: &
       rpattern_sppt,rpattern_shum,rpattern_skeb, rpattern_sfc,rpattern_epbl1,rpattern_epbl2,rpattern_ocnsppt,rpattern_spp
 integer, public :: nepbl=0
 integer, public :: nocnsppt=0
 integer, public :: nsppt=0
 integer, public :: nshum=0
 integer, public :: nskeb=0
 integer, public :: nlndp=0 ! this is the number of different patterns (determined by the tau/lscale input) 
 integer, public :: nspp =0 ! this is the number of different patterns (determined by the tau/lscale input) 
 real*8, public,allocatable :: sl(:)

 real(kind=kind_dbl_prec),public, allocatable :: vfact_sppt(:),vfact_shum(:),vfact_skeb(:),vfact_spp(:)
 real(kind=kind_dbl_prec),public, allocatable :: skeb_vwts(:,:)
 integer                 ,public, allocatable :: skeb_vpts(:,:)
 real(kind=kind_dbl_prec),public, allocatable :: gg_lats(:),gg_lons(:)
 real(kind=kind_dbl_prec),public :: wlon,rnlat,rad2deg
 real(kind=kind_dbl_prec),public, allocatable :: skebu_save(:,:,:),skebv_save(:,:,:)
 integer,public :: INTTYP
 type(stochy_internal_state),public :: gis_stochy,gis_stochy_ocn

 contains
!>@brief The subroutine 'init_stochdata' determins which stochastic physics
!!pattern genertors are needed.
!>@details it reads the nam_stochy namelist and allocates necessary arrays
 subroutine init_stochdata(nlevs,delt,input_nml_file,fn_nml,nlunit,iret)
!\callgraph

! initialize random patterns.  
   use netcdf
   implicit none
   integer, intent(in) :: nlunit,nlevs
   character(len=*),  intent(in) :: input_nml_file(:)
   character(len=64), intent(in) :: fn_nml
   real, intent(in) :: delt
   integer, intent(out) :: iret
   real :: ones(5)

   real :: rnn1
   integer :: nn,k,nm,stochlun,ierr,n
   integer :: locl,indev,indod,indlsod,indlsev
   integer :: l,jbasev,jbasod
   integer :: jcapin,varid1,varid2
   real(kind_dbl_prec),allocatable :: noise_e(:,:),noise_o(:,:)
   include 'function_indlsod'
   include 'function_indlsev'
   include 'netcdf.inc'
   stochlun=99
   levs=nlevs

   iret=0
! read in namelist

   call compns_stochy (mype,size(input_nml_file,1),input_nml_file(:),fn_nml,nlunit,delt,iret)
  
   if (iret/=0) return  ! need to make sure that non-zero irets are being trapped.
   if ( (.NOT. do_sppt) .AND. (.NOT. do_shum) .AND. (.NOT. do_skeb)  .AND. (lndp_type==0) .AND. (.NOT. do_spp)) return

   call initialize_spectral(gis_stochy)
   allocate(noise_e(len_trie_ls,2),noise_o(len_trio_ls,2))
! determine number of random patterns to be used for each scheme.
   do n=1,size(sppt)
     if (sppt(n) > 0) then
        nsppt=nsppt+1
     else
        exit
     endif
   enddo
   if (is_rootpe()) print *,'nsppt = ',nsppt
   do n=1,size(shum)
     if (shum(n) > 0) then
        nshum=nshum+1
     else
        exit
     endif
   enddo
   if (is_rootpe()) print *,'nshum = ',nshum
   do n=1,size(skeb)
     if (skeb(n) > 0) then
        nskeb=nskeb+1
     else
        exit
     endif
   enddo
   if (is_rootpe()) print *,'nskeb = ',nskeb
   ! Draper: nlndp>1 was not properly coded. Hardcode to 1 for now
   !do n=1,size(lndp_z0)
   !  if (lndp_z0(n) > 0 .or. lndp_zt(n)>0 .or. lndp_hc(n)>0 .or. &
   !      lndp_vf(n)>0 .or. lndp_la(n)>0 .or. lndp_al(n)>0) then
   !     nlndp=nlndp+1
   !  else
   !     exit
   !  endif
   !enddo
   if (n_var_lndp>0) nlndp=1
   if (n_var_spp>0) nspp=n_var_spp
   if (is_rootpe())  print *,' nlndp   = ', nlndp
   if (is_rootpe())  print *,' nspp   = ', nspp

   if (nsppt > 0) allocate(rpattern_sppt(nsppt))
   if (nshum > 0) allocate(rpattern_shum(nshum))
   if (nskeb > 0) allocate(rpattern_skeb(nskeb))
   ! mg, sfc perts
   if (nlndp > 0) allocate(rpattern_sfc(nlndp))
   if (nspp > 0) allocate(rpattern_spp(nspp))

!  if stochini is true, then read in pattern from a file
   if (is_rootpe()) then
      if (stochini) then
         print*,'opening stoch_ini'
         !OPEN(stochlun,file='INPUT/atm_stoch.res.bin',form='unformatted',iostat=ierr,status='old')
         ierr=nf90_open('INPUT/atm_stoch.res.nc',nf90_nowrite,ncid=stochlun)
         if (ierr .NE. 0) then
            write(0,*) 'error opening stoch_ini'
            iret = ierr
            return
         end if
         ierr=NF90_GET_ATT(stochlun,NF_GLOBAL,"ntrunc",jcapin)
         if (ierr .NE. 0) then
            write(0,*) 'error getting ntrunc'
            iret = ierr
            return
         end if
         print*,'ntrunc read in',jcapin
      endif
   endif
   ! no spinup needed if initial patterns are defined correctly.
   if (nsppt > 0) then
      if (is_rootpe()) then
         print *, 'Initialize random pattern for SPPT'
         if (stochini) then
            ierr=NF90_INQ_VARID(stochlun,"sppt_seed", varid1)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring SPPT seed'
               iret = ierr
               return
            end if
            ierr=NF90_INQ_VARID(stochlun,"sppt_spec", varid2)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring SPPT spec'
               iret = ierr
               return
            end if
         endif
      endif
     print*,'calling init',lonf,latg,jcap
      call patterngenerator_init(sppt_lscale(1:nsppt),spptint,sppt_tau(1:nsppt),sppt(1:nsppt),iseed_sppt,rpattern_sppt, &
           lonf,latg,jcap,gis_stochy%ls_node,nsppt,1,0,new_lscale)
      do n=1,nsppt
         if (stochini) then
            call read_pattern(rpattern_sppt(n),jcapin,stochlun,1,n,varid1,varid2,.false.,ierr)
            if (ierr .NE. 0) then
               write(0,*) 'error reading SPPT pattern'
               iret = ierr
               return
            end if
         else
            call getnoise(rpattern_sppt(n),noise_e,noise_o)
            do nn=1,len_trie_ls
               rpattern_sppt(n)%spec_e(nn,1,1)=noise_e(nn,1)
               rpattern_sppt(n)%spec_e(nn,2,1)=noise_e(nn,2)
               nm = rpattern_sppt(n)%idx_e(nn)
               if (nm .eq. 0) cycle
               rpattern_sppt(n)%spec_e(nn,1,1) = rpattern_sppt(n)%stdev*rpattern_sppt(n)%spec_e(nn,1,1)*rpattern_sppt(n)%varspectrum(nm)
               rpattern_sppt(n)%spec_e(nn,2,1) = rpattern_sppt(n)%stdev*rpattern_sppt(n)%spec_e(nn,2,1)*rpattern_sppt(n)%varspectrum(nm)
            enddo
            do nn=1,len_trio_ls
               rpattern_sppt(n)%spec_o(nn,1,1)=noise_o(nn,1)
               rpattern_sppt(n)%spec_o(nn,2,1)=noise_o(nn,2)
               nm = rpattern_sppt(n)%idx_o(nn)
               if (nm .eq. 0) cycle
               rpattern_sppt(n)%spec_o(nn,1,1) = rpattern_sppt(n)%stdev*rpattern_sppt(n)%spec_o(nn,1,1)*rpattern_sppt(n)%varspectrum(nm)
               rpattern_sppt(n)%spec_o(nn,2,1) = rpattern_sppt(n)%stdev*rpattern_sppt(n)%spec_o(nn,2,1)*rpattern_sppt(n)%varspectrum(nm)
            enddo
         endif
       enddo
   endif
   if (nshum > 0) then
      if (is_rootpe()) then
         print *, 'Initialize random pattern for SHUM'
         if (stochini) then
            ierr=NF90_INQ_VARID(stochlun,"shum_seed", varid1)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring SHUM seed'
               iret = ierr
               return
            end if
            ierr=NF90_INQ_VARID(stochlun,"shum_spec", varid2)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring SHUM spec'
               iret = ierr
               return
            end if
         endif
      endif
      call patterngenerator_init(shum_lscale(1:nshum),shumint,shum_tau(1:nshum),shum(1:nshum),iseed_shum,rpattern_shum, &
          lonf,latg,jcap,gis_stochy%ls_node,nshum,1,0,new_lscale)
      do n=1,nshum
         if (stochini) then
            call read_pattern(rpattern_shum(n),jcapin,stochlun,1,n,varid1,varid2,.false.,ierr)
            if (ierr .NE. 0) then
               write(0,*) 'error reading SHUM pattern'
               iret = ierr
               return
            end if
          else
             call getnoise(rpattern_shum(n),noise_e,noise_o)
             do nn=1,len_trie_ls
                rpattern_shum(n)%spec_e(nn,1,1)=noise_e(nn,1)
                rpattern_shum(n)%spec_e(nn,2,1)=noise_e(nn,2)
                nm = rpattern_shum(n)%idx_e(nn)
                if (nm .eq. 0) cycle
                rpattern_shum(n)%spec_e(nn,1,1) = rpattern_shum(n)%stdev*rpattern_shum(n)%spec_e(nn,1,1)*rpattern_shum(n)%varspectrum(nm)
                rpattern_shum(n)%spec_e(nn,2,1) = rpattern_shum(n)%stdev*rpattern_shum(n)%spec_e(nn,2,1)*rpattern_shum(n)%varspectrum(nm)
             enddo
             do nn=1,len_trio_ls
                rpattern_shum(n)%spec_o(nn,1,1)=noise_o(nn,1)
                rpattern_shum(n)%spec_o(nn,2,1)=noise_o(nn,2)
                nm = rpattern_shum(n)%idx_o(nn)
                if (nm .eq. 0) cycle
                rpattern_shum(n)%spec_o(nn,1,1) = rpattern_shum(n)%stdev*rpattern_shum(n)%spec_o(nn,1,1)*rpattern_shum(n)%varspectrum(nm)
                rpattern_shum(n)%spec_o(nn,2,1) = rpattern_shum(n)%stdev*rpattern_shum(n)%spec_o(nn,2,1)*rpattern_shum(n)%varspectrum(nm)
             enddo
          endif
       enddo
   endif

   if (nskeb > 0) then
  ! determine number of skeb levels to deal with temperoal/vertical correlations
      skeblevs=nint(skeb_tau(1)/skebint*skeb_vdof)
! backscatter noise.
      if (is_rootpe()) then
         print *, 'Initialize random pattern for SKEB'
         if (stochini) then
            ierr=NF90_INQ_VARID(stochlun,"skeb_seed", varid1)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring SKEB seed'
               iret = ierr
               return
            end if
            ierr=NF90_INQ_VARID(stochlun,"skeb_spec", varid2)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring SKEB spec'
               iret = ierr
               return
            end if
         endif
      endif
      call patterngenerator_init(skeb_lscale(1:nskeb),skebint,skeb_tau(1:nskeb),skeb(1:nskeb),iseed_skeb,rpattern_skeb, &
           lonf,latg,jcap,gis_stochy%ls_node,nskeb,skeblevs,skeb_varspect_opt,new_lscale)
      do n=1,nskeb
         do k=1,skeblevs
            if (stochini) then
               call read_pattern(rpattern_skeb(n),jcapin,stochlun,k,n,varid1,varid2,.true.,ierr)
               if (ierr .NE. 0) then
                  write(0,*) 'error reading SKEB pattern'
                  iret = ierr
                  return
               end if
            else
               call getnoise(rpattern_skeb(n),noise_e,noise_o)
               do nn=1,len_trie_ls
                  rpattern_skeb(n)%spec_e(nn,1,k)=noise_e(nn,1)
                  rpattern_skeb(n)%spec_e(nn,2,k)=noise_e(nn,2)
                  nm = rpattern_skeb(n)%idx_e(nn)
                  if (nm .eq. 0) cycle
                  rpattern_skeb(n)%spec_e(nn,1,k) = rpattern_skeb(n)%stdev*rpattern_skeb(n)%spec_e(nn,1,k)*rpattern_skeb(n)%varspectrum(nm)
                  rpattern_skeb(n)%spec_e(nn,2,k) = rpattern_skeb(n)%stdev*rpattern_skeb(n)%spec_e(nn,2,k)*rpattern_skeb(n)%varspectrum(nm)
               enddo
               do nn=1,len_trio_ls
                  rpattern_skeb(n)%spec_o(nn,1,k)=noise_o(nn,1)
                  rpattern_skeb(n)%spec_o(nn,2,k)=noise_o(nn,2)
                  nm = rpattern_skeb(n)%idx_o(nn)
                  if (nm .eq. 0) cycle
                  rpattern_skeb(n)%spec_o(nn,1,k) = rpattern_skeb(n)%stdev*rpattern_skeb(n)%spec_o(nn,1,k)*rpattern_skeb(n)%varspectrum(nm)
                  rpattern_skeb(n)%spec_o(nn,2,k) = rpattern_skeb(n)%stdev*rpattern_skeb(n)%spec_o(nn,2,k)*rpattern_skeb(n)%varspectrum(nm)
               enddo
            endif
         enddo
      enddo

      gis_stochy%kenorm_e=1.
      gis_stochy%kenorm_o=1. ! used to convert forcing pattern to wind field.
      if (skebnorm==0) then
       do locl=1,ls_max_node
           l = gis_stochy%ls_node(locl,1)
           jbasev = gis_stochy%ls_node(locl,2)
           indev = indlsev(l,l)
           jbasod = gis_stochy%ls_node(locl,3)
           indod = indlsod(l+1,l)
           do n=l,jcap,2
              rnn1 = n*(n+1.)
              gis_stochy%kenorm_e(indev) = rnn1/radius**2
              indev = indev + 1
           enddo
           do n=l+1,jcap,2
              rnn1 = n*(n+1.)
              gis_stochy%kenorm_o(indod) = rnn1/radius**2
              indod = indod + 1
           enddo
        enddo
        if (is_rootpe()) print*,'using streamfunction ',maxval(gis_stochy%kenorm_e(:)),minval(gis_stochy%kenorm_e(:))
      endif
      if (skebnorm==1) then
       do locl=1,ls_max_node
           l = gis_stochy%ls_node(locl,1)
           jbasev = gis_stochy%ls_node(locl,2)
           indev = indlsev(l,l)
           jbasod = gis_stochy%ls_node(locl,3)
           indod = indlsod(l+1,l)
           do n=l,jcap,2
              rnn1 = n*(n+1.)
              gis_stochy%kenorm_e(indev) = sqrt(rnn1)/radius
              indev = indev + 1
           enddo
           do n=l+1,jcap,2
              rnn1 = n*(n+1.)
              gis_stochy%kenorm_o(indod) = sqrt(rnn1)/radius
              indod = indod + 1
           enddo
        enddo
        if (is_rootpe()) print*,'using kenorm ',maxval(gis_stochy%kenorm_e(:)),minval(gis_stochy%kenorm_e(:))
      endif
      ! set the even and odd (n-l) terms of the top row to zero
      do locl=1,ls_max_node
         l = gis_stochy%ls_node(locl,1)
         jbasev = gis_stochy%ls_node(locl,2)
         jbasod = gis_stochy%ls_node(locl,3)
         if (mod(l,2) .eq. mod(jcap+1,2)) then
            gis_stochy%kenorm_e(indlsev(jcap+1,l)) = 0.
         endif
         if (mod(l,2) .ne. mod(jcap+1,2)) then
            gis_stochy%kenorm_o(indlsod(jcap+1,l)) = 0.
         endif
      enddo
      
   endif ! skeb > 0
! mg, sfc-perts
   if (nlndp > 0) then
      if (is_rootpe()) then
         print *, 'Initialize random pattern for SFC-PERTS'
         if (stochini) then
            ierr=NF90_INQ_VARID(stochlun,"sfcpert_seed", varid1)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring SFC-PERTS seed'
               iret = ierr
               return
            end if
            ierr=NF90_INQ_VARID(stochlun,"sfcpert_spec", varid2)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring SFC-PERTS spec'
               iret = ierr
               return
            end if
         endif
      endif
      ones = 1.
      call patterngenerator_init(lndp_lscale(1:nlndp),delt,lndp_tau(1:nlndp),ones(1:nlndp),iseed_lndp,rpattern_sfc, &
                                 lonf,latg,jcap,gis_stochy%ls_node,nlndp,n_var_lndp,0,new_lscale)
      do n=1,nlndp
         if (is_rootpe()) print *, 'Initialize random pattern for LNDP PERTS'
         do k=1,n_var_lndp
            if (stochini) then
               call read_pattern(rpattern_sfc(n),jcapin,stochlun,k,n,varid1,varid2,.true.,ierr)
               if (ierr .NE. 0) then
                  write(0,*) 'error reading SHUM pattern'
                  iret = ierr
                  return
               endif
               if (is_rootpe()) print *, 'lndp pattern read',n,k,minval(rpattern_sfc(n)%spec_o(:,:,k)), maxval(rpattern_sfc(n)%spec_o(:,:,k))
            else
                call getnoise(rpattern_sfc(n),noise_e,noise_o)
                do nn=1,len_trie_ls
                   rpattern_sfc(n)%spec_e(nn,1,k)=noise_e(nn,1)
                   rpattern_sfc(n)%spec_e(nn,2,k)=noise_e(nn,2)
                   nm = rpattern_sfc(n)%idx_e(nn)
                   if (nm .eq. 0) cycle
                   rpattern_sfc(n)%spec_e(nn,1,k) = rpattern_sfc(n)%stdev*rpattern_sfc(n)%spec_e(nn,1,k)*rpattern_sfc(n)%varspectrum(nm)
                   rpattern_sfc(n)%spec_e(nn,2,k) = rpattern_sfc(n)%stdev*rpattern_sfc(n)%spec_e(nn,2,k)*rpattern_sfc(n)%varspectrum(nm)
                enddo
                do nn=1,len_trio_ls
                   rpattern_sfc(n)%spec_o(nn,1,k)=noise_o(nn,1)
                   rpattern_sfc(n)%spec_o(nn,2,k)=noise_o(nn,2)
                   nm = rpattern_sfc(n)%idx_o(nn)
                   if (nm .eq. 0) cycle
                   rpattern_sfc(n)%spec_o(nn,1,k) = rpattern_sfc(n)%stdev*rpattern_sfc(n)%spec_o(nn,1,k)*rpattern_sfc(n)%varspectrum(nm)
                   rpattern_sfc(n)%spec_o(nn,2,k) = rpattern_sfc(n)%stdev*rpattern_sfc(n)%spec_o(nn,2,k)*rpattern_sfc(n)%varspectrum(nm)
                enddo
                if (is_rootpe()) print *, 'lndp pattern initialized, ',n, k, minval(rpattern_sfc(n)%spec_o(:,:,k)), maxval(rpattern_sfc(n)%spec_o(:,:,k))
            endif ! stochini
         enddo ! k, n_var_lndp
      enddo ! n, nlndp
   endif ! nlndp > 0
   if (nspp  > 0) then
      if (is_rootpe()) then
         print *, 'Initialize random pattern for SPP-PERTS'
         if (stochini) then
            ierr=NF90_INQ_VARID(stochlun,"spppert_seed", varid1)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring SPP-PERTS seed'
               iret = ierr
               return
            end if
            ierr=NF90_INQ_VARID(stochlun,"ppcpert_spec", varid2)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring SPP-PERTS spec'
               iret = ierr
               return
            end if
         endif
      endif
      ones = 1.
      call patterngenerator_init(spp_lscale(1:nspp),delt,spp_tau(1:nspp),ones(1:nspp),iseed_spp,rpattern_spp, &
                                 lonf,latg,jcap,gis_stochy%ls_node,nspp,n_var_spp,0,new_lscale)
      do n=1,nspp
         if (is_rootpe()) print *, 'Initialize random pattern for SPP PERTS'
            if (stochini) then
               call read_pattern(rpattern_spp(n),jcapin,stochlun,1,n,varid1,varid2,.true.,ierr)
               if (ierr .NE. 0) then
                  write(0,*) 'error reading SPP  pattern'
                  iret = ierr
                  return
               endif
               if (is_rootpe()) print *, 'spp  pattern read',n,1,minval(rpattern_spp(n)%spec_o(:,:,1)), maxval(rpattern_spp(n)%spec_o(:,:,1))
            else
                call getnoise(rpattern_spp(n),noise_e,noise_o)
                do nn=1,len_trie_ls
                   rpattern_spp(n)%spec_e(nn,1,1)=noise_e(nn,1)
                   rpattern_spp(n)%spec_e(nn,2,1)=noise_e(nn,2)
                   nm = rpattern_spp(n)%idx_e(nn)
                   if (nm .eq. 0) cycle
                   rpattern_spp(n)%spec_e(nn,1,1) = rpattern_spp(n)%stdev*rpattern_spp(n)%spec_e(nn,1,1)*rpattern_spp(n)%varspectrum(nm)
                   rpattern_spp(n)%spec_e(nn,2,1) = rpattern_spp(n)%stdev*rpattern_spp(n)%spec_e(nn,2,1)*rpattern_spp(n)%varspectrum(nm)
                enddo
                do nn=1,len_trio_ls
                   rpattern_spp(n)%spec_o(nn,1,1)=noise_o(nn,1)
                   rpattern_spp(n)%spec_o(nn,2,1)=noise_o(nn,2)
                   nm = rpattern_spp(n)%idx_o(nn)
                   if (nm .eq. 0) cycle
                   rpattern_spp(n)%spec_o(nn,1,1) = rpattern_spp(n)%stdev*rpattern_spp(n)%spec_o(nn,1,1)*rpattern_spp(n)%varspectrum(nm)
                   rpattern_spp(n)%spec_o(nn,2,1) = rpattern_spp(n)%stdev*rpattern_spp(n)%spec_o(nn,2,1)*rpattern_spp(n)%varspectrum(nm)
                enddo
                if (is_rootpe()) print *, 'spp pattern initialized, ',n, 1, minval(rpattern_spp(n)%spec_o(:,:,1)), maxval(rpattern_spp(n)%spec_o(:,:,1))
            endif ! stochini
      enddo ! n, nspp
   endif ! nspp  > 0

   if (is_rootpe() .and. stochini) CLOSE(stochlun)
   deallocate(noise_e,noise_o)
 end subroutine init_stochdata

 subroutine init_stochdata_ocn(nlevs,delt,iret)

 use netcdf
 use compns_stochy_mod, only : compns_stochy_ocn
 use mpp_domains_mod,     only: mpp_broadcast_domain,MPP_DOMAIN_TIME,mpp_domains_init ,mpp_domains_set_stack_size
! initialize random patterns.  A spinup period of spinup_efolds times the
! temporal time scale is run for each pattern.
   integer, intent(in) :: nlevs
   real, intent(in) :: delt
   integer, intent(out) :: iret
   
   integer :: nn,nm,stochlun,n,jcapin
   integer :: l,jbasev,jbasod
   integer :: indev,indod,indlsod,indlsev,varid1,varid2,varid3,varid4,ierr
   
   real(kind_dbl_prec),allocatable :: noise_e(:,:),noise_o(:,:)
   include 'function_indlsod'
   include 'function_indlsev'
   include 'netcdf.inc'
   stochlun=99
   levs=nlevs

   iret=0
   call compns_stochy_ocn (delt,iret)
   if(is_rootpe()) print*,'in init stochdata_ocn'
   if ( (.NOT. pert_epbl) .AND. (.NOT. do_ocnsppt) ) return
   call initialize_spectral(gis_stochy_ocn)
   if (iret/=0) return
   allocate(noise_e(len_trie_ls,2),noise_o(len_trio_ls,2))
! determine number of random patterns to be used for each scheme.
   do n=1,size(epbl)
     if (epbl(n) > 0) then
        nepbl=nepbl+1
     else
        exit
     endif
   enddo

   do n=1,size(ocnsppt)
     if (ocnsppt(n) > 0) then
        nocnsppt=nocnsppt+1
     else
        exit
     endif
   enddo

   if (nepbl > 0) then 
      allocate(rpattern_epbl1(nepbl))
      allocate(rpattern_epbl2(nepbl))
   endif

   if (nocnsppt > 0) allocate(rpattern_ocnsppt(nocnsppt))

!  if stochini is true, then read in pattern from a file
   if (is_rootpe()) then
      if (stochini) then
         print*,'opening stoch_ini'
         ierr=nf90_open('INPUT/ocn_stoch.res.nc',nf90_nowrite,ncid=stochlun)
         if (ierr .NE. 0) then
            write(0,*) 'error opening stoch_ini'
            iret = ierr
            return
         end if
         ierr=NF90_GET_ATT(stochlun,NF_GLOBAL,"ntrunc",jcapin)
         if (ierr .NE. 0) then
            write(0,*) 'error getting ntrunc'
            iret = ierr
            return
         end if
         print*,'ntrunc read in',jcapin
      endif
   endif

   if (nepbl > 0) then
      if (is_rootpe()) then
         print *, 'Initialize random pattern for epbl'
         if (stochini) then
            ierr=NF90_INQ_VARID(stochlun,"epbl1_seed", varid1)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring EPBL1 seed'
               iret = ierr
               return
            end if
            ierr=NF90_INQ_VARID(stochlun,"epbl1_spec", varid2)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring EPBL1 spec'
               iret = ierr
               return
            end if
            ierr=NF90_INQ_VARID(stochlun,"epbl2_seed", varid3)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring EPBL2 seed'
               iret = ierr
               return
            end if
            ierr=NF90_INQ_VARID(stochlun,"epbl2_spec", varid4)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring EPBL2 spec'
               iret = ierr
               return
            end if
         end if
      end if
      call patterngenerator_init(epbl_lscale(1:nepbl),epblint,epbl_tau(1:nepbl),epbl(1:nepbl),iseed_epbl,rpattern_epbl1, &
           lonf,latg,jcap,gis_stochy_ocn%ls_node,nepbl,1,0,new_lscale)
      call patterngenerator_init(epbl_lscale(1:nepbl),epblint,epbl_tau(1:nepbl),epbl(1:nepbl),iseed_epbl2,rpattern_epbl2, &
           lonf,latg,jcap,gis_stochy_ocn%ls_node,nepbl,1,0,new_lscale)
      do n=1,nepbl
         if (stochini) then
            call read_pattern(rpattern_epbl1(n),jcapin,stochlun,1,n,varid1,varid2,.false.,ierr)
            if (ierr .NE. 0) then
               write(0,*) 'error reading EPBL1 pattern'
               iret = ierr
               return
            end if
            call read_pattern(rpattern_epbl2(n),jcapin,stochlun,1,n,varid3,varid4,.false.,ierr)
            if (ierr .NE. 0) then
               write(0,*) 'error reading EPBL1 pattern'
               iret = ierr
               return
            end if
         else
            call getnoise(rpattern_epbl1(n),noise_e,noise_o)
            do nn=1,len_trie_ls
               rpattern_epbl1(n)%spec_e(nn,1,1)=noise_e(nn,1)
               rpattern_epbl1(n)%spec_e(nn,2,1)=noise_e(nn,2)
               rpattern_epbl1(n)%spec_e(nn,1,1)=noise_e(nn,1)
               rpattern_epbl1(n)%spec_e(nn,2,1)=noise_e(nn,2)
               nm = rpattern_epbl1(n)%idx_e(nn)
               if (nm .eq. 0) cycle
               rpattern_epbl1(n)%spec_e(nn,1,1) = rpattern_epbl1(n)%stdev*rpattern_epbl1(n)%spec_e(nn,1,1)*rpattern_epbl1(n)%varspectrum(nm)
               rpattern_epbl1(n)%spec_e(nn,2,1) = rpattern_epbl1(n)%stdev*rpattern_epbl1(n)%spec_e(nn,2,1)*rpattern_epbl1(n)%varspectrum(nm)
            enddo
            do nn=1,len_trio_ls
               rpattern_epbl1(n)%spec_o(nn,1,1)=noise_o(nn,1)
               rpattern_epbl1(n)%spec_o(nn,2,1)=noise_o(nn,2)
               nm = rpattern_epbl1(n)%idx_o(nn)
               if (nm .eq. 0) cycle
               rpattern_epbl1(n)%spec_o(nn,1,1) = rpattern_epbl1(n)%stdev*rpattern_epbl1(n)%spec_o(nn,1,1)*rpattern_epbl1(n)%varspectrum(nm)
               rpattern_epbl1(n)%spec_o(nn,2,1) = rpattern_epbl1(n)%stdev*rpattern_epbl1(n)%spec_o(nn,2,1)*rpattern_epbl1(n)%varspectrum(nm)
            enddo
            call patterngenerator_advance(rpattern_epbl1(n),1,.false.)
         
            call getnoise(rpattern_epbl2(n),noise_e,noise_o)
            do nn=1,len_trie_ls
               rpattern_epbl2(n)%spec_e(nn,1,1)=noise_e(nn,1)
               rpattern_epbl2(n)%spec_e(nn,2,1)=noise_e(nn,2)
               rpattern_epbl2(n)%spec_e(nn,1,1)=noise_e(nn,1)
               rpattern_epbl2(n)%spec_e(nn,2,1)=noise_e(nn,2)
               nm = rpattern_epbl2(n)%idx_e(nn)
               if (nm .eq. 0) cycle
               rpattern_epbl2(n)%spec_e(nn,1,1) = rpattern_epbl2(n)%stdev*rpattern_epbl2(n)%spec_e(nn,1,1)*rpattern_epbl2(n)%varspectrum(nm)
               rpattern_epbl2(n)%spec_e(nn,2,1) = rpattern_epbl2(n)%stdev*rpattern_epbl2(n)%spec_e(nn,2,1)*rpattern_epbl2(n)%varspectrum(nm)
            enddo
            do nn=1,len_trio_ls
               rpattern_epbl2(n)%spec_o(nn,1,1)=noise_o(nn,1)
               rpattern_epbl2(n)%spec_o(nn,2,1)=noise_o(nn,2)
               nm = rpattern_epbl2(n)%idx_o(nn)
               if (nm .eq. 0) cycle
               rpattern_epbl2(n)%spec_o(nn,1,1) = rpattern_epbl2(n)%stdev*rpattern_epbl2(n)%spec_o(nn,1,1)*rpattern_epbl2(n)%varspectrum(nm)
               rpattern_epbl2(n)%spec_o(nn,2,1) = rpattern_epbl2(n)%stdev*rpattern_epbl2(n)%spec_o(nn,2,1)*rpattern_epbl2(n)%varspectrum(nm)
            enddo
            call patterngenerator_advance(rpattern_epbl2(n),1,.false.)
         endif
      enddo
   endif

   if (nocnsppt > 0) then
      if (is_rootpe()) then
         if (stochini) then
            ierr=NF90_INQ_VARID(stochlun,"ocnsppt_seed", varid1)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring OCNSPPT seed'
               iret = ierr
               return
            end if
            ierr=NF90_INQ_VARID(stochlun,"ocnsppt_spec", varid2)
            if (ierr .NE. 0) then
               write(0,*) 'error inquring OCNSPPT spec'
               iret = ierr
               return
            end if
         endif
      endif
      if (is_rootpe()) print *, 'Initialize random pattern for ocnsppt'
      call patterngenerator_init(ocnsppt_lscale(1:nocnsppt),ocnspptint,ocnsppt_tau(1:nocnsppt),ocnsppt(1:nocnsppt),iseed_ocnsppt,rpattern_ocnsppt, &
           lonf,latg,jcap,gis_stochy_ocn%ls_node,nocnsppt,1,0,new_lscale)
      do n=1,nocnsppt
         if (stochini) then
            call read_pattern(rpattern_ocnsppt(n),jcapin,stochlun,1,n,varid1,varid2,.false.,ierr)
            if (ierr .NE. 0) then
               write(0,*) 'error reading SPPT pattern'
               iret = ierr
               return
            end if
         else
            call getnoise(rpattern_ocnsppt(n),noise_e,noise_o)
            do nn=1,len_trie_ls
               rpattern_ocnsppt(n)%spec_e(nn,1,1)=noise_e(nn,1)
               rpattern_ocnsppt(n)%spec_e(nn,2,1)=noise_e(nn,2)
               nm = rpattern_ocnsppt(n)%idx_e(nn)
               if (nm .eq. 0) cycle
               rpattern_ocnsppt(n)%spec_e(nn,1,1) = rpattern_ocnsppt(n)%stdev*rpattern_ocnsppt(n)%spec_e(nn,1,1)*rpattern_ocnsppt(n)%varspectrum(nm)
               rpattern_ocnsppt(n)%spec_e(nn,2,1) = rpattern_ocnsppt(n)%stdev*rpattern_ocnsppt(n)%spec_e(nn,2,1)*rpattern_ocnsppt(n)%varspectrum(nm)
            enddo
            do nn=1,len_trio_ls
               rpattern_ocnsppt(n)%spec_o(nn,1,1)=noise_o(nn,1)
               rpattern_ocnsppt(n)%spec_o(nn,2,1)=noise_o(nn,2)
               nm = rpattern_ocnsppt(n)%idx_o(nn)
            if (nm .eq. 0) cycle
               rpattern_ocnsppt(n)%spec_o(nn,1,1) = rpattern_ocnsppt(n)%stdev*rpattern_ocnsppt(n)%spec_o(nn,1,1)*rpattern_ocnsppt(n)%varspectrum(nm)
              rpattern_ocnsppt(n)%spec_o(nn,2,1) = rpattern_ocnsppt(n)%stdev*rpattern_ocnsppt(n)%spec_o(nn,2,1)*rpattern_ocnsppt(n)%varspectrum(nm)
            enddo
            call patterngenerator_advance(rpattern_ocnsppt(n),1,.false.)
         endif
      enddo
   endif
   deallocate(noise_e,noise_o)
 end subroutine init_stochdata_ocn


!>@brief This subroutine 'read_pattern' will read in the spectral coeffients from a previous run (stored in stoch_ini,
!!turned on by setting STOCHINI=.true.)
!>@details Data read in are flat binary, so the number of stochastic physics patterns running must match previous run
subroutine read_pattern(rpattern,jcapin,lunptn,k,np,varid1,varid2,slice_of_3d,iret)
!\callgraph
   use netcdf
   type(random_pattern), intent(inout) :: rpattern
   integer, intent(in) :: lunptn,np,varid1,varid2,jcapin
   logical, intent(in) :: slice_of_3d
   real(kind_dbl_prec),allocatable  :: pattern2d(:),pattern2din(:)
   real(kind_dbl_prec) :: stdevin,varin
   integer nm,nn,iret,ierr,isize,k,ndimspec2
   integer, allocatable :: isave(:)
   include 'netcdf.inc'
   iret=0
   ndimspec2=2*ndimspec
   allocate(pattern2d(ndimspec2))
   pattern2d=0.
   call random_seed(size=isize,stat=rpattern%rstate)  ! get size of generator state seed array
   allocate(isave(isize))
   ! read only on root process, and send to all tasks
   if (is_rootpe()) then
      allocate(pattern2din((jcapin+1)*(jcapin+2)))
      print*,'reading in random pattern at ',jcapin,ndimspec,size(pattern2din)
      !read(lunptn) pattern2din
      ierr=NF90_GET_VAR(lunptn,varid1,isave,(/1,np/))
      if (ierr .NE. 0) then
         write(0,*) 'error reading seed'
         iret = ierr
         return
      end if
      if (slice_of_3d) then
         ierr=NF90_GET_VAR(lunptn,varid2,pattern2din,(/1,k,np/),(/ndimspec2,1,1/))
      else
         ierr=NF90_GET_VAR(lunptn,varid2,pattern2din,(/1,np/),(/ndimspec2,1/))
      endif
      if (ierr .NE. 0) then
         write(0,*) 'error reading spec var'
         iret = ierr
         return
      end if
      print*,'reading in random pattern (min/max/size/seed)',&
      minval(pattern2din),maxval(pattern2din),size(pattern2din),isave(1:4)
      if (jcapin .eq. ntrunc) then
         pattern2d=pattern2din
      else
         call chgres_pattern(pattern2din,pattern2d,jcapin,ntrunc) ! chgres of spectral files
         ! change the standard deviation of the patterns for a resolution change
         ! needed for SKEB & SHUM
         call computevarspec_r(rpattern,pattern2d,varin)
         print*,'stddev in and out..',sqrt(varin),rpattern%stdev
         stdevin=rpattern%stdev/sqrt(varin)
         pattern2d(:)=pattern2d(:)*stdevin
      endif
      deallocate(pattern2din)
    endif
    call mp_bcst(isave,isize)  ! blast out seed
    call mp_bcst(pattern2d,2*ndimspec)
    call random_seed(put=isave,stat=rpattern%rstate)
   ! subset
   do nn=1,len_trie_ls
      nm = rpattern%idx_e(nn)
      if (nm == 0) cycle
      rpattern%spec_e(nn,1,k) = pattern2d(nm)
      rpattern%spec_e(nn,2,k) = pattern2d(ndimspec+nm)
   enddo
   do nn=1,len_trio_ls
      nm = rpattern%idx_o(nn)
      if (nm == 0) cycle
      rpattern%spec_o(nn,1,k) = pattern2d(nm)
      rpattern%spec_o(nn,2,k) = pattern2d(ndimspec+nm)
   enddo
   !print*,'after scatter...',mype,maxval(pattern2d_e),maxval(pattern2d_o) &
   ! ,minval(pattern2d_e),minval(pattern2d_o)
   deallocate(pattern2d,isave)
 end subroutine read_pattern

end module stochy_data_mod
