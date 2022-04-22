module lndp_apply_perts_mod

    use kinddef, only : kind_dbl_prec

    implicit none

    private

    public :: lndp_apply_perts

    contains

!====================================================================
! lndp_apply_perts
!====================================================================
! Driver for applying perturbations to sprecified land states or parameters
! Draper, July 2020.

    subroutine lndp_apply_perts(blksz, lsm, lsm_noah, lsm_ruc, lsoil,           &
                dtf, kdt, lndp_each_step,                                       &
                n_var_lndp, lndp_var_list, lndp_prt_list,                       &
                sfc_wts, xlon, xlat, stype, smcmax, smcmin, param_update_flag,  &
                smc, slc, stc, vfrac, alvsf, alnsf, alvwf, alnwf, facsf, facwf, &
                snoalb, semis, zorll, ierr)

        implicit none

        ! intent(in)
        integer,                      intent(in) :: blksz(:)
        integer,                      intent(in) :: n_var_lndp, lsoil, kdt
        logical,                      intent(in) :: lndp_each_step
        integer,                      intent(in) :: lsm, lsm_noah, lsm_ruc
        character(len=3),             intent(in) :: lndp_var_list(:)
        real(kind=kind_dbl_prec),     intent(in) :: lndp_prt_list(:)
        real(kind=kind_dbl_prec),     intent(in) :: dtf
        real(kind=kind_dbl_prec),     intent(in) :: sfc_wts(:,:,:)
        real(kind=kind_dbl_prec),     intent(in) :: xlon(:,:)
        real(kind=kind_dbl_prec),     intent(in) :: xlat(:,:)
        logical,                      intent(in) :: param_update_flag
                                        ! true =  parameters have been updated, apply perts
        real(kind=kind_dbl_prec),     intent(in) :: stype(:,:)
        real(kind=kind_dbl_prec),     intent(in) :: smcmax(:)
        real(kind=kind_dbl_prec),     intent(in) :: smcmin(:)

        ! intent(inout)
        real(kind=kind_dbl_prec),     intent(inout) :: smc(:,:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: slc(:,:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: stc(:,:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: vfrac(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: snoalb(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: alvsf(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: alnsf(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: alvwf(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: alnwf(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: facsf(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: facwf(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: semis(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: zorll(:,:)

        ! intent(out)
        integer,                        intent(out) :: ierr

        ! local
        integer         :: nblks, print_i, print_nb, i, nb
        integer         :: this_im, v, soiltyp, k
        logical         :: print_flag

        real(kind=kind_dbl_prec) :: p, min_bound, max_bound, tmp_sic,  pert, factor
        real(kind=kind_dbl_prec), dimension(lsoil) :: zslayer, smc_vertscale, stc_vertscale

        ! decrease in applied pert with depth
        !-- Noah lsm
        real(kind=kind_dbl_prec), dimension(4), parameter  :: smc_vertscale_noah = (/1.0,0.5,0.25,0.125/)
        real(kind=kind_dbl_prec), dimension(4), parameter  :: stc_vertscale_noah = (/1.0,0.5,0.25,0.125/)
        real(kind=kind_dbl_prec), dimension(4), parameter  :: zs_noah = (/0.1, 0.3, 0.6, 1.0/)
        !-- RUC lsm
        real(kind=kind_dbl_prec), dimension(9), parameter :: smc_vertscale_ruc = (/1.0,0.9,0.8,0.6,0.4,0.2,0.1,0.05,0./)
        real(kind=kind_dbl_prec), dimension(9), parameter :: stc_vertscale_ruc = (/1.0,0.9,0.8,0.6,0.4,0.2,0.1,0.05,0./)
        real(kind=kind_dbl_prec), dimension(9), parameter :: zs_ruc = (/0.05, 0.15, 0.20, 0.20, 0.40, 0.60, 0.60, 0.80, 1.00/)

        ierr = 0

        if (lsm/=lsm_noah .and. lsm/=lsm_ruc) then
          write(6,*) 'ERROR: lndp_apply_pert assumes LSM is Noah or RUC,', &
                     ' may need to adapt variable names for a different LSM'
          ierr=10
          return
        endif

        !write (0,*) 'Input to lndp_apply_pert'
        !write (0,*) 'lsm, lsoil, lsm_ruc, lsoil_lsm =', lsm, lsoil, lsm_ruc, lsoil_lsm
        !write (0,*) 'zs_lsm =', zs_lsm
        !write (0,*) 'n_var_lndp, lndp_var_list =', n_var_lndp, lndp_var_list
        !write (0,*) 'smcmin =',smcmin

        ! lndp_prt_list input is per hour, factor converts to per timestep
        ! Do conversion only when variables are perturbed at every time step
        if(lndp_each_step) then
          factor = dtf/3600.
        else
          factor = 1.
        endif

        if (lsm == lsm_noah) then
          do k = 1, lsoil
            zslayer(k) = zs_noah(k)
            smc_vertscale(k) = smc_vertscale_noah(k)
            stc_vertscale(k) = stc_vertscale_noah(k)
          enddo
        elseif (lsm == lsm_ruc) then
          do k = 1, lsoil
            zslayer(k) = zs_ruc(k)
            smc_vertscale(k) = smc_vertscale_ruc(k)
            stc_vertscale(k) = stc_vertscale_ruc(k)
          enddo
        endif

        nblks = size(blksz)

        call  set_printing_nb_i(blksz,xlon,xlat,print_i,print_nb)

        do nb =1,nblks
           do i = 1, blksz(nb)

             !if ( nint(Sfcprop(nb)%slmsk(i)) .NE. 1) cycle ! skip if not land

             !if ( ((isot == 1) .and. (soiltyp == 16)) &
             !  .or.( (isot == 0) .and. (soiltyp  == 9 )) ) cycle ! skip if land-ice

             if ( smc(nb,i,1) .EQ. 1.) cycle ! skip  non-soil (land-ice and non-land)
             ! set printing
             if ( (i==print_i)  .and. (nb==print_nb) ) then
                print_flag = .true.
             else
                print_flag=.false.
             endif

             do v = 1,n_var_lndp
                select case (trim(lndp_var_list(v)))
                !=================================================================
                ! State updates - performed every cycle
                !=================================================================
                case('smc')
                    p=5.
                    soiltyp  = int( stype(nb,i)+0.5 )  ! also need for maxsmc
                    min_bound = smcmin(soiltyp)
                    max_bound = smcmax(soiltyp)

                  if ((lsm /= lsm_ruc) .or. (lsm == lsm_ruc .and. kdt == 2)) then
                  ! with RUC LSM perturb smc only at time step = 2, as in HRRR
                    do k=1,lsoil
                         !store frozen soil moisture
                         tmp_sic= smc(nb,i,k)  - slc(nb,i,k)

                         ! perturb total soil moisture
                         ! factor of sldepth*1000 converts from mm to m3/m3
                         ! lndp_prt_list(v) = 0.3 in input.nml
                         pert = sfc_wts(nb,i,v)*smc_vertscale(k)*lndp_prt_list(v)/(zslayer(k)*1000.)
                         pert = pert*dtf/3600. ! lndp_prt_list input is per hour, convert to per timestep
                                                     ! (necessary for state vars only)
                         call apply_pert('smc',pert,print_flag, smc(nb,i,k),ierr,p,min_bound, max_bound)

                         ! assign all of applied pert to the liquid soil moisture
                         slc(nb,i,k)  =  smc(nb,i,k) -  tmp_sic
                    enddo
                  endif

                case('stc')

                    do k=1,lsoil
                         pert = sfc_wts(nb,i,v)*stc_vertscale(k)*lndp_prt_list(v)
                         pert = pert*dtf/3600. ! lndp_prt_list input is per hour, convert to per timestep
                                                     ! (necessary for state vars only)
                         call apply_pert('stc',pert,print_flag, stc(nb,i,k),ierr)
                    enddo
                !=================================================================
                ! Parameter updates - only if param_update_flag = TRUE
                !=================================================================
                case('vgf')  ! vegetation fraction
                     if (param_update_flag .or. lndp_each_step) then
                         p =5.
                         min_bound=0.
                         max_bound=1.

                         pert = sfc_wts(nb,i,v)*lndp_prt_list(v)
                         pert = pert*factor 
                         call apply_pert ('vfrac',pert,print_flag, vfrac(nb,i), ierr,p,min_bound, max_bound)
                     endif
                case('alb')  ! albedo
                     if (param_update_flag .or. lndp_each_step) then
                         p =5.
                         min_bound=0.0
                         max_bound=0.4

                         pert = sfc_wts(nb,i,v)*lndp_prt_list(v)
                         pert = pert*factor
                         !call apply_pert ('alvsf',pert,print_flag, alvsf(nb,i), ierr,p,min_bound, max_bound)
                         call apply_pert ('alnsf',pert,print_flag, alnsf(nb,i), ierr,p,min_bound, max_bound)
                         !call apply_pert ('alvwf',pert,print_flag, alvwf(nb,i), ierr,p,min_bound, max_bound)
                         call apply_pert ('alnwf',pert,print_flag, alnwf(nb,i), ierr,p,min_bound, max_bound)
                         !call apply_pert ('facsf',pert,print_flag, facsf(nb,i), ierr,p,min_bound, max_bound)
                         !call apply_pert ('facwf',pert,print_flag, facwf(nb,i), ierr,p,min_bound, max_bound)
                     endif
                case('sal')  ! snow albedo
                     if (param_update_flag .or. lndp_each_step) then
                         p =5.
                         min_bound=0.3
                         max_bound=0.85

                         pert = sfc_wts(nb,i,v)*lndp_prt_list(v)
                         pert = pert*factor
                         call apply_pert ('snoalb',pert,print_flag, snoalb(nb,i), ierr,p,min_bound, max_bound)
                     endif
                case('emi')  ! emissivity
                     if (param_update_flag .or. lndp_each_step) then
                         p =5.
                         min_bound=0.8
                         max_bound=1.

                         pert = sfc_wts(nb,i,v)*lndp_prt_list(v)
                         pert = pert*factor
                         call apply_pert ('semis',pert,print_flag, semis(nb,i), ierr,p,min_bound, max_bound)
                     endif
                case('zol')  ! land roughness length
                     if (param_update_flag .or. lndp_each_step) then
                         p =5.
                         min_bound=0.
                         max_bound=300.

                         pert = sfc_wts(nb,i,v)*lndp_prt_list(v)
                         pert = pert*factor
                         call apply_pert ('zol',pert,print_flag, zorll(nb,i), ierr,p,min_bound, max_bound)
                     endif
                case default
                    print*, &
                     'ERROR: unrecognised lndp_prt_list option in lndp_apply_pert, exiting', trim(lndp_var_list(v))
                    ierr = 10
                    return
                end select
             enddo
           enddo
        enddo
    end subroutine  lndp_apply_perts

!====================================================================
! apply_perts
!====================================================================
! Apply perturbations to selected land states or parameters

  subroutine apply_pert(vname,pert,print_flag, state,ierr,p,vmin, vmax)

   ! intent in
    logical, intent(in)                 :: print_flag
    real(kind=kind_dbl_prec), intent(in)    :: pert
    character(len=*), intent(in)        :: vname ! name of variable being perturbed

    real(kind=kind_dbl_prec), optional, intent(in)    :: p ! flat-top paramater, 0 = no flat-top
                                                       ! flat-top function is used for bounded variables
                                                       ! to reduce the magnitude of perturbations  near boundaries.
    real(kind=kind_dbl_prec), optional, intent(in) :: vmin, vmax ! min,max bounds of variable being perturbed

    ! intent (inout)
    real(kind=kind_dbl_prec), intent(inout) :: state

    ! intent out
    integer                             :: ierr

    !local
    real(kind=kind_dbl_prec) :: z

       if ( print_flag ) then
              write(*,*) 'LNDP - applying lndp to ',vname, ', initial value', state
       endif

       ! apply perturbation
       if (present(p) ) then
           if ( .not. (present(vmin) .and. present(vmax) )) then
              ierr=20
              print*, 'error, flat-top function requires min & max to be specified'
           endif

           z = -1. + 2*(state  - vmin)/(vmax - vmin) ! flat-top function
           state =  state  + pert*(1-abs(z**p))
       else
          state =  state  + pert
       endif

       if (present(vmax)) state =  min( state , vmax )
       if (present(vmin)) state =  max( state , vmin )
       !state = max( min( state , vmax ), vmin )

       if ( print_flag ) then
              write(*,*) 'LNDP - applying lndp to ',vname, ', final value', state
       endif

  end subroutine apply_pert


!====================================================================
! set_printing_nb_i
!====================================================================
! routine to turn on print flag for selected location
!
    subroutine set_printing_nb_i(blksz,xlon,xlat,print_i,print_nb)

        implicit none

        ! intent (in)
        integer,                  intent(in) :: blksz(:)
        real(kind=kind_dbl_prec),     intent(in) :: xlon(:,:)
        real(kind=kind_dbl_prec),     intent(in) :: xlat(:,:)


        ! intent (out)
        integer, intent(out) :: print_i, print_nb

        ! local
        integer :: nblks,nb,i
        real, parameter :: plon_trunc =  114.9
        real, parameter :: plat_trunc =  -26.6
        real, parameter  :: delta  = 1.

        nblks = size(blksz)

        print_i = -9
        print_nb = -9
        do nb = 1,nblks
         do i = 1,blksz(nb)
        if ( (xlon(nb,i)*57.29578 > plon_trunc) .and.  (xlon(nb,i)*57.29578 < plon_trunc+delta ) .and. &
           (xlat(nb,i)*57.29578 >  plat_trunc ) .and.  (xlat(nb,i)*57.29578 < plat_trunc+delta ) ) then
                      print_i=i
                      print_nb=nb
                      write(*,*) 'LNDP -print flag is on', xlon(nb,i)*57.29578, xlat(nb,i)*57.29578, nb, i
                      return
         endif
         enddo
        enddo

    end subroutine set_printing_nb_i

end module lndp_apply_perts_mod


