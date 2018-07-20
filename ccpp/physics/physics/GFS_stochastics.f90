!> \file GFS_stochastics.f90
!!  Contains code previously in GFS_stochastics_driver.

    module GFS_stochastics

      contains

      subroutine GFS_stochastics_init ()
      end subroutine GFS_stochastics_init

      subroutine GFS_stochastics_finalize()
      end subroutine GFS_stochastics_finalize

!> \section arg_table_GFS_stochastics_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type              |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | derived type GFS_control_type in FV3                    | DDT           |    0 | GFS_control_type  |           | in     | F        |
!! | Statein        | FV3-GFS_Statein_type                                   | derived type GFS_statein_type in FV3                    | DDT           |    0 | GFS_statein_type  |           | in     | F        |
!! | Stateout       | FV3-GFS_Stateout_type                                  | derived type GFS_stateout_type in FV3                   | DDT           |    0 | GFS_stateout_type |           | in     | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                                   | derived type GFS_sfcprop_type in FV3                    | DDT           |    0 | GFS_sfcprop_type  |           | in     | F        |
!! | Coupling       | FV3-GFS_Coupling_type                                  | derived type GFS_coupling_type in FV3                   | DDT           |    0 | GFS_coupling_type |           | inout  | F        |
!! | Grid           | FV3-GFS_Grid_type                                      | derived type GFS_grid_type in FV3                       | DDT           |    0 | GFS_grid_type     |           | in     | F        |
!! | Tbd            | FV3-GFS_Tbd_type                                       | derived type GFS_tbd_type in FV3                        | DDT           |    0 | GFS_tbd_type      |           | in     | F        |
!! | Cldprop        | FV3-GFS_Cldprop_type                                   | derived type GFS_cldprop_type in FV3                    | DDT           |    0 | GFS_cldprop_type  |           | in     | F        |
!! | Radtend        | FV3-GFS_Radtend_type                                   | derived type GFS_radtend_type in FV3                    | DDT           |    0 | GFS_radtend_type  |           | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                      | derived type GFS_diag_type in FV3                       | DDT           |    0 | GFS_diag_type     |           | inout  | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                | none          |    0 | character         | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                   | flag          |    0 | integer           |           | out    | F        |
!!
!-------------------------------------------------------------------------
! GFS stochastic_driver
!-------------------------------------------------------------------------
!    routine called prior to radiation and physics steps to handle:
!      1) sets up various time/date variables
!      2) sets up various triggers
!      3) defines random seed indices for radiation (in a reproducible way)
!      5) interpolates coefficients for prognostic ozone calculation
!      6) performs surface data cycling via the GFS gcycle routine
!-------------------------------------------------------------------------
      subroutine GFS_stochastics_run (Model, Statein, Stateout, Sfcprop, Coupling, &
                                      Grid, Tbd, Cldprop, Radtend, Diag, errmsg, errflg)

         use machine,               only: kind_phys
         use GFS_typedefs,          only: GFS_control_type, GFS_statein_type,  &
                                          GFS_stateout_type, GFS_sfcprop_type, &
                                          GFS_coupling_type, GFS_grid_type,    &
                                          GFS_tbd_type, GFS_cldprop_type,      &
                                          GFS_radtend_type, GFS_diag_type

         implicit none

         !--- interface variables
         type(GFS_control_type),   intent(in   ) :: Model
         type(GFS_statein_type),   intent(in   ) :: Statein
         type(GFS_stateout_type),  intent(in   ) :: Stateout
         type(GFS_sfcprop_type),   intent(in   ) :: Sfcprop
         type(GFS_coupling_type),  intent(inout) :: Coupling
         type(GFS_grid_type),      intent(in   ) :: Grid
         type(GFS_tbd_type),       intent(in   ) :: Tbd
         type(GFS_cldprop_type),   intent(in   ) :: Cldprop
         type(GFS_radtend_type),   intent(in   ) :: Radtend
         type(GFS_diag_type),      intent(inout) :: Diag
         character(len=*),         intent(out)   :: errmsg
         integer,                  intent(out)   :: errflg

         !--- local variables
         integer :: k, i
         real(kind=kind_phys) :: upert, vpert, tpert, qpert, qnew

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (Model%do_sppt) then
           do k = 1,size(Statein%tgrs,2)
             do i = 1,size(Statein%tgrs,1)

               upert = (Stateout%gu0(i,k)   - Statein%ugrs(i,k))   * Coupling%sppt_wts(i,k)
               vpert = (Stateout%gv0(i,k)   - Statein%vgrs(i,k))   * Coupling%sppt_wts(i,k)
               tpert = (Stateout%gt0(i,k)   - Statein%tgrs(i,k))   * Coupling%sppt_wts(i,k) - Tbd%dtdtr(i,k)
               qpert = (Stateout%gq0(i,k,1) - Statein%qgrs(i,k,1)) * Coupling%sppt_wts(i,k)

               Stateout%gu0(i,k)  = Statein%ugrs(i,k)+upert
               Stateout%gv0(i,k)  = Statein%vgrs(i,k)+vpert

               !negative humidity check
               qnew = Statein%qgrs(i,k,1)+qpert
               if (qnew .GE. 1.0e-10) then
                  Stateout%gq0(i,k,1) = qnew
                  Stateout%gt0(i,k)   = Statein%tgrs(i,k) + tpert + Tbd%dtdtr(i,k)
               endif
             enddo
           enddo

           Diag%totprcp(:)      = Diag%totprcp(:)      + (Coupling%sppt_wts(:,15) - 1.0)*Tbd%dtotprcp(:)
           Diag%cnvprcp(:)      = Diag%cnvprcp(:)      + (Coupling%sppt_wts(:,15) - 1.0)*Tbd%dcnvprcp(:)
           Coupling%rain_cpl(:) = Coupling%rain_cpl(:) + (Coupling%sppt_wts(:,15) - 1.0)*Tbd%drain_cpl(:)
           Coupling%snow_cpl(:) = Coupling%snow_cpl(:) + (Coupling%sppt_wts(:,15) - 1.0)*Tbd%dsnow_cpl(:)
         endif

         if (Model%do_shum) then
           Stateout%gq0(:,:,1) = Stateout%gq0(:,:,1)*(1.0 + Coupling%shum_wts(:,:))
         endif

      end subroutine GFS_stochastics_run

    end module GFS_stochastics
