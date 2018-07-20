module get_prs_fv3

   use machine,  only: kind_phys
   use physcons, only: con_fvirt

!--- public declarations
   public get_prs_fv3_init, get_prs_fv3_run, get_prs_fv3_finalize

!--- local variables
   real(kind=kind_phys), parameter :: zero = 0.0_kind_phys
   real(kind=kind_phys), parameter :: half = 0.5_kind_phys

contains


!! \section arg_table_get_prs_fv3_init Argument Table
!!
   subroutine get_prs_fv3_init()
   end subroutine get_prs_fv3_init


!! \section arg_table_get_prs_fv3_run Argument Table
!! | local_name     | standard_name                                                                     | long_name                                                                           | units      | rank | type      | kind      | intent | optional |
!! |----------------|-----------------------------------------------------------------------------------|-------------------------------------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | ix             | horizontal_dimension                                                              | horizontal dimension                                                                | count      |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                                                | number of vertical layers                                                           | count      |    0 | integer   |           | in     | F        |
!! | phii           | geopotential_at_interface                                                         | interface geopotential                                                              | m2 s-2     |    2 | real      | kind_phys | in     | F        |
!! | prsi           | air_pressure_at_interface                                                         | interface pressure                                                                  | Pa         |    2 | real      | kind_phys | in     | F        |
!! | tgrs           | air_temperature                                                                   | mid-layer temperature                                                               | K          |    2 | real      | kind_phys | in     | F        |
!! | qgrs1          | water_vapor_specific_humidity                                                     | mid-layer specific humidity of water vapor                                          | kg kg-1    |    2 | real      | kind_phys | in     | F        |
!! | del            | air_pressure_difference_between_midlayers                                         | difference between mid-layer pressures                                              | Pa         |    2 | real      | kind_phys | out    | F        |
!! | del_gz         | geopotential_difference_between_midlayers_divided_by_midlayer_virtual_temperature | difference between mid-layer geopotentials divided by mid-layer virtual temperature | m2 s-2 K-1 |    2 | real      | kind_phys | out    | F        |
!! | errmsg         | error_message                                                                     | error message for error handling in CCPP                                            | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                                        | error flag for error handling in CCPP                                               | flag       |    0 | integer   |           | out    | F        |
!!
   subroutine get_prs_fv3_run(ix, levs, phii, prsi, tgrs, qgrs1, del, del_gz, errmsg, errflg)

     implicit none

     ! Interface variables
     integer, intent(in) :: ix, levs
     real(kind=kind_phys), dimension(ix,levs+1),     intent(in)    :: phii
     real(kind=kind_phys), dimension(ix,levs+1),     intent(in)    :: prsi
     real(kind=kind_phys), dimension(ix,levs),       intent(in)    :: tgrs
     real(kind=kind_phys), dimension(ix,levs),       intent(in)    :: qgrs1
     real(kind=kind_phys), dimension(ix,levs),       intent(out)   :: del
     real(kind=kind_phys), dimension(ix,levs+1),     intent(out)   :: del_gz
     character(len=*),                               intent(out)   :: errmsg
     integer,                                        intent(out)   :: errflg

     ! Local variables
     integer :: i, k

     ! Initialize CCPP error handling variables
     errmsg = ''
     errflg = 0

! SJL: Adjust the geopotential height hydrostatically in a way consistent with FV3 discretization
! del_gz is a temp array recording the old info before (t,q) are adjusted
     do k=1,levs
       do i=1,ix
            del(i,k) = prsi(i,k) - prsi(i,k+1)
         del_gz(i,k) = (phii(i,k+1) - phii(i,k)) /                    &
                        (tgrs(i,k)*(1.+con_fvirt*max(zero,qgrs1(i,k))))
       enddo
     enddo

   end subroutine get_prs_fv3_run


!! \section arg_table_get_prs_fv3_finalize Argument Table
!!
   subroutine get_prs_fv3_finalize()
   end subroutine get_prs_fv3_finalize


end module get_prs_fv3



module get_phi_fv3

   use machine,  only: kind_phys
   use physcons, only: con_fvirt

!--- public declarations
   public get_phi_fv3_init, get_phi_fv3_run, get_phi_fv3_finalize

!--- local variables
   real(kind=kind_phys), parameter :: zero = 0.0_kind_phys
   real(kind=kind_phys), parameter :: half = 0.5_kind_phys

contains

!! \section arg_table_get_phi_fv3_init Argument Table
!!
   subroutine get_phi_fv3_init()
   end subroutine get_phi_fv3_init


!! \section arg_table_get_phi_fv3_run Argument Table
!! | local_name     | standard_name                                                                     | long_name                                                                           | units      | rank | type      | kind      | intent | optional |
!! |----------------|-----------------------------------------------------------------------------------|-------------------------------------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | ix             | horizontal_dimension                                                              | horizontal dimension                                                                | count      |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                                                | number of vertical layers                                                           | count      |    0 | integer   |           | in     | F        |
!! | gt0            | air_temperature_updated_by_physics                                                | updated air temperature                                                             | K          |    2 | real      | kind_phys | in     | F        |
!! | gq01           | water_vapor_specific_humidity_updated_by_physics                                  | mid-layer specific humidity of water vapor                                          | kg kg-1    |    2 | real      | kind_phys | in     | F        |
!! | del_gz         | geopotential_difference_between_midlayers_divided_by_midlayer_virtual_temperature | difference between mid-layer geopotentials divided by mid-layer virtual temperature | m2 s-2 K-1 |    2 | real      | kind_phys | inout  | F        |
!! | phii           | geopotential_at_interface                                                         | interface geopotential                                                              | m2 s-2     |    2 | real      | kind_phys | out    | F        |
!! | phil           | geopotential                                                                      | mid-layer geopotential                                                              | m2 s-2     |    2 | real      | kind_phys | out    | F        |
!! | errmsg         | error_message                                                                     | error message for error handling in CCPP                                            | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                                        | error flag for error handling in CCPP                                               | flag       |    0 | integer   |           | out    | F        |
!!
   subroutine get_phi_fv3_run(ix, levs, gt0, gq01, del_gz, phii, phil, errmsg, errflg)

     implicit none

     ! Interface variables
     integer, intent(in) :: ix, levs
     real(kind=kind_phys), dimension(ix,levs),       intent(in)    :: gt0
     real(kind=kind_phys), dimension(ix,levs),       intent(in)    :: gq01
     real(kind=kind_phys), dimension(ix,levs+1),     intent(inout) :: del_gz
     real(kind=kind_phys), dimension(ix,levs+1),     intent(out)   :: phii
     real(kind=kind_phys), dimension(ix,levs),       intent(out)   :: phil
     character(len=*),                               intent(out)   :: errmsg
     integer,                                        intent(out)   :: errflg

     ! Local variables
     integer :: i, k

     ! Initialize CCPP error handling variables
     errmsg = ''
     errflg = 0

! SJL: Adjust the height hydrostatically in a way consistent with FV3 discretization
     do i=1,ix
        phii(i,1) = zero
     enddo
     do k=1,levs
       do i=1,ix
         del_gz(i,k) = del_gz(i,k)*gt0(i,k) *                          &
     &                 (1.+con_fvirt*max(zero,gq01(i,k)))
         phii(i,k+1) = phii(i,k) + del_gz(i,k)
         phil(i,k)   = half*(phii(i,k) + phii(i,k+1))
       enddo
     enddo

   end subroutine get_phi_fv3_run


!! \section arg_table_get_phi_fv3_finalize Argument Table
!!
   subroutine get_phi_fv3_finalize()
   end subroutine get_phi_fv3_finalize


end module get_phi_fv3


