!
!! ! Module: gfs_phy_tracer_config
!
! ! Description: gfs physics tracer configuration module
!
! ! Revision history:
!   Oct 16 2009   Sarah Lu, adopted from dyn fc
!   Nov 21 2009   Sarah Lu, chem tracer specified from ChemRegistry
!   Dec 10 2009   Sarah Lu, add doing_GOCART
!   Jan 12 2010   Sarah Lu, add trcindx
!   Feb 08 2009   Sarah Lu, ri/cpi added to gfs_phy_tracer_type
!   Aug 17 2010   Sarah Lu, remove debug print
!   Oct 16 2010   Sarah Lu, add fscav
!   Aug 08 2011   Jun Wang, remove gocart dependency when not running GOCART
!   Sep 17 2011   Sarah Lu, revise chem tracer initialization
!   Nov 11 2011   Sarah Lu, allocate but not assign value for cpi/ri array
!   Apr 06 2012   Henry Juang, relax hardwire num_tracer, add tracer 4 and 5
!   Apr 23 2012   Jun Wang, remove save attibute for gfs_phy_tracer (already defined)
!   --- -- 2016   Anning Cheng add ntiw,ntlnc,ntinc
!   May 03 2016   S Moorthi add nto, nto2
! -------------------------------------------------------------------------
!
      module gfs_phy_tracer_config
      use machine , only     : kind_phys

      implicit none
      SAVE
!
! tracer specification: add fscav 
!
      type    gfs_phy_tracer_type
        character*20        , pointer      :: chem_name(:) ! chem_tracer name
        character*20        , pointer      :: vname(:)     ! variable name
        real(kind=kind_phys), pointer      :: ri(:)
        real(kind=kind_phys), pointer      :: cpi(:)
        real(kind=kind_phys), pointer      :: fscav(:)    
        integer                  :: ntrac,    ntrac_met, ntrac_chem
        logical                  :: doing_DU, doing_SU,  doing_SS
     &,                             doing_OC, doing_BC,  doing_GOCART
      endtype gfs_phy_tracer_type

      type (gfs_phy_tracer_type) ::  gfs_phy_tracer
!
! misc tracer options
!
      logical                    :: glbsum  = .true.
!

! --- public interface
      public     tracer_config_init, trcindx

      contains

! -------------------------------------------------------------------   
! -------------------------------------------------------------------   
!     subroutine tracer_config_init (gfs_phy_tracer,ntrac,
      subroutine tracer_config_init (ntrac,ntoz,ntcw,ncld,
     &                               ntiw,ntlnc,ntinc,
     &                               fprcp,ntrw,ntsw,ntrnc,ntsnc,
     &                               ntke,nto,nto2,me)

c  
c  This subprogram sets up gfs_phy_tracer
c 
      implicit none
! input
      integer, intent(in)     :: me, ntoz,ntcw,ncld,ntke,
     &                           ntiw,ntlnc,ntinc,nto,nto2,
     &                           fprcp,ntrw,ntsw,ntrnc,ntsnc
! output
!      type (gfs_phy_tracer_type), intent(out)    ::  gfs_phy_tracer
! input/output
      integer, intent(inout)  :: ntrac
! local
      integer                 :: i, j, status, ierr
      character*20            :: rgname

! initialize ntrac_chem (the default is no chemistry)
      gfs_phy_tracer%ntrac_chem   = 0
      gfs_phy_tracer%doing_GOCART = .false.

! initialize chem tracers
      call gocart_tracer_config(me)
!     call gocart_tracer_config(gfs_phy_tracer,me)

! ntrac_met = number of met tracers
!hmhj if ( ntoz < ntcw ) then                       
!hmhj   gfs_phy_tracer%ntrac_met = ntcw + ncld - 1   
!hmhj else                                                           
!hmhj   gfs_phy_tracer%ntrac_met = ntoz                              
!hmhj endif                                          
!hmhj if ( gfs_phy_tracer%ntrac_met /= ntrac ) then
!hmhj   print *,'LU_TRC: ERROR ! inconsistency in ntrac:',
!hmhj&           ntrac, gfs_phy_tracer%ntrac_met
!hmhj   stop 222   
!hmhj endif
! input ntrac is meteorological tracers
      gfs_phy_tracer%ntrac_met = ntrac

! update ntrac = total number of tracers
      gfs_phy_tracer%ntrac = gfs_phy_tracer%ntrac_met +     
     &                       gfs_phy_tracer%ntrac_chem
      ntrac = gfs_phy_tracer%ntrac

      if(me==0) then
       print *, 'LU_TRCp: ntrac_met =',gfs_phy_tracer%ntrac_met
       print *, 'LU_TRCp: ntrac_chem=',gfs_phy_tracer%ntrac_chem
       print *, 'LU_TRCp: ntrac     =',gfs_phy_tracer%ntrac
      endif

! Set up tracer name, cpi, and ri
      if ( gfs_phy_tracer%ntrac > 0 ) then      
        allocate(gfs_phy_tracer%vname(ntrac), stat=status)
        if( status /= 0 ) then
          print *,'LU_TRC: alloc error - gfs_dyn_tracer :',status,me
          return
         endif
        allocate(gfs_phy_tracer%ri(0:ntrac),  stat=status)
        if( status /= 0 ) then
          print *,'LU_TRC: alloc error - gfs_dyn_tracer :',status,me
          return
        endif
        allocate(gfs_phy_tracer%cpi(0:ntrac), stat=status)
        if( status /= 0 ) then
          print *,'LU_TRC: alloc error - gfs_dyn_tracer :',status,me
          return
        endif
        allocate(gfs_phy_tracer%fscav(ntrac), stat=status)
        if( status /= 0 ) then
          print *,'LU_TRC: alloc error - gfs_dyn_tracer :',status,me
          return
        endif

!--- fill in met tracers
                      gfs_phy_tracer%vname(1)     = 'spfh'   
        if(ntoz  > 0) gfs_phy_tracer%vname(ntoz)  = 'o3mr'   
        if(ntcw  > 0) gfs_phy_tracer%vname(ntcw)  = 'clwmr'   
        if(ntiw  > 0) gfs_phy_tracer%vname(ntiw)  = 'climr'   
        if(ntlnc > 0) gfs_phy_tracer%vname(ntlnc) = 'lnc'
        if(ntinc > 0) gfs_phy_tracer%vname(ntinc) = 'inc'
        if(ntrw  > 0) gfs_phy_tracer%vname(ntrw)  = 'rnmr'
        if(ntsw  > 0) gfs_phy_tracer%vname(ntsw)  = 'snwmr'
        if(ntrnc > 0) gfs_phy_tracer%vname(ntrnc) = 'rnc'
        if(ntsnc > 0) gfs_phy_tracer%vname(ntsnc) = 'snc'
        if(ntke  > 0) gfs_phy_tracer%vname(ntke)  = 'tke'   
        if(nto   > 0) gfs_phy_tracer%vname(nto)   = 'o'   
        if(nto2  > 0) gfs_phy_tracer%vname(nto2)  = 'o2'   


        gfs_phy_tracer%fscav(1:gfs_phy_tracer%ntrac_met) = 0.

!--- fill in chem tracers
        if ( gfs_phy_tracer%ntrac_chem > 0 ) then      
          do i = 1,gfs_phy_tracer%ntrac_chem
            j = i + gfs_phy_tracer%ntrac_met
            rgname = trim(gfs_phy_tracer%chem_name(i))
            if(me==0)print *, 'LU_TRC_phy: vname=',j,rgname
                         gfs_phy_tracer%vname(j) = rgname
          enddo
        endif

      endif     !!

      return

      end subroutine tracer_config_init
! -------------------------------------------------------------------
! -------------------------------------------------------------------
      function trcindx( specname, tracer )
      implicit none

      character*(*), intent(in)  ::  specname
      type (gfs_phy_tracer_type), intent(in)    ::  tracer

      character*10  ::  name1, name2
      integer       ::  i, trcindx

! -- set default value
      trcindx = -999

! -- convert specname to upper case
      call fixchar(specname, name1, 1)
      do i = 1, tracer%ntrac
        call fixchar(tracer%vname(i), name2, 1)
        if( name1 == name2 ) then
          trcindx = i 
          exit
        endif
      enddo

      return
      end function trcindx

! -------------------------------------------------------------------
      subroutine fixchar(name_in, name_out, option)
      implicit none

      character*(*), intent(in)   ::  name_in
      character*(*), intent(out)  ::  name_out
      integer, intent(in)         ::  option

      character*10                :: temp
      integer                     :: i, ic

      name_out= '          '
      temp = trim(adjustl(name_in))
      do i = 1, len_trim(temp)
        ic = IACHAR(temp(i:i))
        if(option == 1 ) then             !<--- convert to upper case
          if(ic .ge. 97 .and. ic .le. 122) then
            name_out(i:i) = CHAR( IC-32 )
          else
            name_out(i:i) = temp(i:i)
          endif
        endif
        if(option == 2 ) then             !<--- convert to lower case
          if(ic .ge. 65 .and. ic .le. 90) then
            name_out(i:i) = CHAR( IC+32 )
          else
            name_out(i:i) = temp(i:i)
          endif
        endif

      enddo
      name_out=trim(name_out)
      return

      end subroutine fixchar

! ========================================================================= 

      end module gfs_phy_tracer_config
