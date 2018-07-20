      SUBROUTINE read_o3data (ntoz, me, master)
      use machine,  only: kind_phys
      use ozne_def
!--- in/out
      integer, intent(in) :: ntoz
      integer, intent(in) :: me
      integer, intent(in) :: master
!--- locals
      integer :: i, n, k
      real(kind=4), allocatable, dimension(:) :: oz_lat4, oz_pres4
      real(kind=4), allocatable, dimension(:) :: oz_time4, tempin

      if (ntoz <= 0) then      ! Diagnostic ozone
        rewind (kozc)
        read (kozc,end=101) latsozc, levozc, timeozc, blatc4
  101   if (levozc  < 10 .or. levozc > 100) then
          rewind (kozc)
          levozc  = 17
          latsozc = 18
          blatc   = -85.0
        else
          blatc   = blatc4
        endif
        latsozp   = 2
        levozp    = 1
        timeoz    = 1
        oz_coeff  = 0
        dphiozc = -(blatc+blatc)/(latsozc-1)
        return
      endif

      open(unit=kozpl,file='INPUT/global_o3prdlos.f77', form='unformatted', convert='big_endian')

!--- read in indices
!---
      read (kozpl) oz_coeff, latsozp, levozp, timeoz
      if (me == master) then
        write(*,*) 'Reading in o3data from global_o3prdlos.f77 '
        write(*,*) '      oz_coeff = ', oz_coeff
        write(*,*) '       latsozp = ', latsozp
        write(*,*) '        levozp = ', levozp
        write(*,*) '        timeoz = ', timeoz
      endif

!--- read in data
!---   oz_lat   -  latitude of data        (-90 to 90)
!---   oz_pres  -  vertical pressure level (mb)
!---   oz_time  -  time coordinate         (days)
!---
      allocate (oz_lat(latsozp), oz_pres(levozp),oz_time(timeoz+1))
      allocate (oz_lat4(latsozp), oz_pres4(levozp),oz_time4(timeoz+1))
      rewind (kozpl)
      read (kozpl) oz_coeff, latsozp, levozp, timeoz, oz_lat4, oz_pres4, oz_time4
      oz_pres(:) = oz_pres4(:)
!---  convert pressure levels from mb to ln(Pa)
      oz_pres(:) = log(100.0*oz_pres(:))
      oz_lat(:)  = oz_lat4(:)
      oz_time(:) = oz_time4(:)
      deallocate (oz_lat4, oz_pres4, oz_time4)

!--- read in ozplin which is in order of (lattitudes, ozone levels, coeff number, time)
!--- assume latitudes is on a uniform gaussian grid
!---
      allocate (tempin(latsozp))
      allocate (ozplin(latsozp,levozp,oz_coeff,timeoz))
      DO i=1,timeoz
        DO n=1,oz_coeff
          DO k=1,levozp
            READ(kozpl) tempin
            ozplin(:,k,n,i) = tempin(:)
          ENDDO
        ENDDO
      ENDDO
      deallocate (tempin)

      close(kozpl)

      END SUBROUTINE read_o3data
!
!**********************************************************************
!
      SUBROUTINE setindxoz(npts,dlat,jindx1,jindx2,ddy)
!
      USE MACHINE,  ONLY: kind_phys
      USE OZNE_DEF, ONLY: jo3 => latsozp, oz_lat
!
      implicit none
!
      integer npts, JINDX1(npts),JINDX2(npts)
      real(kind=kind_phys) dlat(npts),DDY(npts)
!
      integer i,j,lat
!
      DO J=1,npts
        jindx2(j) = jo3 + 1
        do i=1,jo3
          if (dlat(j) < oz_lat(i)) then
            jindx2(j) = i
            exit
          endif
        enddo
        jindx1(j) = max(jindx2(j)-1,1)
        jindx2(j) = min(jindx2(j),jo3)
        if (jindx2(j) .ne. jindx1(j)) then
          DDY(j) = (dlat(j)           - oz_lat(jindx1(j))) &
                 / (oz_lat(jindx2(j)) - oz_lat(jindx1(j)))
        else
          ddy(j) = 1.0
        endif
!       print *,' j=',j,' dlat=',dlat(j),' jindx12=',jindx1(j), &
!         jjindx2(j),' oz_lat=',oz_lat(jindx1(j)),              &
!         oz_lat(jindx2(j)),' ddy=',ddy(j)
      ENDDO
 
      RETURN
      END
!
!**********************************************************************
!
      SUBROUTINE ozinterpol(me,npts,IDATE,FHOUR,jindx1,jindx2,ozplout,ddy)
!
      USE MACHINE,  ONLY : kind_phys
      USE OZNE_DEF
      implicit none
      integer             iday,j,j1,j2,l,npts,nc,n1,n2
      real(kind=kind_phys) fhour,tem, tx1, tx2
!
 
      integer  JINDX1(npts), JINDX2(npts)
      integer  me,idate(4)
      integer  IDAT(8),JDAT(8)
!
      real(kind=kind_phys) DDY(npts)
      real(kind=kind_phys) ozplout(npts,levozp,oz_coeff)
      real(kind=kind_phys) RINC(5), rjday
      integer jdow, jdoy, jday
      real(4) rinc4(5)
      integer w3kindreal,w3kindint
!
      IDAT=0
      IDAT(1)=IDATE(4)
      IDAT(2)=IDATE(2)
      IDAT(3)=IDATE(3)
      IDAT(5)=IDATE(1)
      RINC=0.
      RINC(2)=FHOUR
      call w3kind(w3kindreal,w3kindint)
      if(w3kindreal==4) then
        rinc4=rinc
        CALL W3MOVDAT(RINC4,IDAT,JDAT)
      else
        CALL W3MOVDAT(RINC,IDAT,JDAT)
      endif
!
      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(jdat,jdow,jdoy,jday)
      rjday = jdoy + jdat(5) / 24.
      IF (RJDAY .LT. oz_time(1)) RJDAY = RJDAY+365.
!
      n2 = timeoz + 1
      do j=1,timeoz
        if (rjday .lt. oz_time(j)) then
          n2 = j
          exit
        endif
      enddo
      n1 = n2 - 1
      if (n1 <= 0)     n1 = n1 + timeoz
      if (n2 > timeoz) n2 = n2 - timeoz

!
!     if (me .eq. 0) print *,' n1=',n1,' n2=',n2,' rjday=',rjday
!    &,'oz_time=',oz_time(n1),oz_time(n2)
!

      tx1 = (oz_time(n2) - rjday) / (oz_time(n2) - oz_time(n1))
      tx2 = 1.0 - tx1
!
      do nc=1,oz_coeff
        DO L=1,levozp
          DO J=1,npts
            J1  = JINDX1(J)
            J2  = JINDX2(J)
            TEM = 1.0 - DDY(J)
            ozplout(j,L,nc) = & 
            tx1*(TEM*ozplin(J1,L,nc,n1)+DDY(J)*ozplin(J2,L,nc,n1)) & 
          + tx2*(TEM*ozplin(J1,L,nc,n2)+DDY(J)*ozplin(J2,L,nc,n2))
          ENDDO
        ENDDO
      enddo
!
      RETURN
      END
