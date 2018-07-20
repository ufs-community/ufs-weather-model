      subroutine read_h2odata (h2o_phys, me, master)
      use machine,  only: kind_phys
      use h2o_def
!--- in/out
      logical, intent(in) :: h2o_phys
      integer, intent(in) :: me
      integer, intent(in) :: master
!--- locals
      integer :: i, n, k
      real(kind=4), allocatable, dimension(:) :: h2o_lat4, h2o_pres4
      real(kind=4), allocatable, dimension(:) :: h2o_time4, tempin

      if (.not. h2o_phys) then
        latsh2o   = 1 
        levh2o    = 1 
        h2o_coeff = 1
        timeh2o   = 1

        return
      endif

      open(unit=kh2opltc,file='INPUT/global_h2oprdlos.f77', form='unformatted', convert='big_endian')

!--- read in indices
!---
      read (kh2opltc) h2o_coeff, latsh2o, levh2o, timeh2o
      if (me == master) then
        write(*,*) 'Reading in h2odata from global_h2oprdlos.f77 '
        write(*,*) '     h2o_coeff = ', h2o_coeff
        write(*,*) '       latsh2o = ', latsh2o
        write(*,*) '        levh2o = ', levh2o
        write(*,*) '       timeh2o = ', timeh2o
      endif

!--- read in data
!---   h2o_lat   -  latitude of data        (-90 to 90)
!---   h2o_pres  -  vertical pressure level (mb)
!---   h2o_time  -  time coordinate         (days)
!---
      allocate (h2o_lat(latsh2o), h2o_pres(levh2o),h2o_time(timeh2o+1))
      allocate (h2o_lat4(latsh2o), h2o_pres4(levh2o),h2o_time4(timeh2o+1))
      rewind (kh2opltc)
      read (kh2opltc) h2o_coeff, latsh2o, levh2o, timeh2o, h2o_lat4, h2o_pres4, h2o_time4
      h2o_pres(:) = h2o_pres4(:)
!---  convert pressure levels from mb to ln(Pa)
      h2o_pres(:) = log(100.0*h2o_pres(:))
      h2o_lat(:)  = h2o_lat4(:)
      h2o_time(:) = h2o_time4(:)
      deallocate (h2o_lat4, h2o_pres4, h2o_time4)

!--- read in h2oplin which is in order of (lattitudes, water levels, coeff number, time)
!--- assume latitudes is on a uniform gaussian grid
!---
      allocate (tempin(latsh2o))
      allocate (h2oplin(latsh2o,levh2o,h2o_coeff,timeh2o))
      DO i=1,timeh2o
        do n=1,h2o_coeff
          DO k=1,levh2o
            READ(kh2opltc) tempin
            h2oplin(:,k,n,i) = tempin(:)
          ENDDO
        enddo
      ENDDO
      deallocate (tempin)

      close(kh2opltc)

      end subroutine read_h2odata
!
!**********************************************************************
!
      subroutine setindxh2o(npts,dlat,jindx1,jindx2,ddy)
!
! May 2015 Shrinivas Moorthi - Prepare for H2O interpolation
!
      use machine, only: kind_phys
      use h2o_def, only: jh2o => latsh2o, h2o_lat, h2o_time
!
      implicit none
!
      integer                     npts
      integer, dimension(npts) :: jindx1, jindx2
      real(kind=kind_phys)     :: dlat(npts),ddy(npts)
!
      integer i,j,lat
!
      do j=1,npts
        jindx2(j) = jh2o + 1
        do i=1,jh2o
          if (dlat(j) < h2o_lat(i)) then
            jindx2(j) = i
            exit
          endif
        enddo
        jindx1(j) = max(jindx2(j)-1,1)
        jindx2(j) = min(jindx2(j),jh2o)
        if (jindx2(j) /= jindx1(j)) then
          ddy(j) = (dlat(j)            - h2o_lat(jindx1(j))) &
                 / (h2o_lat(jindx2(j)) - h2o_lat(jindx1(j)))
        else
          ddy(j) = 1.0
        endif
!       print *,' j=',j,' dlat=',dlat(j),' jindx12=',jindx1(j), &
!          jindx2(j),' h2o_lat=',h2o_lat(jindx1(j)),          &
!          h2o_lat(jindx2(j)),' ddy=',ddy(j)
      enddo
 
      return
      end
!
!**********************************************************************
!
      subroutine h2ointerpol(me,npts,idate,fhour,jindx1,jindx2,h2oplout,ddy)
!
! May 2015 Shrinivas Moorthi - Prepare for H2O interpolation
!
      use machine , only : kind_phys
      use h2o_def
      implicit none
      integer             j,j1,j2,l,npts,nc,n1,n2
      real(kind=kind_phys) fhour,tem, tx1, tx2
!
 
      integer  jindx1(npts), jindx2(npts)
      integer  me,idate(4)
      integer  idat(8),jdat(8)
!
      real(kind=kind_phys) ddy(npts)
      real(kind=kind_phys) h2oplout(levh2o,npts,h2o_coeff)
      real(kind=kind_phys) rinc(5), rjday
      integer              jdow, jdoy, jday
      real(4)              rinc4(5)
      integer              w3kindreal, w3kindint
!
      idat    = 0
      idat(1) = idate(4)
      idat(2) = idate(2)
      idat(3) = idate(3)
      idat(5) = idate(1)
      rinc    = 0.
      rinc(2) = fhour
      call w3kind(w3kindreal,w3kindint)
      if(w3kindreal==4) then
        rinc4 = rinc
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
      if (rjday < h2o_time(1)) rjday = rjday+365.
!
      n2 = timeh2o + 1
      do j=1,timeh2o
        if (rjday < h2o_time(j)) then
          n2 = j
          exit
        endif
      enddo
      n1 = n2 - 1
      if (n1 <= 0)      n1 = n1 + timeh2o
      if (n2 > timeh2o) n2 = n2 - timeh2o

!
!     if (me .eq. 0) print *,' n1=',n1,' n2=',n2,' rjday=',rjday
!    &,'h2o_time=',h2o_time(n1),h2o_time(n2)
!

      tx1 = (h2o_time(n2) - rjday) / (h2o_time(n2) - h2o_time(n1))
      tx2 = 1.0 - tx1
!
      do nc=1,h2o_coeff
        do l=1,levh2o
          do j=1,npts
            j1  = jindx1(j)
            j2  = jindx2(j)
            tem = 1.0 - ddy(j)
            h2oplout(j,l,nc) = & 
              tx1*(tem*h2oplin(j1,l,nc,n1)+ddy(j)*h2oplin(j2,l,nc,n1)) & 
            + tx2*(tem*h2oplin(j1,l,nc,n2)+ddy(j)*h2oplin(j2,l,nc,n2))
          enddo
        enddo
      enddo
!
      return
      end
