!>\file tridi.f
!! These subroutines are originally internal subroutines in moninedmf.f

!>\ingroup GFS_edmf_main
!!\brief Routine to solve the tridiagonal system to calculate
!!temperature and moisture at \f$ t + \Delta t \f$; part of two-part
!!process to calculate time tendencies due to vertical diffusion.
      subroutine tridi2(l,n,cl,cm,cu,r1,r2,au,a1,a2)
cc
      use machine     , only : kind_phys
      implicit none
      integer             k,n,l,i
      real(kind=kind_phys) fk
cc
      real(kind=kind_phys) cl(l,2:n),cm(l,n),cu(l,n-1),r1(l,n),r2(l,n), &
     &          au(l,n-1),a1(l,n),a2(l,n)
c----------------------------------------------------------------------
      do i=1,l
        fk      = 1./cm(i,1)
        au(i,1) = fk*cu(i,1)
        a1(i,1) = fk*r1(i,1)
        a2(i,1) = fk*r2(i,1)
      enddo
      do k=2,n-1
        do i=1,l
          fk      = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k) = fk*cu(i,k)
          a1(i,k) = fk*(r1(i,k)-cl(i,k)*a1(i,k-1))
          a2(i,k) = fk*(r2(i,k)-cl(i,k)*a2(i,k-1))
        enddo
      enddo
      do i=1,l
        fk      = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk*(r1(i,n)-cl(i,n)*a1(i,n-1))
        a2(i,n) = fk*(r2(i,n)-cl(i,n)*a2(i,n-1))
      enddo
      do k=n-1,1,-1
        do i=1,l
          a1(i,k) = a1(i,k)-au(i,k)*a1(i,k+1)
          a2(i,k) = a2(i,k)-au(i,k)*a2(i,k+1)
        enddo
      enddo
c-----------------------------------------------------------------------
      return
      end

c-----------------------------------------------------------------------
!>\ingroup GFS_edmf_main
!!  \brief Routine to solve the tridiagonal system to calculate u- and
!! v-momentum at \f$ t + \Delta t \f$; part of two-part process to
!! calculate time tendencies due to vertical diffusion.
      subroutine tridin(l,n,nt,cl,cm,cu,r1,r2,au,a1,a2)
cc
      use machine     , only : kind_phys
      implicit none
      integer             is,k,kk,n,nt,l,i
      real(kind=kind_phys) fk(l)
cc
      real(kind=kind_phys) cl(l,2:n), cm(l,n), cu(l,n-1),               &
     &                     r1(l,n),   r2(l,n*nt),                       &
     &                     au(l,n-1), a1(l,n), a2(l,n*nt),              &
     &                     fkk(l,2:n-1)
c-----------------------------------------------------------------------
      do i=1,l
        fk(i)   = 1./cm(i,1)
        au(i,1) = fk(i)*cu(i,1)
        a1(i,1) = fk(i)*r1(i,1)
      enddo
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          a2(i,1+is) = fk(i) * r2(i,1+is)
        enddo
      enddo
      do k=2,n-1
        do i=1,l
          fkk(i,k) = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)  = fkk(i,k)*cu(i,k)
          a1(i,k)  = fkk(i,k)*(r1(i,k)-cl(i,k)*a1(i,k-1))
        enddo
      enddo
      do kk = 1, nt
        is = (kk-1) * n
        do k=2,n-1
          do i=1,l
            a2(i,k+is) = fkk(i,k)*(r2(i,k+is)-cl(i,k)*a2(i,k+is-1))
          enddo
        enddo
      enddo
      do i=1,l
        fk(i)   = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk(i)*(r1(i,n)-cl(i,n)*a1(i,n-1))
      enddo
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          a2(i,n+is) = fk(i)*(r2(i,n+is)-cl(i,n)*a2(i,n+is-1))
        enddo
      enddo
      do k=n-1,1,-1
        do i=1,l
          a1(i,k) = a1(i,k) - au(i,k)*a1(i,k+1)
        enddo
      enddo
      do kk = 1, nt
        is = (kk-1) * n
        do k=n-1,1,-1
          do i=1,l
            a2(i,k+is) = a2(i,k+is) - au(i,k)*a2(i,k+is+1)
          enddo
        enddo
      enddo
c-----------------------------------------------------------------------
      return
      end
