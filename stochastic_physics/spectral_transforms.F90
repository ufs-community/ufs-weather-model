!>@brief The module 'spectral_transforms' contains the subroutines spec_to_four and four_to_grid
module spectral_transforms

 use kinddef
 use mpi_wrapper, only : mp_alltoall,mype,npes
 use stochy_internal_state_mod, only : stochy_internal_state
 use stochy_namelist_def

      private 
      public :: spec_to_four, four_to_grid,dozeuv_stochy,dezouv_stochy
      public :: initialize_spectral,stochy_la2ga

      integer, public :: ls_dim,            &
              ls_max_node,       &
              len_trie_ls,       &
              len_trio_ls,       &
              jcap,latg,latg2,   &
              skeblevs,levs,lnt, &
              lonf,lonfx

!
      integer, public, allocatable :: lat1s_a(:), lon_dims_a(:)
      real, public, allocatable, dimension(:) ::  colrad_a, wgt_a, rcs2_a, &
                                       sinlat_a, coslat_a


      contains

!>@brief The subrountine 'spec_to_four' converts the spherical harmonics to fourier coefficients
!>@details This code is taken from the legacy spectral GFS
      subroutine spec_to_four(flnev,flnod,plnev,plnod, &
                               ls_node,                &
                               workdim,four_gr,        &
                               ls_nodes,max_ls_nodes,  &
                               lats_nodes,global_lats, &
                               lats_node,ipt_lats_node, &
                               nvars )
!

      implicit none
!
      external esmf_dgemm
!     
      integer, intent(in)     :: nvars
      real(kind=kind_dbl_prec) flnev(len_trie_ls,2*nvars)
      real(kind=kind_dbl_prec) flnod(len_trio_ls,2*nvars)
!
      real(kind=kind_dbl_prec) plnev(len_trie_ls,latg2)
      real(kind=kind_dbl_prec) plnod(len_trio_ls,latg2)
!
      integer              ls_node(ls_dim,3)
!
!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!
!    local scalars
!    -------------
!
      integer              j, l, lat, lat1, n, kn, n2,indev,indod
!
!    local arrays
!    ------------
!
      real(kind=kind_dbl_prec), dimension(nvars*2,latg2) ::  apev, apod
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      integer              workdim, lats_node, ipt_lats_node
!
      real(kind=kind_dbl_prec) four_gr(lonf+2,nvars,workdim)
!
      integer              ls_nodes(ls_dim,npes)
      integer, dimension(npes) :: max_ls_nodes, lats_nodes
      integer, dimension(latg)  :: global_lats
      real(kind=4),target,dimension(2,nvars,ls_dim*workdim,npes)::  workr,works
      real(kind=4),pointer:: work1dr(:),work1ds(:)
      integer, dimension(npes) :: kpts, kptr, sendcounts, recvcounts,  sdispls
!
      integer              ilat,ipt_ls, lmax,lval,jj,nv
      integer              node,arrsz,my_pe,nvar
      integer              ilat_list(npes)              !    for OMP buffer copy
!
!    statement functions
!    -------------------
!
      integer              indlsev, jbasev, indlsod, jbasod
!
      include 'function_indlsev'
      include 'function_indlsod'
!
      real(kind=kind_dbl_prec), parameter ::  cons0=0.0d0, cons1=1.0d0
!
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      n2=2*nvars
      arrsz=n2*ls_dim*workdim*npes
      kpts   = 0
!
      do j = 1, ls_max_node   ! start of do j loop #####################
!
        l = ls_node(j,1)
        jbasev = ls_node(j,2)
        jbasod = ls_node(j,3)

        indev  = indlsev(l,l)
        indod  = indlsod(l+1,l)
!
        lat1 = lat1s_a(l)

!           compute the even and odd components of the fourier coefficients
!
!           compute the sum of the even real      terms for each level
!           compute the sum of the even imaginary terms for each level
!
        call esmf_dgemm('t', 'n', n2, latg2-lat1+1, (jcap+3-l)/2, &
                        cons1, flnev(indev,1), len_trie_ls, plnev(indev,lat1), &
                        len_trie_ls, cons0,  apev(1,lat1), n2 )
!
!           compute the sum of the odd real      terms for each level
!           compute the sum of the odd imaginary terms for each level
!
        call esmf_dgemm('t', 'n', n2, latg2-lat1+1, (jcap+2-l)/2,  &
                        cons1, flnod(indod,1), len_trio_ls, plnod(indod,lat1), &
                        len_trio_ls, cons0, apod(1,lat1), n2 )
!
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!       compute the fourier coefficients for each level
!       -----------------------------------------------
!
        ilat_list(1) = 0
        do node = 1, npes - 1
          ilat_list(node+1) = ilat_list(node) + lats_nodes(node)
        end do
!$omp parallel do private(node,jj,ilat,lat,ipt_ls,nvar,kn,n2)
        do node=1,npes
          do jj=1,lats_nodes(node)
            ilat  = ilat_list(node) + jj
            lat    = global_lats(ilat)
            ipt_ls = min(lat,latg-lat+1)
            if ( ipt_ls >= lat1s_a(ls_nodes(j,mype+1)) ) then
              kpts(node) = kpts(node) + 1
              kn = kpts(node)
!
              if ( lat <= latg2 ) then
!                                                northern hemisphere
                do nvar=1,nvars
                  n2 = nvar + nvar
                  works(1,nvar,kn,node) = apev(n2-1,ipt_ls) + apod(n2-1,ipt_ls)
                  works(2,nvar,kn,node) = apev(n2,ipt_ls)  + apod(n2,ipt_ls)
                enddo
              else
!                                                southern hemisphere
                do nvar=1,nvars
                  n2 = nvar + nvar
                  works(1,nvar,kn,node) = apev(n2-1,ipt_ls) - apod(n2-1,ipt_ls)
                  works(2,nvar,kn,node) = apev(n2,ipt_ls)   - apod(n2,ipt_ls)
                enddo
              endif
            endif
          enddo
        enddo
!
      enddo   ! end of do j loop #######################################
!
      kptr = 0
      do node=1,npes
         do l=1,max_ls_nodes(node)
            lval = ls_nodes(l,node)+1
            do j=1,lats_node
               lat = global_lats(ipt_lats_node-1+j)
               if ( min(lat,latg-lat+1) >= lat1s_a(lval-1) ) then
                  kptr(node) = kptr(node) + 1
               endif
            enddo
         enddo
      enddo
!
!
!$omp parallel do private(node)
      do node=1,npes
         sendcounts(node) = kpts(node) * n2
         recvcounts(node) = kptr(node) * n2
         sdispls(node)    = (node-1)   * n2 * ls_dim * workdim
      end do
      work1dr(1:arrsz)=>workr
      work1ds(1:arrsz)=>works
      call mp_alltoall(work1ds, sendcounts, sdispls, &
                        work1dr,recvcounts,sdispls)
      nullify(work1dr)
      nullify(work1ds)
!$omp parallel do private(j,lat,lmax,nvar,lval,nv)
      do j=1,lats_node
         lmax = min(jcap,lonf/2)
         n2   = lmax + lmax + 3
         if ( n2 <= lonf+2 ) then
            do nv=1,nvars
               do lval = n2, lonf+2
                  four_gr(lval,nv,j) = cons0
               enddo
            enddo
         endif
      enddo
!
      kptr = 0
!!
!$omp parallel do private(node,l,lval,j,lat,nvar,kn,n2)
      do node=1,npes
        do l=1,max_ls_nodes(node)
          lval = ls_nodes(l,node)+1
          n2   = lval + lval
          do j=1,lats_node
            lat = global_lats(ipt_lats_node-1+j)
            if ( min(lat,latg-lat+1) >= lat1s_a(lval-1) ) then
              kptr(node) = kptr(node) + 1
              kn = kptr(node)
              do nv=1,nvars
                four_gr(n2-1,nv,j) = workr(1,nv,kn,node)
                four_gr(n2,  nv,j) = workr(2,nv,kn,node)
              enddo
            endif
          enddo
        enddo
      enddo
!
      return
      end subroutine spec_to_four 

!>@brief The subroutine 'four_to_grid' calculate real values form fourrier coefficients
!>@details This code is taken from the legacy spectral GFS
      subroutine four_to_grid(syn_gr_a_1,syn_gr_a_2, lon_dim_coef,nvars)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
!!
      integer, intent(in)  :: nvars
      real(kind=kind_dbl_prec)     syn_gr_a_1(lon_dim_coef,nvars)
      real(kind=kind_dbl_prec)     syn_gr_a_2(lonf,nvars)
      integer                  lon_dim_coef
!________________________________________________________
      real(kind=kind_dbl_prec) aux1crs(44002)
      integer                  init
!________________________________________________________


      init      = 1

      call dcrft_stochy(init,                  &
      syn_gr_a_1(:,:)   ,lon_dim_coef,         &
      syn_gr_a_2(:,:)   ,lonf,         &
      lonf, nvars,                             &
      aux1crs,22000,                           &
      aux1crs(22001),20000)

      init = 0
      call dcrft_stochy(init,                  &
           syn_gr_a_1(:,:)   ,lon_dim_coef,    &
           syn_gr_a_2(:,:)   ,lonf,    &
           lonf, nvars,                        &
           aux1crs,22000,                      &
           aux1crs(22001),20000)

      return
      end


      SUBROUTINE dcrft_stochy(init,x,ldx,y,ldy,n,nvars, table,n1,wrk,n2)
 
      implicit none
      integer ,intent(in) :: ldx,ldy,n,nvars
      integer init,n1,n2,i,j
      real x(ldx,nvars),y(ldy,nvars),table(44002),wrk
 
      IF (init.ne.0) THEN
         CALL rffti_stochy(n,table)
      ELSE
         DO j=1,nvars
            y(1,j)=x(1,j)
            DO i=2,n
              y(i,j)=x(i+1,j)
            ENDDO
            CALL rfftb_stochy(n,y(:,j),table)
         ENDDO
      ENDIF
 
      RETURN
      END

!     ******************************************************************
!     ******************************************************************
!     ******                                                      ******
!     ******                        FFTPACK                       ******
!     ******                                                      ******
!     ******************************************************************
!     ******************************************************************
!
      SUBROUTINE RFFTB_STOCHY (N,R,WSAVE)

      implicit none

      real, intent(inout) :: R(:)  
      real, intent(inout) :: WSAVE(44002)

      integer :: N

      IF (N .EQ. 1) RETURN
      CALL RFFTB1_STOCHY (N,R,WSAVE,WSAVE(N+1:),WSAVE(2*N+1:))
      RETURN
      END

      SUBROUTINE RFFTI_STOCHY (N,WSAVE)

      implicit none

      REAL, intent(inout) ::    WSAVE(44002)
      integer :: N

      IF (N .EQ. 1) RETURN
      CALL RFFTI1_STOCHY (N,WSAVE(N+1:),WSAVE(2*N+1:))
      RETURN
      END


      SUBROUTINE RFFTB1_STOCHY (N,C,CH,WA,RFAC)

      implicit none

      integer, intent(in) :: N
      real, intent(inout) :: CH(44002)
      real, intent(inout) :: C(:)
      real, intent(inout) :: WA(:)
      real, intent(inout) :: RFAC(:)

      integer :: NF,NA,L1,IW,IP,L2,IDO,IDL1,IX2,IX3,IX4
      integer :: K1,I

      NF = INT(RFAC(2))
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = INT(RFAC(K1+2))
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL RADB4_STOCHY (IDO,L1,C(1:4*IDO*L1),CH(1:4*IDO*L1),WA(IW:),WA(IX2:),WA(IX3:))
         GO TO 102
  101    CALL RADB4_STOCHY (IDO,L1,CH(1:4*IDO*L1),C(1:4*IDO*L1),WA(IW:),WA(IX2:),WA(IX3:))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL RADB2_STOCHY (IDO,L1,C,CH,WA(IW:))
         GO TO 105
  104    CALL RADB2_STOCHY (IDO,L1,CH,C,WA(IW:))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
         CALL RADB3_STOCHY (IDO,L1,C,CH,WA(IW:),WA(IX2:))
         GO TO 108
  107    CALL RADB3_STOCHY (IDO,L1,CH,C,WA(IW:),WA(IX2:))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
         CALL RADB5_STOCHY (IDO,L1,C,CH,WA(IW:),WA(IX2:),WA(IX3:),WA(IX4:))
         GO TO 111
  110    CALL RADB5_STOCHY (IDO,L1,CH,C,WA(IW:),WA(IX2:),WA(IX3:),WA(IX4:))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL RADBG_STOCHY (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW:))
         GO TO 114
  113    CALL RADBG_STOCHY (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW:))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      DO 117 I=1,N
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END


      SUBROUTINE RFFTI1_STOCHY (N,WA,RFAC)
      
      implicit none
      
      integer, intent(in) :: N
      REAL, intent(inout) :: WA(:)
      REAL, intent(inout) :: RFAC(:) 

      integer :: NTRYH(4)
      integer :: NL,NF, I, J, NQ,NR,LD,FI,IS,ID,L1,L2,IP
      integer :: NTRY, NFM1, K1,II, IB, IDO, IPM, IC
      REAL, parameter :: TPI=6.28318530717959
      real    :: ARG,ARGLD,ARGH, TI2,TI4

      DATA NTRYH(:) /4,2,3,5/


      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF ( (J-4) .LE. 0) THEN
         GOTO 102
      ELSE
         GOTO 103
      ENDIF
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR.EQ.0) THEN
         GO TO 105
      ELSE
         GO TO 101
      ENDIF
  105 NF = NF+1
      RFAC(NF+2) = FLOAT(NTRY)
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         RFAC(IB+2) = RFAC(IB+1)
  106 CONTINUE
      RFAC(3) = 2.
  107 IF (NL .NE. 1) GO TO 104
      RFAC(1) = FLOAT(N)
      RFAC(2) = FLOAT(NF)
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
!OCL NOVREC
      DO 110 K1=1,NFM1
         IP = INT(RFAC(K1+2))
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = FLOAT(LD)*ARGH
            FI = 0
!OCL SCALAR
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END


      SUBROUTINE RADB2_STOCHY (IDO,L1,CC,CH,WA1)

      implicit none
     
      integer, intent(in) :: IDO
      integer, intent(in) :: L1
      real, intent(inout) :: CC(IDO,2,L1)
      real, intent(inout) :: CH(IDO,L1,2)
      real, intent(inout) :: WA1(:)

      integer :: K,I,IC,IDP2
      real    :: TR2,TI2
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(IDO,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(IDO,2,K)
  101 CONTINUE
      IF ( (IDO-2) .LT. 0) THEN
         GO TO 107
      ELSE IF (( IDO-2).EQ. 0)THEN
         GO TO 105
      ELSE
        GO TO 102
      ENDIF
  102 IDP2 = IDO+2
!OCL NOVREC
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
            TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
            CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
            TI2 = CC(I,1,K)+CC(IC,2,K)
            CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
            CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(IDO,K,1) = CC(IDO,1,K)+CC(IDO,1,K)
         CH(IDO,K,2) = -(CC(1,2,K)+CC(1,2,K))
  106 CONTINUE
  107 RETURN
      END


      SUBROUTINE RADB3_STOCHY (IDO,L1,CC,CH,WA1,WA2)

      implicit none

      integer, intent(in) :: IDO,L1
      real, intent(inout) :: CC(IDO,3,L1)
      real, intent(inout) :: CH(IDO,L1,3)
      real, intent(inout) :: WA1(:)
      real, intent(inout) :: WA2(:)

      REAL, parameter :: TAUR= -.5
      REAL, parameter :: TAUI=.866025403784439
      integer         :: I,K,IDP2,IC
      real            :: TR2,CR2,TI1,CI2,CR3,CI3,DR2,DR3,DI2,DI3
      real            :: TI2,TI4


      DO 101 K=1,L1
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
!OCL NOVREC
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,3,K)-CC(IC,2,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
  102    CONTINUE
  103 CONTINUE
      RETURN
      END


      SUBROUTINE RADB4_STOCHY (IDO,L1,CC,CH,WA1,WA2,WA3)

      implicit none

      integer, intent(in) :: IDO,L1
      real, intent(inout) :: CC(IDO,4,L1)
      real, intent(inout) :: CH(IDO,L1,4)
      real, intent(inout) :: WA1(:)
      real, intent(inout) :: WA2(:)
      real, intent(inout) :: WA3(:)

      REAL, parameter :: SQRT2=1.414213562373095
      integer         :: I,K,IDP2,IC
      real            :: TR1,TR2,TR3,TR4,TI1,TI2,TI3,TI4
      real            :: CI2,CI3,CI4,CR2,CR3,CR4
      DO 101 K=1,L1
         TR1 = CC(1,1,K)-CC(IDO,4,K)
         TR2 = CC(1,1,K)+CC(IDO,4,K)
         TR3 = CC(IDO,2,K)+CC(IDO,2,K)
         TR4 = CC(1,3,K)+CC(1,3,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,2) = TR1-TR4
         CH(1,K,3) = TR2-TR3
         CH(1,K,4) = TR1+TR4
  101 CONTINUE
      IF ( (IDO-2) .LT.0 ) THEN
          GO TO 107
      ELSE IF ( (IDO-2) .EQ.0 ) THEN
          GO TO 105
      ELSE
          GO TO 102
      ENDIF
  102 IDP2 = IDO+2
!OCL NOVREC
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TI1 = CC(I,1,K)+CC(IC,4,K)
            TI2 = CC(I,1,K)-CC(IC,4,K)
            TI3 = CC(I,3,K)-CC(IC,2,K)
            TR4 = CC(I,3,K)+CC(IC,2,K)
            TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
            TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
            TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1-TR4
            CR4 = TR1+TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
            CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
            CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
            CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
            CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
            CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         TI1 = CC(1,2,K)+CC(1,4,K)
         TI2 = CC(1,4,K)-CC(1,2,K)
         TR1 = CC(IDO,1,K)-CC(IDO,3,K)
         TR2 = CC(IDO,1,K)+CC(IDO,3,K)
         CH(IDO,K,1) = TR2+TR2
         CH(IDO,K,2) = SQRT2*(TR1-TI1)
         CH(IDO,K,3) = TI2+TI2
         CH(IDO,K,4) = -SQRT2*(TR1+TI1)
  106 CONTINUE
  107 RETURN
      END


      SUBROUTINE RADB5_STOCHY (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      DIMENSION  CC(IDO,5,L1), CH(IDO,L1,5), WA1(*), WA2(*), WA3(*), WA4(*)
      REAL, parameter :: TR11=0.309016994374947
      REAL, parameter :: TI11= 0.951056516295154
      REAL, parameter :: TR12=-0.809016994374947
       REAL, parameter :: TI12=0.587785252292473
      DO 101 K=1,L1
         TI5 = CC(1,3,K)+CC(1,3,K)
         TI4 = CC(1,5,K)+CC(1,5,K)
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         TR3 = CC(IDO,4,K)+CC(IDO,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI5 = TI11*TI5+TI12*TI4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(1,K,5) = CR2+CI5
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TI5 = CC(I,3,K)+CC(IC,2,K)
            TI2 = CC(I,3,K)-CC(IC,2,K)
            TI4 = CC(I,5,K)+CC(IC,4,K)
            TI3 = CC(I,5,K)-CC(IC,4,K)
            TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
            TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
  102    CONTINUE
  103 CONTINUE
      RETURN
      END


      SUBROUTINE RADBG_STOCHY (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      DIMENSION CH(IDO,L1,IP), CC(IDO,IP,L1), C1(IDO,L1,IP), C2(IDL1,IP), &
                CH2(IDL1,IP) , WA(*)
      REAL, parameter :: TPI=6.28318530717959
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO K=1,L1
         DO I=1,IDO
            CH(I,K,1) = CC(I,1,K)
         ENDDO
      ENDDO
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
!OCL NOVREC
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            CH(1,K,J) = CC(IDO,J2-2,K)+CC(IDO,J2-2,K)
            CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
!OCL NOVREC
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
            DO 109 I=3,IDO,2
               IC = IDP2-I
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1.
      AI1 = 0.
!OCL NOVREC
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)
            C2(IK,LC) = AI1*CH2(IK,IP)
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
!OCL NOVREC
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
!OCL NOVREC
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  121    CONTINUE
  122 CONTINUE
!OCL NOVREC
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)
            CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
!OCL NOVREC
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
            DO 125 I=3,IDO,2
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            C1(1,K,J) = CH(1,K,J)
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
!OCL NOVREC
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 RETURN
      END


!>@brief The subroutine 'dozeuv_stochy' caculates odd u and even v winds harmonics from the odd harmonics
! of divergence and even harmonics of vorticty
!>@details This code is taken from the legacy spectral GFS
      subroutine dozeuv_stochy(dod,zev,uod,vev,epsedn,epsodn, snnp1ev,snnp1od,ls_node)


      implicit none
      real(kind_dbl_prec), intent(in)  :: dod(len_trio_ls,2)
      real(kind_dbl_prec), intent(in)  :: zev(len_trie_ls,2)
      real(kind_dbl_prec), intent(out) :: uod(len_trio_ls,2)
      real(kind_dbl_prec), intent(out) :: vev(len_trie_ls,2)
      real(kind_dbl_prec), intent(in)  :: epsedn(len_trie_ls)
      real(kind_dbl_prec), intent(in)  :: epsodn(len_trio_ls)
      real(kind_dbl_prec), intent(in)  :: snnp1ev(len_trie_ls)
      real(kind_dbl_prec), intent(in)  :: snnp1od(len_trio_ls)
      integer,             intent(in)  :: ls_node(ls_dim,3)
!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
! locaals
      integer              l,locl,n
      integer              indev,indev1,indev2
      integer              indod,indod1,indod2
      integer              inddif
      real(kind_dbl_prec) rl
      real(kind_dbl_prec) cons0     !constant
      integer              indlsev,jbasev
      integer              indlsod,jbasod
      real(kind_evod)  rerth

      include 'function2'


!......................................................................
      cons0 = 0.d0     !constant
      rerth  =6.3712e+6      ! radius of earth (m)

      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         vev(indlsev(l,l),1) = cons0     !constant
         vev(indlsev(l,l),2) = cons0     !constant
      enddo
!......................................................................
      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,L)
         else
            indev2 = indlsev(jcap  ,L)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
         do indev = indev1 , indev2
            uod(indev-inddif,1) = -epsodn(indev-inddif) * zev(indev,1)
            uod(indev-inddif,2) = -epsodn(indev-inddif) * zev(indev,2)
         enddo
      enddo
!......................................................................
      do locl=1,ls_max_node
          l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L) + 1
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,L)
         else
            indev2 = indlsev(jcap  ,L)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
         do indev = indev1 , indev2
            vev(indev,1) = epsedn(indev) * dod(indev-inddif,1)
            vev(indev,2) = epsedn(indev) * dod(indev-inddif,2)
         enddo
      enddo

!......................................................................
      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasod=ls_node(locl,3)
         indod1 = indlsod(L+1,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indod2 = indlsod(jcap  ,L)
         else
            indod2 = indlsod(jcap+1,L) - 1
         endif
         if ( l .ge. 1 ) then
              rl = l
            do indod = indod1 , indod2
!              u(l,n)=-i*l*d(l,n)/(n*(n+1))
               uod(indod,1) = uod(indod,1) + rl * dod(indod,2) / snnp1od(indod)
               uod(indod,2) = uod(indod,2) - rl * dod(indod,1) / snnp1od(indod)
            enddo
         endif
      enddo
!......................................................................
      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         indev1 = indlsev(L,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,L)
         else
            indev2 = indlsev(jcap  ,L)
         endif
         if ( l .ge. 1 ) then
              rl = l
            do indev = indev1 , indev2
!              u(l,n)=-i*l*d(l,n)/(n*(n+1))
               vev(indev,1) = vev(indev,1) + rl * zev(indev,2) / snnp1ev(indev)
               vev(indev,2) = vev(indev,2) - rl * zev(indev,1) / snnp1ev(indev)
            enddo
         endif
      enddo
!......................................................................
      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L) + 1
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,L)
         else
            indev2 = indlsev(jcap  ,L)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
         do indev = indev1 , indev2
            uod(indev-inddif,1) = uod(indev-inddif,1) + epsedn(indev) * zev(indev,1)
            uod(indev-inddif,2) = uod(indev-inddif,2) + epsedn(indev) * zev(indev,2)
         enddo
      enddo
!......................................................................
      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,L) - 1
         else
            indev2 = indlsev(jcap  ,L) - 1
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
         do indev = indev1 , indev2
             vev(indev,1) = vev(indev,1) - epsodn(indev-inddif) * dod(indev-inddif,1)
             vev(indev,2) = vev(indev,2) - epsodn(indev-inddif) * dod(indev-inddif,2)
         enddo
      enddo
!......................................................................
      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L)
         indod1 = indlsod(L+1,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,L)
            indod2 = indlsod(jcap  ,L)
         else
            indev2 = indlsev(jcap  ,L)
            indod2 = indlsod(jcap+1,L)
         endif
         do indod = indod1 , indod2
            uod(indod,1) = uod(indod,1) * rerth
            uod(indod,2) = uod(indod,2) * rerth
         enddo
         do indev = indev1 , indev2
            vev(indev,1) = vev(indev,1) * rerth
            vev(indev,2) = vev(indev,2) * rerth
         enddo
      enddo
      return
      end

!>@brief The subroutine 'dezouv_stochy' caculates even u and odd v winds harmonics from the even harmonics
! of divergence and odd harmonics of vorticty
!>@details This code is taken from the legacy spectral GFS
      subroutine dezouv_stochy(dev,zod,uev,vod,epsedn,epsodn,snnp1ev,snnp1od,ls_node)


      implicit none

      real(kind_dbl_prec), intent(in)  :: dev(len_trie_ls,2)
      real(kind_dbl_prec), intent(in)  :: zod(len_trio_ls,2)
      real(kind_dbl_prec), intent(out) :: uev(len_trie_ls,2)
      real(kind_dbl_prec), intent(out) :: vod(len_trio_ls,2)
      real(kind_dbl_prec), intent(in)  :: epsedn(len_trie_ls)
      real(kind_dbl_prec), intent(in)  :: epsodn(len_trio_ls)
      real(kind_dbl_prec), intent(in)  :: snnp1ev(len_trie_ls)
      real(kind_dbl_prec), intent(in)  :: snnp1od(len_trio_ls)
      integer,             intent(in)  :: ls_node(ls_dim,3)

!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod

      integer              l,locl,n
      integer              indev,indev1,indev2
      integer              indod,indod1,indod2
      integer              inddif

      real(kind_dbl_prec) rl
      real(kind_dbl_prec) cons0     !constant

      integer              indlsev,jbasev
      integer              indlsod,jbasod
      real(kind_evod)  rerth

      include 'function2'
!......................................................................
      cons0 = 0.d0     !constant
      rerth  =6.3712e+6      ! radius of earth (m)

      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         uev(indlsev(l,l),1) = cons0     !constant
         uev(indlsev(l,l),2) = cons0     !constant
      enddo

!......................................................................

      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L) + 1
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,L)
         else
            indev2 = indlsev(jcap  ,L)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
         do indev = indev1 , indev2
            uev(indev,1) = -epsedn(indev) * zod(indev-inddif,1)
            uev(indev,2) = -epsedn(indev) * zod(indev-inddif,2)
         enddo
      enddo

      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,L)
         else
            indev2 = indlsev(jcap  ,L)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1

         do indev = indev1 , indev2
            vod(indev-inddif,1) = epsodn(indev-inddif) * dev(indev,1)
            vod(indev-inddif,2) = epsodn(indev-inddif) * dev(indev,2)
         enddo
      enddo

!......................................................................

      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         indev1 = indlsev(L,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,L)
         else
            indev2 = indlsev(jcap  ,L)
         endif
         if ( l .ge. 1 ) then
              rl = l
            do indev = indev1 , indev2
               uev(indev,1) = uev(indev,1) + rl * dev(indev,2) / snnp1ev(indev)
               uev(indev,2) = uev(indev,2) - rl * dev(indev,1) / snnp1ev(indev)
            enddo
         endif
      enddo

!......................................................................

      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasod=ls_node(locl,3)
         indod1 = indlsod(L+1,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indod2 = indlsod(jcap  ,L)
         else
            indod2 = indlsod(jcap+1,L) - 1
         endif
         if ( l .ge. 1 ) then
              rl = l
            do indod = indod1 , indod2
               vod(indod,1) = vod(indod,1) + rl * zod(indod,2) / snnp1od(indod)
               vod(indod,2) = vod(indod,2) - rl * zod(indod,1) / snnp1od(indod)
            enddo
         endif
      enddo

!......................................................................

      do locl=1,ls_max_node
         l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,L) - 1
         else
            indev2 = indlsev(jcap  ,L) - 1
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1

         do indev = indev1 , indev2
            uev(indev,1) = uev(indev, 1) + epsodn(indev-inddif) * zod(indev-inddif,1)
            uev(indev,2) = uev(indev, 2) + epsodn(indev-inddif) * zod(indev-inddif,2)
         enddo
      enddo
!......................................................................
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L) + 1
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,L)
         else
            indev2 = indlsev(jcap  ,L)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
         do indev = indev1 , indev2
                 vod(indev-inddif,1) = vod(indev-inddif,1) - epsedn(indev) * dev(indev, 1)
                 vod(indev-inddif,2) = vod(indev-inddif,2) - epsedn(indev) * dev(indev, 2)
         enddo
      enddo
!......................................................................
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L)
         indod1 = indlsod(L+1,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,L)
            indod2 = indlsod(jcap  ,L)
         else
            indev2 = indlsev(jcap  ,L)
            indod2 = indlsod(jcap+1,L)
         endif
         do indev = indev1 , indev2
            uev(indev,1) = uev(indev,1) * rerth
            uev(indev,2) = uev(indev,2) * rerth
         enddo

         do indod = indod1 , indod2
            vod(indod,1) = vod(indod,1) * rerth
            vod(indod,2) = vod(indod,2) * rerth
         enddo
      enddo

      return
      end

   !  interpolation from lat/lon or gaussian grid to other lat/lon grid
   !
!>@brief The subroutine 'stochy_la2ga' intepolates from the global gaussian grid
!! to the cubed sphere points
!>@details This code is taken from the legacy spectral GFS
   subroutine stochy_la2ga(regin,imxin,jmxin,rinlon,rinlat,rlon,rlat, &
                           gauout,len,outlat, outlon)
      implicit none
      ! interface variables
      real (kind=kind_io8), intent(in)  :: regin(imxin,jmxin)
      integer,              intent(in)  :: imxin
      integer,              intent(in)  :: jmxin
      real (kind=kind_io8), intent(in)  :: rinlon(imxin)
      real (kind=kind_io8), intent(in)  :: rinlat(jmxin)
      real (kind=kind_io8), intent(in)  :: rlon
      real (kind=kind_io8), intent(in)  :: rlat
      real (kind=kind_io8), intent(out) :: gauout(len)
      integer,              intent(in)  :: len
      real (kind=kind_io8), intent(in)  :: outlat(len)
      real (kind=kind_io8), intent(in)  :: outlon(len)
      ! local variables
      real (kind=kind_io8) :: sum2,sum1,sum3,sum4
      real (kind=kind_io8) :: wsum,wsumiv,sums,sumn,wi2j2,x,y,wi1j1
      real (kind=kind_io8) :: wi1j2,wi2j1,aphi,rnume,alamd,denom
      integer              :: i,j,jq,jx
      integer              :: j1,j2,ii,i1,i2
      integer              :: iindx1(len)
      integer              :: iindx2(len)
      integer              :: jindx1(len)
      integer              :: jindx2(len)
      real(kind=kind_io8)  :: ddx(len)
      real(kind=kind_io8)  :: ddy(len)
      real(kind=kind_io8)  :: wrk(len)
!
!
!       find i-index for interpolation
        do i=1,len
          alamd = outlon(i)
          if (alamd .lt. rlon)   alamd = alamd + 360.0
          if (alamd .gt. 360.0+rlon) alamd = alamd - 360.0
          wrk(i)    = alamd
          iindx1(i) = imxin
        enddo
        do i=1,len
          do ii=1,imxin
            if(wrk(i) .ge. rinlon(ii)) iindx1(i) = ii
          enddo
        enddo
        do i=1,len
          i1 = iindx1(i)
          if (i1 .lt. 1) i1 = imxin
          i2 = i1 + 1
          if (i2 .gt. imxin) i2 = 1
          iindx1(i) = i1
          iindx2(i) = i2
          denom     = rinlon(i2) - rinlon(i1)
          if(denom.lt.0.) denom = denom + 360.
          rnume = wrk(i) - rinlon(i1)
          if(rnume.lt.0.) rnume = rnume + 360.
          ddx(i) = rnume / denom
        enddo
!
!  find j-index for interplation
!
        if(rlat.gt.0.) then
          do j=1,len
            jindx1(j)=0
          enddo
          do jx=1,jmxin
            do j=1,len
              if(outlat(j).le.rinlat(jx)) jindx1(j) = jx
            enddo
          enddo
          do j=1,len
            jq = jindx1(j)
            aphi=outlat(j)
            if(jq.ge.1 .and. jq .lt. jmxin) then
              j2=jq+1
              j1=jq
             ddy(j)=(aphi-rinlat(j1))/(rinlat(j2)-rinlat(j1))
            elseif (jq .eq. 0) then
              j2=1
              j1=1
              if(abs(90.-rinlat(j1)).gt.0.001) then
                ddy(j)=(aphi-rinlat(j1))/(90.-rinlat(j1))
              else
                ddy(j)=0.0
              endif
            else
              j2=jmxin
              j1=jmxin
              if(abs(-90.-rinlat(j1)).gt.0.001) then
                ddy(j)=(aphi-rinlat(j1))/(-90.-rinlat(j1))
              else
                ddy(j)=0.0
              endif
            endif
            jindx1(j)=j1
            jindx2(j)=j2
          enddo
        else
          do j=1,len
            jindx1(j) = jmxin+1
          enddo
          do jx=jmxin,1,-1
            do j=1,len
              if(outlat(j).le.rinlat(jx)) jindx1(j) = jx
            enddo
          enddo
          do j=1,len
            jq = jindx1(j)
            aphi=outlat(j)
            if(jq.gt.1 .and. jq .le. jmxin) then
              j2=jq
              j1=jq-1
              ddy(j)=(aphi-rinlat(j1))/(rinlat(j2)-rinlat(j1))
            elseif (jq .eq. 1) then
              j2=1
              j1=1
              if(abs(-90.-rinlat(j1)).gt.0.001) then
                 ddy(j)=(aphi-rinlat(j1))/(-90.-rinlat(j1))
              else
                 ddy(j)=0.0
              endif
            else
              j2=jmxin
              j1=jmxin
              if(abs(90.-rinlat(j1)).gt.0.001) then
                 ddy(j)=(aphi-rinlat(j1))/(90.-rinlat(j1))
              else
                 ddy(j)=0.0
              endif
            endif
            jindx1(j)=j1
            jindx2(j)=j2
          enddo
        endif
!
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
        do i=1,imxin
          sum1 = sum1 + regin(i,1)
          sum2 = sum2 + regin(i,jmxin)
        enddo
        sum1 = sum1 / imxin
        sum2 = sum2 / imxin
        sum3 = sum1
        sum4 = sum2
!
!  quasi-bilinear interpolation
!
        do i=1,len
          y  = ddy(i)
          j1 = jindx1(i)
          j2 = jindx2(i)
          x  = ddx(i)
          i1 = iindx1(i)
          i2 = iindx2(i)
!
          wi1j1 = (1.-x) * (1.-y)
          wi2j1 =     x  *( 1.-y)
          wi1j2 = (1.-x) *      y
          wi2j2 =     x  *      y
!
          wsum   = wi1j1 + wi2j1 + wi1j2 + wi2j2
          wrk(i) = wsum
          if(wsum.ne.0.) then
            wsumiv = 1./wsum
            if(j1.ne.j2) then
              gauout(i) = (wi1j1*regin(i1,j1) + wi2j1*regin(i2,j1) + &
                           wi1j2*regin(i1,j2) + wi2j2*regin(i2,j2))  &
                        *wsumiv
            else
              if (rlat .gt. 0.0) then
                sumn = sum3
                sums = sum4
                if( j1 .eq. 1) then
                  gauout(i) = (wi1j1*sumn        +wi2j1*sumn        + &
                               wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2)) &
                            * wsumiv
                elseif (j1 .eq. jmxin) then
                  gauout(i) = (wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+ &
                               wi1j2*sums        +wi2j2*sums        ) &
                            * wsumiv
                endif
              else
                sums = sum3
                sumn = sum4
                if( j1 .eq. 1) then
                  gauout(i) = (wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+ &
                               wi1j2*sums        +wi2j2*sums        ) &
                            * wsumiv
                elseif (j1 .eq. jmxin) then
                  gauout(i) = (wi1j1*sumn        +wi2j1*sumn        + &
                               wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2)) &
                            * wsumiv
                endif
              endif
            endif            ! if j1 .ne. j2
          endif
        enddo
        do i=1,len
          j1 = jindx1(i)
          j2 = jindx2(i)
          i1 = iindx1(i)
          i2 = iindx2(i)
          if(wrk(i) .eq. 0.0) then
            write(6,*) ' la2ga: error'
            call sleep(2)
            stop
          endif
        enddo
      return
!
   end subroutine stochy_la2ga

!>@brief The subroutine 'initialize_spectral' initializes the
!gridded component of the stochastic physics pattern
!>@details This code is taken from the legacy spectral GFS
      subroutine initialize_spectral(gis_stochy)

! this subroutine set up the internal state variables,
! allocate internal state arrays for initializing the gfs system.
!----------------------------------------------------------------
!
      implicit none
!
!      type(stochy_internal_state), pointer, intent(inout) :: gis_stochy
      type(stochy_internal_state), intent(inout) :: gis_stochy
      integer           :: i, l, locl

!-------------------------------------------------------------------

! set up gfs internal state dimension and values for dynamics etc
!-------------------------------------------------------------------
      gis_stochy%lon_dim_a = lon_s + 2
      jcap=ntrunc
      latg   = lat_s
      latg2  = latg/2
      lonf   = lon_s

      allocate(lat1s_a(0:jcap))
      allocate(lon_dims_a(latg))

      allocate(wgt_a(latg2))
      allocate(rcs2_a(latg2))

      ls_dim = (jcap)/gis_stochy%nodes+1
!!
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      allocate (      gis_stochy%ls_node (ls_dim,3) )
      allocate (      gis_stochy%ls_nodes(ls_dim,gis_stochy%nodes) )
      allocate (  gis_stochy%max_ls_nodes(gis_stochy%nodes) )
!
      allocate (  gis_stochy%lats_nodes_a(gis_stochy%nodes) )
      allocate ( gis_stochy%global_lats_a(latg) )
!



!---------------------------------------------------
!
      call get_ls_node_stochy( gis_stochy%mype, gis_stochy%ls_node(:,1), ls_max_node, gis_stochy%nodes)
!
!
      len_trie_ls = 0
      len_trio_ls = 0
      do locl=1,ls_max_node
         gis_stochy%ls_node(locl,2) = len_trie_ls
         gis_stochy%ls_node(locl,3) = len_trio_ls
         l = gis_stochy%ls_node(locl,1)
         len_trie_ls = len_trie_ls+(jcap+3-l)/2
         len_trio_ls = len_trio_ls+(jcap+2-l)/2
      enddo
!
      allocate ( gis_stochy%epse  (len_trie_ls) )
      allocate ( gis_stochy%epso  (len_trio_ls) )
      allocate ( gis_stochy%epsedn(len_trie_ls) )
      allocate ( gis_stochy%epsodn(len_trio_ls) )
      allocate ( gis_stochy%kenorm_e(len_trie_ls) )
      allocate ( gis_stochy%kenorm_o(len_trio_ls) )
!
      allocate ( gis_stochy%snnp1ev(len_trie_ls) )
      allocate ( gis_stochy%snnp1od(len_trio_ls) )
!
      allocate ( gis_stochy%plnev_a(len_trie_ls,latg2) )
      allocate ( gis_stochy%plnod_a(len_trio_ls,latg2) )
      allocate ( gis_stochy%plnew_a(len_trie_ls,latg2) )
      allocate ( gis_stochy%plnow_a(len_trio_ls,latg2) )

      allocate(colrad_a(latg2))
      allocate(sinlat_a(latg))
      allocate(coslat_a(latg))
!!
      call getcon_spectral(gis_stochy)
!
      gis_stochy%lats_node_a     = gis_stochy%lats_nodes_a(gis_stochy%mype+1)


      allocate ( gis_stochy%trie_ls (len_trie_ls,2,gis_stochy%lotls) )
      allocate ( gis_stochy%trio_ls (len_trio_ls,2,gis_stochy%lotls) )


      end subroutine initialize_spectral


!>@brief The subroutine 'get_ls_node_stochy' calculates the decomposition of the spherical harmonics based on the processor layout
      subroutine get_ls_node_stochy(me_fake,ls_node,ls_max_node_fake, nodes)
!>@details This code is taken from the legacy spectral GFS
!
      implicit none
!
      integer   me_fake, ls_max_node_fake, nodes
      integer   ls_node(ls_dim)

      integer   ijk, jptls, l, node, nodesio, jcap1
!
      nodesio = nodes

      ls_node = -1
      jcap1=jcap+1
!
      jptls =  0
      l = 0
!.............................................
      do ijk=1,jcap1
!
         do node=1,nodesio
            if (node == me_fake+1) then
               jptls = jptls + 1
               ls_node(jptls) = l
            endif
            l = l + 1
            if (l > jcap) go to 200
         enddo
!
         do node=nodesio,1,-1
            if (node == me_fake+1) then
               jptls = jptls + 1
               ls_node(jptls) = l
            endif
            l = l + 1
            if (l > jcap) go to 200
         enddo
!
      enddo
!.............................................
!
  200 continue
!
!.............................................
!
      ls_max_node_fake = 0
      do ijk=1,ls_dim
         if(ls_node(ijk) >= 0) then
            ls_max_node_fake = ijk
          endif
      enddo
!
      return
      end

!>@brief The subroutine 'getcon_spectral' gets various constants for the spectral and related gaussian grid
!! and caluated the assoicate legendre polynomials
!>@details This code is taken from the legacy spectral GFS
      subroutine getcon_spectral ( gis_stochy)

      implicit none
!
      integer              i,j,l,lat,n
      integer              ls_node(ls_dim,3)
!
!     ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!     ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!     ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!
      type(stochy_internal_state), intent(inout) :: gis_stochy
!
      integer       locl,node, indev, indod, indlsev,jbasev,indlsod,jbasod
!
      integer gl_lats_index
      integer global_time_sort_index_a(latg)
!
      include 'function2'
!
      real(kind=kind_dbl_prec), parameter :: cons0 = 0.d0, cons0p5  = 0.5d0,&
                                         cons1 = 1.d0, cons0p92 = 0.92d0
!
      gl_lats_index = 0
      gis_stochy%global_lats_a = -1
      global_time_sort_index_a=lonf

      do node=1,gis_stochy%nodes
          call get_lats_node_a_stochy( node-1, gis_stochy%global_lats_a,gis_stochy%lats_nodes_a(node),&
                               gl_lats_index,global_time_sort_index_a, gis_stochy%nodes)
      enddo
      call setlats_a_stochy(gis_stochy)

      do node=1,gis_stochy%nodes
         call get_ls_node_stochy( node-1, gis_stochy%ls_nodes(1,node),gis_stochy%max_ls_nodes(node), gis_stochy%nodes )
      enddo
!
      gis_stochy%lats_dim_a = 0
      do node=1,gis_stochy%nodes
         gis_stochy%lats_dim_a = max(gis_stochy%lats_dim_a,gis_stochy%lats_nodes_a(node))
      enddo

      gis_stochy%ipt_lats_node_a   = 1
      if ( gis_stochy%mype > 0 ) then
        do node=1,gis_stochy%mype
         gis_stochy%ipt_lats_node_a = gis_stochy%ipt_lats_node_a + gis_stochy%lats_nodes_a(node)
        enddo
      endif

!
      call glats_stochy(latg2,colrad_a,wgt_a,rcs2_a)
      call epslon_stochy(gis_stochy)
      call pln2eo_a_stochy(gis_stochy,latg2)
      call gozrineo_a_stochy(gis_stochy,latg2)
!
!
      do locl=1,ls_max_node
         l = gis_stochy%ls_node(locl,1)
         jbasev = gis_stochy%ls_node(locl,2)
         indev  = indlsev(l,l)
         do n = l, jcap, 2
            gis_stochy%snnp1ev(indev) = n*(n+1)
            indev                     = indev+1
         end do
      end do
!
      do locl=1,ls_max_node
         l = gis_stochy%ls_node(locl,1)
         jbasod = gis_stochy%ls_node(locl,3)
         if ( l <= jcap-1 ) then
            indod = indlsod(l+1,l)
            do n = l+1, jcap, 2
               gis_stochy%snnp1od(indod) = n*(n+1)
               indod                     = indod+1
            end do
         end if
      end do
!
!
      do locl=1,ls_max_node
         l = gis_stochy%ls_node(locl,1)
         jbasev = gis_stochy%ls_node(locl,2)
         jbasod = gis_stochy%ls_node(locl,3)
         if (mod(L,2) == mod(jcap+1,2)) then ! set even (n-l) terms of top row to zero
            gis_stochy%snnp1ev(indlsev(jcap+1,l)) = cons0
         else                                ! set odd (n-l) terms of top row to zero
            gis_stochy%snnp1od(indlsod(jcap+1,l)) = cons0
         endif
      enddo
!
      do j=1,latg
        if( j <= latg2 ) then
          sinlat_a(j) =  cos(colrad_a(j))
        else
          sinlat_a(j) = -cos(colrad_a(latg+1-j))
        endif
        coslat_a(j) = sqrt(1.-sinlat_a(j)*sinlat_a(j))
      enddo
!
      do L=0,jcap
         do lat = 1, latg2
            if ( L <= min(jcap,lonf/2) ) then
               lat1s_a(L) = lat
               go to 200
            endif
         end do
  200    continue
      end do

      do j=1,gis_stochy%lats_node_a
            lon_dims_a(j) = lonfx
      enddo
      return
      end

!>@brief The subroutine 'get_lats_node_a_stochy' calculates the decomposition of the gaussian grid based on the processor layout
!>@details This code is taken from the legacy spectral GFS
      subroutine get_lats_node_a_stochy(me_fake,global_lats_a, &
                lats_nodes_a_fake,gl_lats_index,               &
                global_time_sort_index,nodes)
!
      implicit none

      integer,intent(in)  :: me_fake
      integer,intent(in)  :: nodes
      integer,intent(in)  :: lats_nodes_a_fake
      integer,intent(inout) :: gl_lats_index
      integer,intent(inout) :: global_lats_a(latg)
      integer, intent(in) :: global_time_sort_index(latg)

      integer :: ijk
      integer :: jptlats
      integer :: lat
      integer :: node

      lat = 1

!.............................................
      do ijk=1,latg
         do node=1,nodes
            if (node.eq.me_fake+1) then
               gl_lats_index=gl_lats_index+1
               global_lats_a(gl_lats_index) = global_time_sort_index(lat)
            endif
            lat = lat + 1
            if (lat .gt. latg) go to 200
         enddo

         do node=nodes,1,-1
            if (node.eq.me_fake+1) then
               gl_lats_index=gl_lats_index+1
               global_lats_a(gl_lats_index) = global_time_sort_index(lat)
            endif
            lat = lat + 1
            if (lat .gt. latg) go to 200
         enddo
      enddo
  200 continue
      return
      end

!>@brief The subroutine 'gozrineo_a_stochy' calculates the deriviates of assoicate legendre polynomials
!>@details This code is taken from the legacy spectral GFS
      subroutine gozrineo_a_stochy(gis_stochy, num_lat)

      implicit none

      type(stochy_internal_state), intent(inout) :: gis_stochy
      integer,  intent(in)                       :: num_lat

      integer                  l,lat,locl,n
      integer                  indev,indev1,indev2
      integer                  indod,indod1,indod2
      integer                  inddif

      real(kind=kind_dbl_prec) rn,rnp1,wcsa

      real(kind=kind_dbl_prec) cons0     !constant
      real(kind=kind_dbl_prec) cons2     !constant
      real  rerth

      integer                  indlsev,jbasev
      integer                  indlsod,jbasod

      include 'function2'


      cons0 = 0.d0     !constant
      cons2 = 2.d0     !constant
      rerth  =6.3712e+6      ! radius of earth (m)


      do lat=1,num_lat

         wcsa=rcs2_a(lat)/rerth

         do locl=1,ls_max_node
                 l=gis_stochy%ls_node(locl,1)
            jbasev=gis_stochy%ls_node(locl,2)
            jbasod=gis_stochy%ls_node(locl,3)
            indev1 = indlsev(L,L)
            indod1 = indlsod(L+1,L)
            if (mod(L,2).eq.mod(jcap+1,2)) then
               indev2 = indlsev(jcap+1,L)
               indod2 = indlsod(jcap  ,L)
            else
               indev2 = indlsev(jcap  ,L)
               indod2 = indlsod(jcap+1,L)
            endif
            do indev = indev1 , indev2
               gis_stochy%plnew_a(indev,lat) = gis_stochy%plnev_a(indev,lat) * wgt_a(lat)
            enddo

            do indod = indod1 , indod2
               gis_stochy%plnow_a(indod,lat) = gis_stochy%plnod_a(indod,lat) * wgt_a(lat)
            enddo
         enddo
      enddo
      return
      end

!>@brief The subroutine 'setlats_a_stochy' selects the latitude points on this task
!>@details This code is taken from the legacy spectral GFS
      subroutine setlats_a_stochy(gis_stochy)
!
      implicit none
!
      type(stochy_internal_state), intent(inout) :: gis_stochy

      integer :: nodesio,                       &
                 jcount,jpt,lat,lats_sum,node,i,ii,  &
                 ngrptg,ngrptl,ipe,irest,idp,        &
                 ngrptgh,nodesioh           
!
      integer,allocatable :: lats_hold(:,:)
!
      allocate ( lats_hold(latg,gis_stochy%nodes) )
!
      gis_stochy%lats_nodes_a = 0
      nodesio = gis_stochy%nodes
!
      ngrptg = 0
      do lat=1,latg
         do i=1,lonf
           ngrptg = ngrptg + 1
         enddo
      enddo

!
!   ngrptg contains total number of grid points.
!
!     distribution of the grid
      nodesioh = nodesio / 2

      if (nodesioh*2 /= nodesio) then
        ngrptl = 0
        ipe    = 0
        irest  = 0
        idp    = 1

        do lat=1,latg
          ngrptl = ngrptl + lonf

          if (ngrptl*nodesio <= ngrptg+irest) then
            gis_stochy%lats_nodes_a(ipe+1)  = gis_stochy%lats_nodes_a(ipe+1) + 1
            lats_hold(idp,ipe+1) = lat
            idp = idp + 1
          else
            ipe = ipe + 1
            if (ipe <= nodesio) lats_hold(1,ipe+1) = lat
            idp    = 2
            irest  = irest + ngrptg - (ngrptl-lonf)*nodesio
            ngrptl = lonf
            gis_stochy%lats_nodes_a(ipe+1) = gis_stochy%lats_nodes_a(ipe+1) + 1
          endif
        enddo
      else
        nodesioh = nodesio/2
        ngrptgh  = ngrptg/2
        ngrptl = 0
        ipe    = 0
        irest  = 0
        idp    = 1

        do lat=1,latg/2
          ngrptl = ngrptl + lonf

          if (ngrptl*nodesioh <= ngrptgh+irest .or. lat == latg/2) then
            gis_stochy%lats_nodes_a(ipe+1)  = gis_stochy%lats_nodes_a(ipe+1) + 1
            lats_hold(idp,ipe+1) = lat
            idp = idp + 1
          else
            ipe = ipe + 1
            if (ipe <= nodesioh) then
              lats_hold(1,ipe+1) = lat
            endif
            idp    = 2
            irest  = irest + ngrptgh - (ngrptl-lonf)*nodesioh
            ngrptl = lonf
            gis_stochy%lats_nodes_a(ipe+1) = gis_stochy%lats_nodes_a(ipe+1) + 1
          endif
        enddo
        do node=1, nodesioh
          ii = nodesio-node+1
          jpt = gis_stochy%lats_nodes_a(node)
          gis_stochy%lats_nodes_a(ii) = jpt
          do i=1,jpt
            lats_hold(jpt+1-i,ii) = latg+1-lats_hold(i,node)
          enddo
        enddo


      endif
!!
!!........................................................
!!
      jpt = 0
      do node=1,nodesio
        if ( gis_stochy%lats_nodes_a(node) > 0 ) then
          do jcount=1,gis_stochy%lats_nodes_a(node)
            gis_stochy%global_lats_a(jpt+jcount) = lats_hold(jcount,node)
          enddo
        endif
        jpt = jpt + gis_stochy%lats_nodes_a(node)
      enddo

      deallocate (lats_hold)

      return
      end

!>@brief The subroutine 'glats_stochy' calculate the latitudes for the gaussian grid 
!>@details This code is taken from the legacy spectral GFS
      subroutine glats_stochy(lgghaf,colrad,wgt,rcs2)
!
! Jan 2013   Henry Juang  increase precision by kind_qdt_prec=16
!                         to help wgt (Gaussian weighting)
      implicit none
      integer                  iter,k,k1,l2,lgghaf
!
! increase precision for more significant digit to help wgt
      real(kind=kind_qdt_prec) drad,dradz,p1,p2,phi,pi,rad,rc
      real(kind=kind_qdt_prec) rl2,scale,si,sn,w,x
      real(kind=kind_dbl_prec), dimension(lgghaf) ::  colrad, wgt, rcs2
!
      real(kind=kind_dbl_prec), parameter :: cons0 = 0.d0, cons1 = 1.d0, &
                                             cons2 = 2.d0, cons4 = 4.d0, &
                                             cons180 = 180.d0, &
                                             cons0p25 = 0.25d0
#ifdef NO_QUAD_PRECISION
      real(kind=kind_qdt_prec), parameter :: eps = 1.d-12
#else
      real(kind=kind_qdt_prec), parameter :: eps = 1.d-20
#endif

!
! for better accuracy to select smaller number
!     eps = 1.d-12
!     eps = 1.d-20
!
      si    = cons1
      l2    = 2*lgghaf
      rl2   = l2
      scale = cons2/(rl2*rl2)
      k1    = l2-1
      pi    = atan(si)*cons4

!  for better accuracy to start iteration
      dradz = pi / float(lgghaf) / 200.0
      rad   = cons0
      do k=1,lgghaf
        iter = 0
        drad = dradz
1       call poly(l2,rad,p2)
2       p1 = p2
        iter = iter + 1
        rad = rad + drad
        call poly(l2,rad,p2)
        if(sign(si,p1) == sign(si,p2)) go to 2
        if(drad < eps)go to 3
        rad  = rad-drad
        drad = drad * cons0p25
        go to 1
3       continue
        colrad(k) = rad
        phi = rad * cons180 / pi
        call poly(k1,rad,p1)
        x        = cos(rad)
        w        =  scale * (cons1 - x*x)/ (p1*p1)
        wgt(k)   = w
        sn       = sin(rad)
        w        = w/(sn*sn)
        rc       = cons1/(sn*sn)
        rcs2(k)  = rc
        call poly(l2,rad,p1)
      enddo
!
      return
      end

!>@brief The subroutine 'poly' does something with latitudes
!>@details This code is taken from the legacy spectral GFS
      subroutine poly(n,rad,p)
!
      implicit none
      integer                  i,n
!
! increase precision for more significant digit to help wgt
      real(kind=kind_qdt_prec) floati,g,p,rad,x,y1,y2,y3
!
      real(kind=kind_dbl_prec), parameter ::  cons1 = 1.d0
!
      x  = cos(rad)
      y1 = cons1
      y2 = x
      do i=2,n
        g = x*y2
        floati = i
        y3 = g - y1 + g - (g-y1)/floati
        y1 = y2
        y2 = y3
      enddo
      p = y3
      return
      end

!>@brief The subroutine 'pln2eo_a_stochy' calculates the assoicated legendre polynomials
!>@details This code is taken from the legacy spectral GFS
      subroutine pln2eo_a_stochy(gis_stochy,num_lat)
!
! use x-number method to archieve accuracy due to recursive to avoid
! underflow and overflow if necessary by henry juang 2012 july
!
      implicit none
!
! define x number constant for real8 start
      type(stochy_internal_state), intent(inout) :: gis_stochy
      integer, intent(in)                        :: num_lat
      integer,  parameter :: in_f = 960 , in_h = in_f/2
      real(kind=kind_dbl_prec), parameter :: bb_f = 2.d0 ** ( in_f )
      real(kind=kind_dbl_prec), parameter :: bs_f = 2.d0 ** (-in_f )
      real(kind=kind_dbl_prec), parameter :: bb_h = 2.d0 ** ( in_h )
      real(kind=kind_dbl_prec), parameter :: bs_h = 2.d0 ** (-in_h )
! define x number constant end

!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
      integer                  l,lat,locl,max_l,n
      integer                  indev
      integer                  indod
! need index for alp to be x-number
      integer                  id, ialp1, ialp2, ialp3, iprod
      integer                  ialp10(0:jcap)
      real(kind=kind_dbl_prec) aa, bb, w

      real(kind=kind_dbl_prec) alp1,alp2,alp3
      real(kind=kind_dbl_prec) cos2,fl,prod,sinlat,coslat
      real(kind=kind_dbl_prec) alp10(0:jcap)
      real(kind=kind_dbl_prec) cons0,cons0p5,cons1,cons2,cons3    !constant
      integer                  indlsev,jbasev
      integer                  indlsod,jbasod

      include 'function2'

      cons0=0.0d0       !constant
      cons0p5=0.5d0     !constant
      cons1=1.0d0       !constant
      cons2=2.0d0       !constant
      cons3=3.0d0       !constant

      max_l=-1
      do locl=1,ls_max_node
         max_l = max ( max_l, gis_stochy%ls_node(locl,1) )
      enddo

      do lat=1,num_lat

         sinlat = cos(colrad_a(lat))
         cos2=cons1-sinlat*sinlat           !constant
         coslat = sqrt(cos2)

! use x number for alp10
         alp10(0) = sqrt(0.5)
         ialp10(0) = 0

         do l=1,max_l
            fl = l
            prod=coslat*sqrt(cons1+cons1/(cons2*fl))
            iprod=0
            w = abs(prod)
            if( w.ge.bb_h ) then
              prod = prod * bs_f
              iprod = iprod + 1
            elseif( w.lt.bs_h ) then
              prod = prod * bb_f
              iprod = iprod - 1
            endif
            alp10(l)=alp10(l-1)*prod
            ialp10(l)=ialp10(l-1)+iprod
            w = abs(alp10(l))
            if( w.ge.bb_h ) then
              alp10(l) = alp10(l) * bs_f
              ialp10(l) = ialp10(l) + 1
            elseif( w.lt.bs_h ) then
              alp10(l) = alp10(l) * bb_f
              ialp10(l) = ialp10(l) - 1
            endif
         enddo

         do locl=1,ls_max_node
                 l=gis_stochy%ls_node(locl,1)
            jbasev=gis_stochy%ls_node(locl,2)
            jbasod=gis_stochy%ls_node(locl,3)
            n=l
            fl=l
! get m=normalized x number for alp1 start
            alp1=alp10(l)
            ialp1=ialp10(l)

            indev=indlsev(n  ,l)
            indod=indlsod(n+1,l)
! x2f start
            if( ialp1.eq.0 ) then
              gis_stochy%plnev_a(indev     ,lat)=alp1
            elseif( ialp1.eq.-1 ) then
              gis_stochy%plnev_a(indev     ,lat)=alp1 * bs_f
            elseif( ialp1.lt.-1 ) then
              gis_stochy%plnev_a(indev     ,lat)=0.0
            else
              gis_stochy%plnev_a(indev     ,lat)=alp1 * bb_f
            endif
! x2f end

! xltime    alp2=sqrt(cons2*fl+cons3)*sinlat*alp1     !constant
! xltime start
            prod=sqrt(cons2*fl+cons3)*sinlat
            iprod=0
            w = abs(prod)
            if( w.ge.bb_h ) then
              prod = prod * bs_f
              iprod = iprod + 1
            elseif( w.lt.bs_h ) then
              prod = prod * bb_f
              iprod = iprod - 1
            endif
            alp2=alp1*prod
            ialp2 = ialp1 + iprod
! xltime end
! norm alp2 start
            w = abs(alp2)
            if( w.ge.bb_h ) then
              alp2 = alp2*bs_f
              ialp2 = ialp2 + 1
            elseif( w.lt.bs_h ) then
              alp2 = alp2*bb_f
              ialp2 = ialp2 - 1
            endif
! norm alp2 end

! x2f start
            if( ialp2.eq.0 ) then
              gis_stochy%plnod_a(indod       ,lat)=alp2
            elseif( ialp2.eq.-1 ) then
              gis_stochy%plnod_a(indod       ,lat)=alp2 * bs_f
            elseif( ialp2.lt.-1 ) then
              gis_stochy%plnod_a(indod       ,lat)=0.0
            else
              gis_stochy%plnod_a(indod       ,lat)=alp2 * bb_f
            endif
! x2f end

            do n=l+2,jcap+1
               if(mod(n+l,2).eq.0) then
                  indev=indev+1
! xlsum2 start
                  aa = sinlat / gis_stochy%epse(indev)
                  bb = gis_stochy%epso(indod) / gis_stochy%epse(indev)
                  id = ialp2 - ialp1
                  if( id.eq.0 ) then
                    alp3 = aa*alp2 - bb*alp1
                    ialp3 = ialp1
                  elseif( id.eq.1 ) then
                    alp3 = aa*alp2 - bb*alp1*bs_f
                    ialp3 = ialp2
                  elseif( id.eq.-1 ) then
                    alp3 = aa*alp2*bs_f - bb*alp1
                    ialp3 = ialp1
                  elseif( id.gt.1 ) then
                    alp3 = aa*alp2
                    ialp3 = ialp2
                  else
                    alp3 = - bb*alp1
                    ialp3 = ialp1
                  endif
! xlsum2 end
! xnorm alp3 start
                  w = abs(alp3)
                  if( w.ge.bb_h ) then
                    alp3 = alp3*bs_f
                    ialp3 = ialp3 + 1
                  elseif( w.lt.bs_h ) then
                    alp3 = alp3*bb_f
                    ialp3 = ialp3 - 1
                  endif
! xnorm alp3 end

! x2f alp3 start
                  if( ialp3.eq.0 ) then
                    gis_stochy%plnev_a(indev,lat)=alp3
                  elseif( ialp3.eq.-1 ) then
                    gis_stochy%plnev_a(indev,lat)=alp3 * bs_f
                  elseif( ialp3.lt.-1 ) then
                    gis_stochy%plnev_a(indev,lat)=0.0
                  else
                    gis_stochy%plnev_a(indev,lat)=alp3 * bb_f
                  endif
! x2f alp3 end

               else
                  indod=indod+1

! xlsum2 start
                  aa = sinlat / gis_stochy%epso(indod)
                  bb = gis_stochy%epse(indev) / gis_stochy%epso(indod)
                  id = ialp2 - ialp1
                  if( id.eq.0 ) then
                    alp3 = aa*alp2 - bb*alp1
                    ialp3 = ialp1
                  elseif( id.eq.1 ) then
                    alp3 = aa*alp2 - bb*alp1*bs_f
                    ialp3 = ialp2
                  elseif( id.eq.-1 ) then
                    alp3 = aa*alp2*bs_f - bb*alp1
                    ialp3 = ialp1
                  elseif( id.gt.1 ) then
                    alp3 = aa*alp2
                    ialp3 = ialp2
                  else
                    alp3 = - bb*alp1
                    ialp3 = ialp1
                  endif
! xlsum2 end
! xnorm alp3 start
                  w = abs(alp3)
                  if( w.ge.bb_h ) then
                    alp3 = alp3*bs_f
                    ialp3 = ialp3 + 1
                  elseif( w.lt.bs_h ) then
                    alp3 = alp3*bb_f
                    ialp3 = ialp3 - 1
                  endif
! xnorm alp3 end

! x2f alp3 start
                  if( ialp3.eq.0 ) then
                    gis_stochy%plnod_a(indod,lat)=alp3
                  elseif( ialp3.eq.-1 ) then
                    gis_stochy%plnod_a(indod,lat)=alp3 * bs_f
                  elseif( ialp3.lt.-1 ) then
                    gis_stochy%plnod_a(indod,lat)=0.0
                  else
                    gis_stochy%plnod_a(indod,lat)=alp3 * bb_f
                  endif
! x2f alp3 end
               endif
               alp1=alp2
               alp2=alp3
               ialp1 = ialp2
               ialp2 = ialp3
            enddo
         enddo
      enddo

      return
      end

!>@brief The subroutine 'epslon_stochy' calculate coeffients for use in spectral space
!>@details This code is taken from the legacy spectral GFS
      subroutine epslon_stochy(gis_stochy)

      implicit none

      type(stochy_internal_state), intent(inout) :: gis_stochy

      integer                  ls_node(ls_dim,3)

!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod

      integer                  l,locl,n

      integer                  indev
      integer                  indod

      real(kind_dbl_prec) f1,f2,rn,val

      real(kind_dbl_prec) cons0     !constant

      integer                  indlsev,jbasev
      integer                  indlsod,jbasod

      include 'function2'

      cons0=0.0d0     !constant
!c
!c......................................................................
!c
      do locl=1,ls_max_node
              l=gis_stochy%ls_node(locl,1)
         jbasev=gis_stochy%ls_node(locl,2)
         indev=indlsev(l,l)
         gis_stochy%epse  (indev)=cons0     !constant
         gis_stochy%epsedn(indev)=cons0     !constant
          indev=indev+1

         do n=l+2,jcap+1,2
            rn=n
            f1=n*n-l*l
            f2=4*n*n-1
            val=sqrt(f1/f2)
            gis_stochy%epse  (indev)=val
            gis_stochy%epsedn(indev)=val/rn
             indev=indev+1
         enddo
      enddo
      do locl=1,ls_max_node
              l=gis_stochy%ls_node(locl,1)
         jbasod=gis_stochy%ls_node(locl,3)
         indod=indlsod(l+1,l)

         do n=l+1,jcap+1,2
            rn=n
            f1=n*n-l*l
            f2=4*n*n-1
            val=sqrt(f1/f2)
            gis_stochy%epso  (indod)=val
            gis_stochy%epsodn(indod)=val/rn
             indod=indod+1
         enddo
      enddo

      return
      end
 end module spectral_transforms
