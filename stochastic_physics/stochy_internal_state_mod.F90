!>@brief The module 'stochy_internal_state_mod' contains the spherical harmonic definitions and arrays
!! describing the target gaussian grid
!
! !module: stochy_internal_state_mod
!                         --- internal state definition of the
!                             gridded component of the spectral random patterns
!
! !description:  define the spectral internal state used to
!                                             create the internal state.
!---------------------------------------------------------------------------
! !revision history:
!
!  Oct 11 2016     P Pegion port of gfs_dynamics_interal_state
!
! !interface:
!

      module stochy_internal_state_mod

!!uses:
!------
!      use spectral_layout_mod


      implicit none
      private

! -----------------------------------------------
      type,public::stochy_internal_state		! start type define
! -----------------------------------------------

      integer                           :: nodes
      integer                           ::  lats_node_a
      integer                           ::  mype
      integer                           :: lon_dim_a
      integer                           :: lats_dim_a
      integer                           :: ipt_lats_node_a

      integer              ,allocatable ::      ls_node    (:,:)
      integer              ,allocatable ::      ls_nodes   (:,:)
      integer              ,allocatable ::  max_ls_nodes   (:)

      integer              ,allocatable ::  lats_nodes_a   (:)
      integer              ,allocatable ::  global_lats_a  (:)
      integer              ,allocatable ::  global_lats_h  (:)
      integer                           :: xhalo,yhalo

      real,allocatable ::        epse  (:)
      real,allocatable ::        epso  (:)
      real,allocatable ::        epsedn(:)
      real,allocatable ::        epsodn(:)
      real,allocatable ::        kenorm_e(:)
      real,allocatable ::        kenorm_o(:)

      real,allocatable ::       snnp1ev(:)
      real,allocatable ::       snnp1od(:)

      real,allocatable ::       plnev_a(:,:)
      real,allocatable ::       plnod_a(:,:)
      real,allocatable ::       plnew_a(:,:)
      real,allocatable ::       plnow_a(:,:)


      real,allocatable ::       trie_ls(:,:,:)
      real,allocatable ::       trio_ls(:,:,:)

      integer              lotls
      integer              nlunit

      integer              k,l,locl,n
      integer              lan,lat
      integer              nx,ny,nz
      integer, allocatable :: len(:)
      real, allocatable :: parent_lons(:,:),parent_lats(:,:)


!
! -----------------------------------------------------
      end type stochy_internal_state		! end type define
! -----------------------------------------------------

! this state is supported by c pointer not f90 pointer, thus
! need this wrap.
!-----------------------------------------------------------
      type stochy_wrap		! begin type define
          type (stochy_internal_state), pointer :: int_state
      end type stochy_wrap	! end type define

      end module stochy_internal_state_mod
