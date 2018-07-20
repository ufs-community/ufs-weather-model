!$$$  Module Documentation Block
!
! Module: mersenne_twister   Modern random number generator
!   Prgmmr: Iredell          Org: W/NX23     date: 2005-06-14
!
! Abstract: This module calculates random numbers using the Mersenne twister.
!   (It has been adapted to a Fortran 90 module from open source software.
!   The comments from the original software are given below in the remarks.)
!   The Mersenne twister (aka MT19937) is a state-of-the-art random number
!   generator based on Mersenne primes and originally developed in 1997 by
!   Matsumoto and Nishimura. It has a period before repeating of 2^19937-1,
!   which certainly should be good enough for geophysical purposes. :-)
!   Considering the algorithm's robustness, it runs fairly speedily.
!   (Some timing statistics are given below in the remarks.)
!   This adaptation uses the standard Fortran 90 random number interface,
!   which can generate an arbitrary number of random numbers at one time.
!   The random numbers generated are uniformly distributed between 0 and 1.
!   The module also can generate random numbers from a Gaussian distribution
!   with mean 0 and standard deviation 1, using a Numerical Recipes algorithm.
!   The module also can generate uniformly random integer indices.
!   There are also thread-safe versions of the generators in this adaptation,
!   necessitating the passing of generator states which must be kept private.
!
! Program History Log:
!   2005-06-14  Mark Iredell
!
! Usage: 
!   The module can be compiled with 4-byte reals or with 8-byte reals, but
!   4-byte integers are required. The module should be endian-independent.
!   The Fortran 90 interfaces random_seed and random_number are overloaded
!   and can be used as in the standard by adding the appropriate use statement
!     use mersenne_twister
!   In the below use cases, harvest is a real array of arbitrary size,
!   and iharvest is an integer array of arbitrary size.
!   To generate uniformly distributed random numbers between 0 and 1,
!     call random_number(harvest)
!   To generate Gaussian distributed random numbers with 0 mean and 1 sigma,
!     call random_gauss(harvest)
!   To generate uniformly distributed random integer indices between 0 and n,
!     call random_index(n,iharvest)
!   In standard "saved" mode, the random number generator can be used without
!   setting a seed. But to set a seed, only 1 non-zero integer is required, e.g.
!     call random_setseed(4357) ! set default seed
!   The full generator state can be set via the standard interface random_seed,
!   but it is recommended to use this method only to restore saved states, e.g.
!     call random_seed(size=lsave)  ! get size of generator state seed array
!     allocate isave(lsave)         ! allocate seed array
!     call random_seed(get=isave)   ! fill seed array (then maybe save to disk)
!     call random_seed(put=isave)   ! restore state (after read from disk maybe)
!   Locally kept generator states can also be saved in a seed array, e.g.
!     type(random_stat):: stat
!     call random_seed(get=isave,stat=stat)  ! fill seed array
!     call random_seed(put=isave,stat=stat)  ! restore state
!   To generate random numbers in a threaded region, the "thread-safe" mode
!   must be used where generator states of type random_state are passed, e.g.
!     type(random_stat):: stat(8)
!     do i=1,8                               ! threadable loop
!       call random_setseed(7171*i,stat(i))  ! thread-safe call
!     enddo
!     do i=1,8                               ! threadable loop
!       call random_number(harvest,stat(i))  ! thread-safe call
!     enddo
!     do i=1,8                               ! threadable loop
!       call random_gauss(harvest,stat(i))   ! thread-safe call
!     enddo
!     do i=1,8                               ! threadable loop
!       call random_index(n,iharvest,stat(i))! thread-safe call
!     enddo
!   There is also a relatively inefficient "interactive" mode available, where
!   setting seeds and generating random numbers are done in the same call.
!   There is also a functional mode available, returning one value at a time.
!   
! Public Defined Types:
!   random_stat       Generator state (private contents)
!   
! Public Subprograms:
!   random_seed         determine size or put or get state
!     size              optional integer output size of seed array
!     put               optional integer(:) input seed array
!     get               optional integer(:) output seed array
!     stat              optional type(random_stat) (thread-safe mode)
!   random_setseed      set seed (thread-safe mode)
!     inseed            integer seed input
!     stat              type(random_stat) output
!   random_setseed      set seed (saved mode)
!     inseed            integer seed input
!   random_number       get mersenne twister random numbers (thread-safe mode)
!     harvest           real(:) numbers output
!     stat              type(random_stat) input
!   random_number       get mersenne twister random numbers (saved mode)
!     harvest           real(:) numbers output
!   random_number       get mersenne twister random numbers (interactive mode)
!     harvest           real(:) numbers output
!     inseed            integer seed input
!   random_number_f     get mersenne twister random number (functional mode)
!     harvest           real number output
!   random_gauss        get gaussian random numbers (thread-safe mode)
!     harvest           real(:) numbers output
!     stat              type(random_stat) input
!   random_gauss        get gaussian random numbers (saved mode)
!     harvest           real(:) numbers output
!   random_gauss        get gaussian random numbers (interactive mode)
!     harvest           real(:) numbers output
!     inseed            integer seed input
!   random_gauss_f      get gaussian random number (functional mode)
!     harvest           real number output
!   random_index        get random indices (thread-safe mode)
!     imax              integer maximum index input
!     iharvest          integer(:) numbers output
!     stat              type(random_stat) input
!   random_index        get random indices (saved mode)
!     imax              integer maximum index input
!     iharvest          integer(:) numbers output
!   random_index        get random indices (interactive mode)
!     imax              integer maximum index input
!     iharvest          integer(:) numbers output
!     inseed            integer seed input
!   random_index_f      get random index (functional mode)
!     imax              integer maximum index input
!     iharvest          integer number output
!
! Remarks:
!   (1) Here are the comments in the original open source code:
!     A C-program for MT19937: Real number version
!     genrand() generates one pseudorandom real number (double)
!     which is uniformly distributed on [0,1]-interval, for each
!     call. sgenrand(seed) set initial values to the working area
!     of 624 words. Before genrand(), sgenrand(seed) must be
!     called once. (seed is any 32-bit integer except for 0).
!     Integer generator is obtained by modifying two lines.
!     Coded by Takuji Nishimura, considering the suggestions by
!     Topher Cooper and Marc Rieffel in July-Aug. 1997.
!     This library is free software; you can redistribute it and/or
!     modify it under the terms of the GNU Library General Public
!     License as published by the Free Software Foundation; either
!     version 2 of the License, or (at your option) any later
!     version.
!     This library is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!     See the GNU Library General Public License for more details.
!     You should have received a copy of the GNU Library General
!     Public License along with this library; if not, write to the
!     Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
!     02111-1307  USA
!     Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
!     When you use this, send an email to: matumoto@math.keio.ac.jp
!     with an appropriate reference to your work.
!     Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   (2) On a single IBM Power4 processor on the NCEP operational cluster (2005)
!     each Mersenne twister random number takes less than 30 ns, about 3 times
!     slower than the default random number generator, and each random number
!     from a Gaussian distribution takes less than 150 ns.
!     
! Attributes:
!   Language: Fortran 90
!
!$$$
      module mersenne_twister
        private
!  Public declarations
        public random_stat
        public random_seed
        public random_setseed
        public random_number
        public random_number_f
        public random_gauss
        public random_gauss_f
        public random_index
        public random_index_f
!  Parameters
        integer,parameter:: n=624
        integer,parameter:: m=397
        integer,parameter:: mata=-1727483681 ! constant vector a
        integer,parameter:: umask=-2147483648 ! most significant w-r bits
        integer,parameter:: lmask =2147483647 ! least significant r bits
        integer,parameter:: tmaskb=-1658038656 ! tempering parameter
        integer,parameter:: tmaskc=-272236544 ! tempering parameter
        integer,parameter:: mag01(0:1)=(/0,mata/)
        integer,parameter:: iseed=4357
        integer,parameter:: nrest=n+6
!  Defined types
        type random_stat
          private
          integer:: mti=n+1
          integer:: mt(0:n-1)
          integer:: iset
          real:: gset
        end type
!  Saved data
        type(random_stat),save:: sstat
!  Overloaded interfaces
        interface random_setseed
          module procedure random_setseed_s
          module procedure random_setseed_t
        end interface
        interface random_number
          module procedure random_number_i
          module procedure random_number_s
          module procedure random_number_t
        end interface
        interface random_gauss
          module procedure random_gauss_i
          module procedure random_gauss_s
          module procedure random_gauss_t
        end interface
        interface random_index
          module procedure random_index_i
          module procedure random_index_s
          module procedure random_index_t
        end interface
!  All the subprograms
      contains
!  Subprogram random_seed
!  Sets and gets state; overloads Fortran 90 standard.
        subroutine random_seed(size,put,get,stat)
          implicit none
          integer,intent(out),optional:: size
          integer,intent(in),optional:: put(nrest)
          integer,intent(out),optional:: get(nrest)
          type(random_stat),intent(inout),optional:: stat
          if(present(size)) then     ! return size of seed array
!     if(present(put).or.present(get))&
!      call errmsg('RANDOM_SEED: more than one option set - some ignored')
            size=nrest
          elseif(present(put)) then  ! restore from seed array
!     if(present(get))&
!      call errmsg('RANDOM_SEED: more than one option set - some ignored')
            if(present(stat)) then
              stat%mti=put(1)
              stat%mt=put(2:n+1)
              stat%iset=put(n+2)
              stat%gset=transfer(put(n+3:nrest),stat%gset)
              if(stat%mti.lt.0.or.stat%mti.gt.n.or.any(stat%mt.eq.0).or.
     &         stat%iset.lt.0.or.stat%iset.gt.1) then
                call random_setseed_t(iseed,stat)
!         call errmsg('RANDOM_SEED: invalid seeds put - default seeds used')
              endif
            else
              sstat%mti=put(1)
              sstat%mt=put(2:n+1)
              sstat%iset=put(n+2)
              sstat%gset=transfer(put(n+3:nrest),sstat%gset)
              if(sstat%mti.lt.0.or.sstat%mti.gt.n.or.any(sstat%mt.eq.0)
     &         .or.sstat%iset.lt.0.or.sstat%iset.gt.1) then
                call random_setseed_t(iseed,sstat)
!         call errmsg('RANDOM_SEED: invalid seeds put - default seeds used')
              endif
            endif
          elseif(present(get)) then  ! save to seed array
            if(present(stat)) then
              if(stat%mti.eq.n+1) call random_setseed_t(iseed,stat)
              get(1)=stat%mti
              get(2:n+1)=stat%mt
              get(n+2)=stat%iset
              get(n+3:nrest)=transfer(stat%gset,get,nrest-(n+3)+1)
            else
              if(sstat%mti.eq.n+1) call random_setseed_t(iseed,sstat)
              get(1)=sstat%mti
              get(2:n+1)=sstat%mt
              get(n+2)=sstat%iset
              get(n+3:nrest)=transfer(sstat%gset,get,nrest-(n+3)+1)
            endif
          else                       ! reset default seed
            if(present(stat)) then
              call random_setseed_t(iseed,stat)
            else
              call random_setseed_t(iseed,sstat)
            endif
          endif
        end subroutine
!  Subprogram random_setseed_s
!  Sets seed in saved mode.
        subroutine random_setseed_s(inseed)
          implicit none
          integer,intent(in):: inseed
          call random_setseed_t(inseed,sstat)
        end subroutine
!  Subprogram random_setseed_t
!  Sets seed in thread-safe mode.
        subroutine random_setseed_t(inseed,stat)
          implicit none
          integer,intent(in):: inseed
          type(random_stat),intent(out):: stat
          integer ii,mti
          ii=inseed
          if(ii.eq.0) ii=iseed
          stat%mti=n
          stat%mt(0)=iand(ii,-1)
          do mti=1,n-1
            stat%mt(mti)=iand(69069*stat%mt(mti-1),-1)
          enddo
          stat%iset=0
          stat%gset=0.
        end subroutine
!  Subprogram random_number_f
!  Generates random numbers in functional mode.
        function random_number_f() result(harvest)
          implicit none
          real:: harvest
          real h(1)
          if(sstat%mti.eq.n+1) call random_setseed_t(iseed,sstat)
          call random_number_t(h,sstat)
          harvest=h(1)
        end function
!  Subprogram random_number_i
!  Generates random numbers in interactive mode.
        subroutine random_number_i(harvest,inseed)
          implicit none
          real,intent(out):: harvest(:)
          integer,intent(in):: inseed
          type(random_stat) stat
          call random_setseed_t(inseed,stat)
          call random_number_t(harvest,stat)
        end subroutine
!  Subprogram random_number_s
!  Generates random numbers in saved mode; overloads Fortran 90 standard.
        subroutine random_number_s(harvest)
          implicit none
          real,intent(out):: harvest(:)
          if(sstat%mti.eq.n+1) call random_setseed_t(iseed,sstat)
          call random_number_t(harvest,sstat)
        end subroutine
!  Subprogram random_number_t
!  Generates random numbers in thread-safe mode.
        subroutine random_number_t(harvest,stat)
          implicit none
          real,intent(out):: harvest(:)
          type(random_stat),intent(inout):: stat
          integer j,kk,y
          integer tshftu,tshfts,tshftt,tshftl
          tshftu(y)=ishft(y,-11)
          tshfts(y)=ishft(y,7)
          tshftt(y)=ishft(y,15)
          tshftl(y)=ishft(y,-18)
          do j=1,size(harvest)
            if(stat%mti.ge.n) then
              do kk=0,n-m-1
                y=ior(iand(stat%mt(kk),umask),iand(stat%mt(kk+1),lmask))
                stat%mt(kk)=ieor(ieor(stat%mt(kk+m),ishft(y,-1)),
     &   mag01(iand(y,1)))
              enddo
              do kk=n-m,n-2
                y=ior(iand(stat%mt(kk),umask),iand(stat%mt(kk+1),lmask))
                stat%mt(kk)=ieor(ieor(stat%mt(kk+(m-n)),ishft(y,-1)),
     &   mag01(iand(y,1)))
              enddo
              y=ior(iand(stat%mt(n-1),umask),iand(stat%mt(0),lmask))
              stat%mt(n-1)=ieor(ieor(stat%mt(m-1),ishft(y,-1)),
     &   mag01(iand(y,1)))
              stat%mti=0
            endif
            y=stat%mt(stat%mti)
            y=ieor(y,tshftu(y))
            y=ieor(y,iand(tshfts(y),tmaskb))
            y=ieor(y,iand(tshftt(y),tmaskc))
            y=ieor(y,tshftl(y))
            if(y.lt.0) then
              harvest(j)=(real(y)+2.0**32)/(2.0**32-1.0)
            else
              harvest(j)=real(y)/(2.0**32-1.0)
            endif
            stat%mti=stat%mti+1
          enddo
        end subroutine
!  Subprogram random_gauss_f
!  Generates Gaussian random numbers in functional mode.
        function random_gauss_f() result(harvest)
          implicit none
          real:: harvest
          real h(1)
          if(sstat%mti.eq.n+1) call random_setseed_t(iseed,sstat)
          call random_gauss_t(h,sstat)
          harvest=h(1)
        end function
!  Subprogram random_gauss_i
!  Generates Gaussian random numbers in interactive mode.
        subroutine random_gauss_i(harvest,inseed)
          implicit none
          real,intent(out):: harvest(:)
          integer,intent(in):: inseed
          type(random_stat) stat
          call random_setseed_t(inseed,stat)
          call random_gauss_t(harvest,stat)
        end subroutine
!  Subprogram random_gauss_s
!  Generates Gaussian random numbers in saved mode.
        subroutine random_gauss_s(harvest)
          implicit none
          real,intent(out):: harvest(:)
          if(sstat%mti.eq.n+1) call random_setseed_t(iseed,sstat)
          call random_gauss_t(harvest,sstat)
        end subroutine
!  Subprogram random_gauss_t
!  Generates Gaussian random numbers in thread-safe mode.
        subroutine random_gauss_t(harvest,stat)
          implicit none
          real,intent(out):: harvest(:)
          type(random_stat),intent(inout):: stat
          integer mx,my,mz,j
          real r2(2),r,g1,g2
          mz=size(harvest)
          if(mz.le.0) return
          mx=0
          if(stat%iset.eq.1) then
            mx=1
            harvest(1)=stat%gset
            stat%iset=0
          endif
          my=(mz-mx)/2*2+mx
          do
            call random_number_t(harvest(mx+1:my),stat)
            do j=mx,my-2,2
              call rgauss(harvest(j+1),harvest(j+2),r,g1,g2)
              if(r.lt.1.) then
                harvest(mx+1)=g1
                harvest(mx+2)=g2
                mx=mx+2
              endif
            enddo
            if(mx.eq.my) exit
          enddo
          if(my.lt.mz) then
            do
              call random_number_t(r2,stat)
              call rgauss(r2(1),r2(2),r,g1,g2)
              if(r.lt.1.) exit
            enddo
            harvest(mz)=g1
            stat%gset=g2
            stat%iset=1
          endif
        contains
!  Numerical Recipes algorithm to generate Gaussian random numbers.
          subroutine rgauss(r1,r2,r,g1,g2)
            real,intent(in):: r1,r2
            real,intent(out):: r,g1,g2
            real v1,v2,fac
            v1=2.*r1-1.
            v2=2.*r2-1.
            r=v1**2+v2**2
            if(r.lt.1.) then
              fac=sqrt(-2.*log(r)/r)
              g1=v1*fac
              g2=v2*fac
            endif
          end subroutine
        end subroutine
!  Subprogram random_index_f
!  Generates random indices in functional mode.
        function random_index_f(imax) result(iharvest)
          implicit none
          integer,intent(in):: imax
          integer:: iharvest
          integer ih(1)
          if(sstat%mti.eq.n+1) call random_setseed_t(iseed,sstat)
          call random_index_t(imax,ih,sstat)
          iharvest=ih(1)
        end function
!  Subprogram random_index_i
!  Generates random indices in interactive mode.
        subroutine random_index_i(imax,iharvest,inseed)
          implicit none
          integer,intent(in):: imax
          integer,intent(out):: iharvest(:)
          integer,intent(in):: inseed
          type(random_stat) stat
          call random_setseed_t(inseed,stat)
          call random_index_t(imax,iharvest,stat)
        end subroutine
!  Subprogram random_index_s
!  Generates random indices in saved mode.
        subroutine random_index_s(imax,iharvest)
          implicit none
          integer,intent(in):: imax
          integer,intent(out):: iharvest(:)
          if(sstat%mti.eq.n+1) call random_setseed_t(iseed,sstat)
          call random_index_t(imax,iharvest,sstat)
        end subroutine
!  Subprogram random_index_t
!  Generates random indices in thread-safe mode.
        subroutine random_index_t(imax,iharvest,stat)
          implicit none
          integer,intent(in):: imax
          integer,intent(out):: iharvest(:)
          type(random_stat),intent(inout):: stat
          integer,parameter:: mh=n
          integer i1,i2,mz
          real h(mh)
          mz=size(iharvest)
          do i1=1,mz,mh
            i2=min((i1-1)+mh,mz)
            call random_number_t(h(:i2-(i1-1)),stat)
            iharvest(i1:i2)=max(ceiling(h(:i2-(i1-1))*imax),1)
          enddo
        end subroutine
      end module
