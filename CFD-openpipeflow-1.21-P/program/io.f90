!**************************************************************************
!  IN/OUT 
!
!**************************************************************************
#include "../parallel.h"
 module io
!**************************************************************************
   use velocity
   use mpif
   use netcdf

   implicit none
   save

   character(200)        :: io_statefile
   integer               :: io_save1, io_save2 , io_save3
   integer,     private  :: io_KE, io_ID, io_dt, io_pt, io_fr , io_st
   type (coll), private  :: c1, c2, c3
   type (phys), private  :: pr, p1, p2
   type (lumesh), private ::  LNps(0:i_pH1)
 contains
 
!--------------------------------------------------------------------------
!  initialiser fn
!--------------------------------------------------------------------------
   subroutine io_precompute()
      io_statefile = 'state.cdf.in'
      io_save1 = 0
      io_save2 = 0
      io_dt    = 10
      io_KE    = 0
      io_ID    = 0
      io_pt    = 0
      io_fr    = 50
   end subroutine io_precompute 
 

!--------------------------------------------------------------------------
!  Open files written to every ... steps runtime
!--------------------------------------------------------------------------
   subroutine io_openfiles()
      character(10), save :: s = 'unknown', a = 'sequential'
      if(mpi_rnk/=0) return
      if(io_dt/=0)  open(io_dt,status=s,access=a, file='tim_step.dat')
      if(io_KE/=0)  open(io_KE,status=s,access=a, file='vel_energy.dat')
      if(io_ID/=0)  open(io_ID,status=s,access=a, file='vel_totEID.dat')
      if(io_pt/=0)  open(io_pt,status=s,access=a, file='vel_point.dat')
      if(io_fr/=0)  open(io_fr,status=s,access=a, file='vel_friction.dat')
!      s = 'old'
      a = 'append'
   end subroutine io_openfiles


!--------------------------------------------------------------------------
!  Close files written to during runtime
!--------------------------------------------------------------------------
   subroutine io_closefiles()
      if(mpi_rnk/=0) return
      if(io_dt/=0) close(io_dt)
      if(io_KE/=0) close(io_KE)
      if(io_ID/=0) close(io_ID)
      if(io_pt/=0) close(io_pt)
      if(io_fr/=0) close(io_fr)
   end subroutine io_closefiles


!--------------------------------------------------------------------------
!  Write to files
!--------------------------------------------------------------------------
   subroutine io_write2files()

      if(modulo(tim_step,i_save_rate1)==0) then
         call io_save_state_p()
         !call io_write_statistics()
         !call io_save_spectrum()
         !call io_save_meanprof()
         io_save1 = io_save1+1
      endif

      if(modulo(tim_step,i_save_rate2)==0) then
         if(io_KE/=0) call io_write_energy()
         if(io_ID/=0) call io_write_totEID()
         if(io_pt/=0) call io_write_pointvel()
         if(io_fr/=0) call io_write_friction()
         if(io_dt/=0 .and. d_timestep>0d0) call io_write_timestep()
         io_save2 = io_save2+1
      end if

      if(modulo(tim_step,i_save_rate3)==0) then        
         !call io_write_statistics()
         call io_write_velocity_statistics()
         call io_write_energy_statistics()
         io_save3 = io_save3+1
      endif

      if(io_dt/=0 .and. tim_new_dt) call io_write_timestep()

      if(modulo(tim_step,i_save_rate2*50)==0) then
         call io_closefiles()
         call io_openfiles()
      end if

   end subroutine io_write2files

!--------------------------------------------------------------------------
!  Load state - start from previous solution
!--------------------------------------------------------------------------
   subroutine io_load_state_p()
      integer :: e, f, i, rd
      integer :: n, N_,istp
      integer :: HIN(2, 0:i_pH1), PHM
      double precision :: d
      integer, parameter :: lda = 2*i_KL+1
      double precision :: A(lda,lda,i_N)
      double precision, allocatable :: r(:)
      !double precision     :: Dtmp(i_N, 0:i_pH1)
      double precision, allocatable :: Dtmp(:,:)
      logical :: interp
      integer :: K1,M1, nhi,nh,nh_,nh__,k,m
      integer :: K__, M__, Mp__  
      
      allocate( Dtmp(i_N, 0:i_pH1))
      e=nf90_open(io_statefile, IOR(nf90_nowrite,NF90_MPIIO), f,&
                comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)      
      if(e/=nf90_noerr) then
         if(mpi_rnk==0) open(99, file='PRECOMPUTING')
         if(mpi_rnk==0) close(99, status='delete')
         if(mpi_rnk==0) print*, 'state file not found!: '//io_statefile
#ifdef _MPI
         call mpi_barrier(mpi_comm_world, mpi_er)
         call mpi_finalize(mpi_er)
#endif
         stop 'io_load_state: file not found!'
      end if

      e=nf90_get_att(f,nf90_global,'t', d)
      if(d_time<0d0) tim_t = d
      if(mpi_rnk==0 .and. dabs(tim_t-d)>1d-8)  &
         print*,' t    :',d,' --> ', tim_t
      
      if(i_startstep<0) then 
        e=nf90_get_att(f,nf90_global,'istep', istp)
        tim_step=istp
        tim_steps=istp
        if(mpi_rnk==0) print*,' tim_step    :',istp,' --> ', tim_step
      end if

      e=nf90_get_att(f,nf90_global,'Re', d)
      if(mpi_rnk==0 .and. dabs(d_Re-d)>1d-8)  &
         print*,' Re   :',d,' --> ', d_Re
      e=nf90_get_att(f,nf90_global,'alpha', d)
      if(mpi_rnk==0 .and. dabs(d_alpha-d)>1d-8)  &
         print*,' alpha:',d,' --> ', d_alpha

      e=nf90_inq_dimid(f,'r', rd)
      e=nf90_inquire_dimension(f,rd, len=N_)

      allocate( r(1-i_KL:N_))
      e=nf90_inq_varid(f,'r', i)
      e=nf90_get_var(f,i, r(1:N_))
      r(1-i_KL:0) = -r(i_KL:1:-1)
      interp = (N_/=i_N)
      if(.not.interp) interp = (maxval(dabs(mes_D%r(:,1)-r(1:N_)))>1d-8)
      if(interp) then
         if(mpi_rnk==0) print*,' N    :', N_, ' --> ',i_N 
         stop 'interp not implmented!'
      end if

      e=nf90_inq_varid(f,'dt', i)
      e=nf90_get_var(f,i, d)
      if(d_timestep<0d0) tim_dt = d
      if(d_timestep>0d0) tim_dt = d_timestep
      e=nf90_inq_varid(f,'dtcor', i)
      if(e==nf90_noerr)  e=nf90_get_var(f,i, tim_corr_dt)
      e=nf90_inq_varid(f,'dtcfl', i)
      if(e==nf90_noerr)  e=nf90_get_var(f,i, tim_cfl_dt)

 
      HIN =0
      e=nf90_inq_varid(f,'Ur', i)
      if(e/=nf90_noerr)  print*, 'Field '//'Ur'//' not found!'
      if(e/=nf90_noerr)  stop 'io_load_coll'
      e=nf90_get_att(f,i, 'K',  K__)
      e=nf90_get_att(f,i, 'M',  M__)
      e=nf90_get_att(f,i,'Mp', Mp__)
      if(mpi_rnk==0) then
         if(K__ /=i_K)  print*, 'Ur', ' K :', K__, ' --> ',i_K 
         if(M__ /=i_M)  print*, 'Ur', ' M :', M__, ' --> ',i_M 
         if(Mp__/=i_Mp) print*, 'Ur', ' Mp:', Mp__,' --> ',i_Mp 
      end if
      if(Mp__/=i_Mp .and. M__>2) stop 'io_load_coll: Mp /= i_Mp'


      K1 = max(K__,i_K)-1
      M1 = max(M__,i_M)-1
      
      nhi= -1
      nh__ = -1
      nh = -1
      do m = 0, M1
         do k = -K1, K1
            if(k<0 .and. m==0) cycle
            if(abs(k)>=i_K .or. m>=i_M)  nh__ = nh__ + 1
            if(abs(k)>=i_K .or. m>=i_M)  cycle
            nh   = nh + 1
            if(abs(k)>=K__ .or. m>=M__)  cycle
            nh__ = nh__ + 1
            if(nh<var_H%pH0 .or. nh>var_H%pH0+var_H%pH1)  cycle
            nh_  = nh-var_H%pH0
            nhi   = nhi + 1
            HIN(1,nhi)=nh_
            HIN(2,nhi)=nh__
         end do
      end do
      do k = 0, _Np
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif
      if(mpi_rnk .eq. k) then
        !do m=0,i_pH1
        !print*, mpi_rnk,nhi,m,HIN(1,m),HIN(2,m)
        !end do
        print*, mpi_rnk,nhi, 0, maxval(HIN(1,:)),HIN(2,0),HIN(2,maxval(HIN(1,:)))
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif
      endif
      end do


      e=nf90_var_par_access(f, i, nf90_collective)
      !if (nhi.ge.0) then 
      e=nf90_get_var(f,i, Dtmp(1:N_,0:nhi),start=(/1,HIN(2,0)+1,1/))
      do m = 0, nhi
         vel_ur%Re(1:N_,HIN(1,m)) =Dtmp(1:N_,m)
      enddo 
      e=nf90_get_var(f,i, Dtmp(1:N_,0:nhi),start=(/1,HIN(2,0)+1,2/))
      do m = 0, nhi
         vel_ur%Im(1:N_,HIN(1,m)) =Dtmp(1:N_,m)
      enddo
      !endif 
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif

      e=nf90_inq_varid(f,'Ut', i)
      if(e/=nf90_noerr)  print*, 'Field '//'Ut'//' not found!'
      if(e/=nf90_noerr)  stop 'io_load_coll'
      e=nf90_var_par_access(f, i, nf90_collective)
      !if (nhi.ge.0) then 
      e=nf90_get_var(f,i, Dtmp(1:N_,0:nhi),start=(/1,HIN(2,0)+1,1/))
      do m = 0, nhi
         vel_ut%Re(1:N_,HIN(1,m)) =Dtmp(1:N_,m)
      enddo 
      e=nf90_get_var(f,i, Dtmp(1:N_,0:nhi),start=(/1,HIN(2,0)+1,2/))
      do m = 0, nhi
         vel_ut%Im(1:N_,HIN(1,m)) =Dtmp(1:N_,m)
      enddo
      !endif 
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif



      e=nf90_inq_varid(f,'Uz', i)
      if(e/=nf90_noerr)  print*, 'Field '//'Uz'//' not found!'
      if(e/=nf90_noerr)  stop 'io_load_coll'
      e=nf90_var_par_access(f, i, nf90_collective)
      !if (nhi.ge.0) then 
      e=nf90_get_var(f,i, Dtmp(1:N_,0:nhi),start=(/1,HIN(2,0)+1,1/))
      do m = 0, nhi
         vel_uz%Re(1:N_,HIN(1,m)) =Dtmp(1:N_,m)
      enddo 
      e=nf90_get_var(f,i, Dtmp(1:N_,0:nhi),start=(/1,HIN(2,0)+1,2/))
      do m = 0, nhi
         vel_uz%Im(1:N_,HIN(1,m)) =Dtmp(1:N_,m)
      enddo 
      !endif
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif
   
      deallocate(r)
      deallocate(Dtmp)


      e=nf90_close(f)
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif

   end subroutine io_load_state_p

!--------------------------------------------------------------------------
!  Load state - start from previous solution
!--------------------------------------------------------------------------
   subroutine io_load_state()
      integer :: e, f, i, rd
      integer :: n, N_,istp
      double precision :: d
      integer, parameter :: lda = 2*i_KL+1
      double precision :: A(lda,lda,i_N)
      double precision, allocatable :: r(:)
      logical :: interp

      e=nf90_open(io_statefile, nf90_nowrite, f)      
      if(e/=nf90_noerr) then
         if(mpi_rnk==0) open(99, file='PRECOMPUTING')
         if(mpi_rnk==0) close(99, status='delete')
         if(mpi_rnk==0) print*, 'state file not found!: '//io_statefile
#ifdef _MPI
         call mpi_barrier(mpi_comm_world, mpi_er)
         call mpi_finalize(mpi_er)
#endif
         stop 'io_load_state: file not found!'
      end if

      e=nf90_get_att(f,nf90_global,'t', d)
      if(d_time<0d0) tim_t = d
      if(mpi_rnk==0 .and. dabs(tim_t-d)>1d-8)  &
         print*,' t    :',d,' --> ', tim_t
      
      if(i_startstep<0) then 
	e=nf90_get_att(f,nf90_global,'istep', istp)
        tim_step=istp
        tim_steps=istp
	if(mpi_rnk==0) print*,' tim_step    :',istp,' --> ', tim_step
      end if

      e=nf90_get_att(f,nf90_global,'Re', d)
      if(mpi_rnk==0 .and. dabs(d_Re-d)>1d-8)  &
         print*,' Re   :',d,' --> ', d_Re
      e=nf90_get_att(f,nf90_global,'alpha', d)
      if(mpi_rnk==0 .and. dabs(d_alpha-d)>1d-8)  &
         print*,' alpha:',d,' --> ', d_alpha

      e=nf90_inq_dimid(f,'r', rd)
      e=nf90_inquire_dimension(f,rd, len=N_)

      allocate( r(1-i_KL:N_))
      e=nf90_inq_varid(f,'r', i)
      e=nf90_get_var(f,i, r(1:N_))
      r(1-i_KL:0) = -r(i_KL:1:-1)
      interp = (N_/=i_N)
      if(.not.interp) interp = (maxval(dabs(mes_D%r(:,1)-r(1:N_)))>1d-8)
      if(interp) then
         if(mpi_rnk==0) print*,' N    :', N_, ' --> ',i_N 
         call io_interp_wts(i_KL+N_,r,i_N,mes_D%r(1,1), A)
      end if

      e=nf90_inq_varid(f,'dt', i)
      e=nf90_get_var(f,i, d)
      if(d_timestep<0d0) tim_dt = d
      if(d_timestep>0d0) tim_dt = d_timestep
      e=nf90_inq_varid(f,'dtcor', i)
      if(e==nf90_noerr)  e=nf90_get_var(f,i, tim_corr_dt)
      e=nf90_inq_varid(f,'dtcfl', i)
      if(e==nf90_noerr)  e=nf90_get_var(f,i, tim_cfl_dt)

      call io_load_coll(f,'Ur',interp,N_,r,A,1, vel_ur)
      call io_load_coll(f,'Ut',interp,N_,r,A,1, vel_ut)
      call io_load_coll(f,'Uz',interp,N_,r,A,0, vel_uz)      

      deallocate(r)

      e=nf90_close(f)
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif

   end subroutine io_load_state


!--------------------------------------------------------------------------
!  Load coll variable
!--------------------------------------------------------------------------
   subroutine io_load_coll(f,nm,interp,N__,r,W,S, a)
      integer,          intent(in)  :: f
      character(*),     intent(in)  :: nm
      logical,          intent(in)  :: interp
      integer,          intent(in)  :: N__
      double precision, intent(in)  :: r(i_KL+N__)
      double precision, intent(in)  :: W(2*i_KL+1,2*i_KL+1,*)
      integer,          intent(in)  :: S
      type (coll),      intent(out) :: a
      double precision :: fn(1-i_KL:N__)
      integer :: K1,M1, nh,nh_,nh__,n,k,m
      integer :: K__, M__, Mp__      
      integer :: e,i
          
      e=nf90_inq_varid(f,nm, i)
      if(e/=nf90_noerr)  print*, 'Field '//nm//' not found!'
      if(e/=nf90_noerr)  stop 'io_load_coll'
      e=nf90_get_att(f,i, 'K',  K__)
      e=nf90_get_att(f,i, 'M',  M__)
      e=nf90_get_att(f,i,'Mp', Mp__)
      if(mpi_rnk==0) then
         if(K__ /=i_K)  print*, nm, ' K :', K__, ' --> ',i_K 
         if(M__ /=i_M)  print*, nm, ' M :', M__, ' --> ',i_M 
         if(Mp__/=i_Mp) print*, nm, ' Mp:', Mp__,' --> ',i_Mp 
      end if
      if(Mp__/=i_Mp .and. M__>2) stop 'io_load_coll: Mp /= i_Mp'

      a%Re = 0d0
      a%Im = 0d0

      K1 = max(K__,i_K)-1
      M1 = max(M__,i_M)-1

      nh__ = -1
      nh = -1
      do m = 0, M1
         do k = -K1, K1
            if(k<0 .and. m==0) cycle
            if(abs(k)>=i_K .or. m>=i_M)  nh__ = nh__ + 1
            if(abs(k)>=i_K .or. m>=i_M)  cycle
            nh   = nh + 1
            if(abs(k)>=K__ .or. m>=M__)  cycle
            nh__ = nh__ + 1
            if(nh<var_H%pH0 .or. nh>var_H%pH0+var_H%pH1)  cycle
            nh_  = nh-var_H%pH0
            if(interp) then
               e=nf90_get_var(f,i, fn(1:N__),start=(/1,nh__+1,1/))
               fn(:0) = fn(i_KL:1:-1)
               if(modulo(m*i_Mp+S,2)==1)  fn(:0) = -fn(:0)
               call io_interp(i_KL+N__,r,fn,W,i_N,mes_D%r(1,1), a%Re(1,nh_))
               e=nf90_get_var(f,i, fn(1:N__),start=(/1,nh__+1,2/))
               fn(:0) = fn(i_KL:1:-1)
               if(modulo(m*i_Mp+S,2)==1)  fn(:0) = -fn(:0)
               call io_interp(i_KL+N__,r,fn,W,i_N,mes_D%r(1,1), a%Im(1,nh_))
            else
               e=nf90_get_var(f,i, a%Re(1:N__,nh_),start=(/1,nh__+1,1/))
               e=nf90_get_var(f,i, a%Im(1:N__,nh_),start=(/1,nh__+1,2/))
            end if
         end do
      end do

   end subroutine io_load_coll


!--------------------------------------------------------------------------
!  interpolate; return weights
!--------------------------------------------------------------------------
   subroutine io_interp_wts(ni,xi,no,xo, A)
      integer, parameter :: lda = 2*i_KL+1
      integer,          intent(in)  :: ni, no
      double precision, intent(in)  :: xi(ni), xo(no)
      double precision, intent(out) :: A(lda,lda,no) 
      integer :: n,nn,i,j,l,r
      
      do n = 1, no
         j = 1
         do while(xi(j)<xo(n)-1d-8 .and. j<ni)
           j = j+1
         end do
         l = max(1,j-i_KL)
         r = min(j+i_KL,ni)
         nn = r-l+1
         do i = 1, nn
            A(i,1,n) = 1d0
         end do
         do j = 2, nn
            do i = 1, nn
               A(i,j,n) = A(i,j-1,n) * (xi(l+i-1)-xo(n)) / dble(j-1) 
            end do
         end do
         call mes_mat_invert(nn,A(1,1,n),lda)
      end do

   end subroutine io_interp_wts


!--------------------------------------------------------------------------
!  interpolate, given weights from io_interp_wts()
!--------------------------------------------------------------------------
   subroutine io_interp(ni,xi,fi,A,no,xo, fo)
      integer, parameter :: lda = 2*i_KL+1
      integer,          intent(in)  :: ni, no
      double precision, intent(in)  :: xi(ni), fi(ni), xo(no)
      double precision, intent(in)  :: A(lda,lda,no) 
      double precision, intent(out) :: fo(no)
      integer :: n,nn,i,j,l,r
      
      do n = 1, no
         j = 1
         do while(xi(j)<xo(n)-1d-8 .and. j<ni)
           j = j+1
         end do
         l = max(1,j-i_KL)
         r = min(j+i_KL,ni)
         nn = r-l+1
         fo(n) = dot_product(A(1,1:nn,n),fi(l:r))
      end do

   end subroutine io_interp


!--------------------------------------------------------------------------
!  Save state
!--------------------------------------------------------------------------
   subroutine io_save_state()
      character(4) :: cnum
      integer :: e, f
      integer :: rd, Hd, ReImd, dims(3)
      integer :: r,dt,dtcor,dtcfl, Ur,Ut,Uz

      write(cnum,'(I4.4)') io_save1

      if(mpi_rnk==0) then
         !print*, ' saving state'//cnum//'  t=', tim_t
         !e=nf90_create('state'//cnum//'.cdf.dat', nf90_clobber, f)
         e=nf90_create('state'//cnum//'.cdf.dat', NF90_64BIT_OFFSET, f)

         e=nf90_put_att(f, nf90_global, 't', tim_t)
         e=nf90_put_att(f, nf90_global, 'Re', d_Re)
         e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)

         e=nf90_def_dim(f, 'r', i_N, rd)
         e=nf90_def_dim(f, 'H', i_H1+1, Hd)
         e=nf90_def_dim(f, 'ReIm', 2, ReImd)

         e=nf90_def_var(f, 'r',     nf90_double, (/rd/), r)
         e=nf90_def_var(f, 'dt',    nf90_double, dt)
         e=nf90_def_var(f, 'dtcor', nf90_double, dtcor)
         e=nf90_def_var(f, 'dtcfl', nf90_double, dtcfl)

         dims = (/rd,Hd,ReImd/)
         call io_define_coll(f, 'Ur', dims, Ur)
         call io_define_coll(f, 'Ut', dims, Ut)         
         call io_define_coll(f, 'Uz', dims, Uz)         

         e=nf90_enddef(f)

         e=nf90_put_var(f, r, mes_D%r(1:i_N,1))
         e=nf90_put_var(f, dt, tim_dt)
         e=nf90_put_var(f, dtcor, tim_corr_dt)
         e=nf90_put_var(f, dtcfl, tim_cfl_dt)
      end if

      call io_save_coll(f,Ur, vel_ur)
      call io_save_coll(f,Ut, vel_ut)
      call io_save_coll(f,Uz, vel_uz)

      if(mpi_rnk==0)  &
         e=nf90_close(f)

   end subroutine io_save_state
 

!--------------------------------------------------------------------------
!  Save coll variable
!--------------------------------------------------------------------------
   subroutine io_define_coll(f,nm,dims, id)
      integer,      intent(in) :: f, dims(3)
      character(*), intent(in) :: nm
      integer, intent(out) :: id
      integer :: e
      e=nf90_def_var(f, nm, nf90_double, dims, id)
      e=nf90_put_att(f, id,  'K', i_K)      
      e=nf90_put_att(f, id,  'M', i_M)
      e=nf90_put_att(f, id, 'Mp', i_Mp)
   end subroutine io_define_coll

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine io_save_coll(f,id,a)
      integer,     intent(in) :: f, id
      type (coll), intent(in) :: a
      integer :: e
      
#ifndef _MPI
      e=nf90_put_var(f,id,a%Re(1:i_N,0:i_H1), start=(/1,1,1/))
      e=nf90_put_var(f,id,a%Im(1:i_N,0:i_H1), start=(/1,1,2/))

#else
      integer :: r, pH0,pH1

      if(mpi_rnk==0) then
         e=nf90_put_var(f,id,a%Re(1:i_N,0:var_H%pH1), start=(/1,1,1/))
         e=nf90_put_var(f,id,a%Im(1:i_N,0:var_H%pH1), start=(/1,1,2/))         
         do r = 1, mpi_sze-1
            pH0 = var_H%pH0_(r)
            pH1 = var_H%pH1_(r)
            mpi_tg = r
            call mpi_recv( c1%Re, i_N*(pH1+1), mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            call mpi_recv( c1%Im, i_N*(pH1+1), mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            e=nf90_put_var(f,id,c1%Re(1:i_N,0:pH1), start=(/1,pH0+1,1/))
            e=nf90_put_var(f,id,c1%Im(1:i_N,0:pH1), start=(/1,pH0+1,2/))
         end do
      else
         mpi_tg = mpi_rnk
         call mpi_send( a%Re, i_N*(var_H%pH1+1), mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
         call mpi_send( a%Im, i_N*(var_H%pH1+1), mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
      end if
#endif      
   end subroutine io_save_coll

!--------------------------------------------------------------------------
!  Save state parallel
!--------------------------------------------------------------------------
   subroutine io_save_state_p()
      character(8) :: cnum
      integer :: e, f
      integer :: rd, Hd, ReImd, dims(3),chunk_size(2)
      integer :: r,dt,dtcor,dtcfl, Ur,Ut,Uz
      integer :: pH0,pH1
      
      pH0 = var_H%pH0_(mpi_rnk)
      pH1 = var_H%pH1_(mpi_rnk)
      !print*, mpi_rnk,PH0,PH1
      !write(cnum,'(I4.4)') io_save1
      write(cnum,'(I8.8)') tim_step
 

      if(mpi_rnk==0) then
         print*, ' saving state'//cnum//'  t=', tim_t
         !print*, 'TOtal H',i_H1+1
      end if
         
         e=nf90_create('state'//cnum//'.cdf.dat', IOR(NF90_NETCDF4, NF90_MPIIO), f,&
                comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
         e=nf90_put_att(f, nf90_global, 'istep', tim_step)
         e=nf90_put_att(f, nf90_global, 't', tim_t)
         e=nf90_put_att(f, nf90_global, 'Re', d_Re)
         e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
         

         e=nf90_def_dim(f, 'r', i_N, rd)
         e=nf90_def_dim(f, 'H', i_H1+1, Hd)
         e=nf90_def_dim(f, 'ReIm', 2, ReImd)

         e=nf90_def_var(f, 'r',     nf90_double, (/rd/), r)
         e=nf90_def_var(f, 'dt',    nf90_double, dt)
         e=nf90_def_var(f, 'dtcor', nf90_double, dtcor)
         e=nf90_def_var(f, 'dtcfl', nf90_double, dtcfl)
         

         dims = (/rd,Hd,ReImd/)

         e=nf90_def_var(f, 'Ur', nf90_double, dims, Ur)
         e=nf90_put_att(f, Ur,  'K', i_K)      
         e=nf90_put_att(f, Ur,  'M', i_M)
         e=nf90_put_att(f, Ur, 'Mp', i_Mp)
         !if(mpi_rnk==0) print*, ' defining Ur '
         
         e=nf90_def_var(f, 'Ut', nf90_double, dims, Ut)
         e=nf90_put_att(f, Ut,  'K', i_K)      
         e=nf90_put_att(f, Ut,  'M', i_M)
         e=nf90_put_att(f, Ut, 'Mp', i_Mp)
         !if(mpi_rnk==0) print*, ' defining  Ut '
         
         e=nf90_def_var(f, 'Uz', nf90_double, dims, Uz)
         e=nf90_put_att(f, Uz,  'K', i_K)      
         e=nf90_put_att(f, Uz,  'M', i_M)
         e=nf90_put_att(f, Uz, 'Mp', i_Mp)
         !if(mpi_rnk==0) print*, ' defining  Uz '
        
         e=nf90_enddef(f)
         !if(mpi_rnk==0) print*, ' finish defining '
         
         e=nf90_put_var(f, r, mes_D%r(1:i_N,1))
         e=nf90_put_var(f, dt, tim_dt)
         e=nf90_put_var(f, dtcor, tim_corr_dt)
         e=nf90_put_var(f, dtcfl, tim_cfl_dt)
         !if(mpi_rnk==0) print*, ' sving   r '
         
         e=nf90_var_par_access(f, Ur, nf90_collective)
         e=nf90_put_var(f,Ur,vel_ur%Re,&
                start=(/1,pH0+1,1/),count=(/i_N,PH1+1/))
         e=nf90_put_var(f,Ur,vel_ur%Im,&
                 start=(/1,pH0+1,2/),count=(/i_N,PH1+1/))
         !if(mpi_rnk==0) print*, ' saving Ur '
   
         e=nf90_var_par_access(f, Ut, nf90_collective)
         e=nf90_put_var(f,Ut,vel_ut%Re,&
                start=(/1,pH0+1,1/),count=(/i_N,PH1+1/))
         e=nf90_put_var(f,Ut,vel_ut%Im,&
                 start=(/1,pH0+1,2/),count=(/i_N,PH1+1/))
         !if(mpi_rnk==0) print*, ' saving Ut '
         
         e=nf90_var_par_access(f, Uz, nf90_collective)
         e=nf90_put_var(f,Uz,vel_uz%Re,&
                start=(/1,pH0+1,1/),count=(/i_N,PH1+1/))
         e=nf90_put_var(f,Uz,vel_uz%Im,&
                 start=(/1,pH0+1,2/),count=(/i_N,PH1+1/))
         !if(mpi_rnk==0) print*, ' saving Uz '
  
      e=nf90_close(f)
   end subroutine io_save_state_p
!--------------------------------------------------------------------------
!  Save state parallel-final
!--------------------------------------------------------------------------
   subroutine io_save_state_pf()

      character(8) :: cnum
      integer :: e, f
      integer :: rd, Hd, ReImd, dims(3),chunk_size(2)
      integer :: r,dt,dtcor,dtcfl, Ur,Ut,Uz
      integer :: pH0,pH1

      pH0 = var_H%pH0_(mpi_rnk)
      pH1 = var_H%pH1_(mpi_rnk)
      !print*, mpi_rnk,PH0,PH1
      !write(cnum,'(I4.4)') io_save1
      write(cnum,'(I8.8)') tim_step

      if(mpi_rnk==0) then
         print*, ' saving state'//cnum//'  t=', tim_t
         !print*, 'TOtal H',i_H1+1
      end if

         e=nf90_create('state.cdf.in', IOR(NF90_NETCDF4, NF90_MPIIO), f,&
                comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
         e=nf90_put_att(f, nf90_global, 'istep', tim_step)
         e=nf90_put_att(f, nf90_global, 't', tim_t)
         e=nf90_put_att(f, nf90_global, 'Re', d_Re)
         e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)


         e=nf90_def_dim(f, 'r', i_N, rd)
         e=nf90_def_dim(f, 'H', i_H1+1, Hd)
         e=nf90_def_dim(f, 'ReIm', 2, ReImd)

         e=nf90_def_var(f, 'r',     nf90_double, (/rd/), r)
         e=nf90_def_var(f, 'dt',    nf90_double, dt)
         e=nf90_def_var(f, 'dtcor', nf90_double, dtcor)
         e=nf90_def_var(f, 'dtcfl', nf90_double, dtcfl)


                 dims = (/rd,Hd,ReImd/)

         e=nf90_def_var(f, 'Ur', nf90_double, dims, Ur)
         e=nf90_put_att(f, Ur,  'K', i_K)
         e=nf90_put_att(f, Ur,  'M', i_M)
         e=nf90_put_att(f, Ur, 'Mp', i_Mp)
         !if(mpi_rnk==0) print*, ' defining Ur '

         e=nf90_def_var(f, 'Ut', nf90_double, dims, Ut)
         e=nf90_put_att(f, Ut,  'K', i_K)
         e=nf90_put_att(f, Ut,  'M', i_M)
         e=nf90_put_att(f, Ut, 'Mp', i_Mp)
         !if(mpi_rnk==0) print*, ' defining  Ut '

         e=nf90_def_var(f, 'Uz', nf90_double, dims, Uz)
         e=nf90_put_att(f, Uz,  'K', i_K)
         e=nf90_put_att(f, Uz,  'M', i_M)
         e=nf90_put_att(f, Uz, 'Mp', i_Mp)
         !if(mpi_rnk==0) print*, ' defining  Uz '

         e=nf90_enddef(f)
         !if(mpi_rnk==0) print*, ' finish defining '

         e=nf90_put_var(f, r, mes_D%r(1:i_N,1))
         e=nf90_put_var(f, dt, tim_dt)
         e=nf90_put_var(f, dtcor, tim_corr_dt)
         e=nf90_put_var(f, dtcfl, tim_cfl_dt)
         !if(mpi_rnk==0) print*, ' sving   r '
         e=nf90_var_par_access(f, Ur, nf90_collective)
         e=nf90_put_var(f,Ur,vel_ur%Re,&
                start=(/1,pH0+1,1/),count=(/i_N,PH1+1/))
         e=nf90_put_var(f,Ur,vel_ur%Im,&
                 start=(/1,pH0+1,2/),count=(/i_N,PH1+1/))
         !if(mpi_rnk==0) print*, ' saving Ur '
         e=nf90_var_par_access(f, Ut, nf90_collective)
         e=nf90_put_var(f,Ut,vel_ut%Re,&
                start=(/1,pH0+1,1/),count=(/i_N,PH1+1/))
         e=nf90_put_var(f,Ut,vel_ut%Im,&
                 start=(/1,pH0+1,2/),count=(/i_N,PH1+1/))
         !if(mpi_rnk==0) print*, ' saving Ut '

         e=nf90_var_par_access(f, Uz, nf90_collective)
         e=nf90_put_var(f,Uz,vel_uz%Re,&
                start=(/1,pH0+1,1/),count=(/i_N,PH1+1/))
         e=nf90_put_var(f,Uz,vel_uz%Im,&
                 start=(/1,pH0+1,2/),count=(/i_N,PH1+1/))
         !if(mpi_rnk==0) print*, ' saving Uz '

      e=nf90_close(f)
   end subroutine io_save_state_pf

!--------------------------------------------------------------------------
!  save spectrum
!--------------------------------------------------------------------------
   subroutine io_save_spectrum()
      double precision :: n_(1:i_N), k_(0:i_K1), m_(0:i_M1)
      double precision :: n__(1:i_N), k__(0:i_K1), m__(0:i_M1)
      double precision,save :: TM(i_N,i_N), x(i_N)
      double precision :: d(i_N), dRe(i_N), dIm(i_N)
      logical, save :: set=.false.
      character(4) :: cnum
      integer :: i, n,kp
      _loop_km_vars
   10 format(i4,1e20.12)
      
      if(.not.set) then
         set =.true.
         do n = 0, i_N-1
            x(n+1) = 0.5d0 * ( 1d0 + dcos(d_PI*(i_N-n)/dble(i_N)) )
         end do
         do n = 1, i_N
            call cheby(n-1, 0, x, i_N, TM(1,n))
         end do
         call mes_mat_invert(i_N,TM,i_N)
         TM = transpose(TM)
      end if

      n_ = 0d0
      k_ = 0d0
      m_ = 0d0
      _loop_km_begin
         dRe = matmul(vel_ur%Re(:,nh), TM)
         dIm = matmul(vel_ur%Im(:,nh), TM)
         d = dRe*dRe+dIm*dIm
         dRe = matmul(vel_ut%Re(:,nh), TM)
         dIm = matmul(vel_ut%Im(:,nh), TM)
         d = max(d, dRe*dRe+dIm*dIm)
         dRe = matmul(vel_uz%Re(:,nh), TM)
         dIm = matmul(vel_uz%Im(:,nh), TM)
         d = max(d, dRe*dRe+dIm*dIm)
         d = dsqrt(d)
         kp = abs(k)
         do n = 1, i_N
             n_(n)  = max(d(n), n_(n))
             k_(kp) = max(d(n), k_(kp))
             m_(m)  = max(d(n), m_(m))            
         end do
      _loop_km_end
#ifdef _MPI
      call mpi_allreduce(n_, n__, i_N, mpi_double_precision,  &
         mpi_max, mpi_comm_world, mpi_er)
      n_ = n__
      call mpi_allreduce(k_, k__, i_K, mpi_double_precision,  &
         mpi_max, mpi_comm_world, mpi_er)
      k_ = k__
      call mpi_allreduce(m_, m__, i_M, mpi_double_precision,  &
         mpi_max, mpi_comm_world, mpi_er)
      m_ = m__
#endif

      if(mpi_rnk/=0) return
      write(cnum,'(I4.4)') io_save1
      open(11, status='unknown', file='vel_spec'//cnum//'.dat')
      write(11,*) '# t = ', tim_t
      write(11,*) '# n'
      do i = 1, i_N
         write(11,10) i, n_(i)      
      end do
      write(11,*)
      write(11,*) '# k'
      do i = 0, i_K1
         write(11,10) i, k_(i)      
      end do
      write(11,*)
      write(11,*) '# m'
      do i = 0, i_M1
         write(11,10) i*i_Mp, m_(i)      
      end do
      close(11)

   end subroutine io_save_spectrum


!--------------------------------------------------------------------------
!  save axial mean flow profile
!--------------------------------------------------------------------------
   subroutine io_save_meanprof()
      character(4) :: cnum
      integer :: n

      if(mpi_rnk/=0) return
      write(cnum,'(I4.4)') io_save1
      open(11, status='unknown', file='vel_prof'//cnum//'.dat')
      write(11,*) '# t = ', tim_t
      write(11,*) '# r  uz(r)'
      do n = 1, i_N
         write(11,'(2e20.12)')  mes_D%r(n,1),  &
            vel_uz%Re(n,0) + 1d0-mes_D%r(n,2)      
      end do
      close(11)

   end subroutine io_save_meanprof

 
!--------------------------------------------------------------------------
!  write to energy file
!--------------------------------------------------------------------------
   subroutine io_write_energy()
      double precision :: E,Ek0,Em0,E_, Ek(0:i_K1), Em(0:i_M1)

      call var_coll_norm(vel_ur, E,Ek,Em)
      Ek0 = Ek(0)
      Em0 = Em(0)
      call var_coll_norm(vel_ut, E_,Ek,Em)
      E   = E + E_
      Ek0 = Ek0 + Ek(0)
      Em0 = Em0 + Em(0)
      call var_coll_norm(vel_uz, E_,Ek,Em)
      E   = E + E_
      Ek0 = Ek0 + Ek(0)
      Em0 = Em0 + Em(0)
      
      if(mpi_rnk/=0) return
      write(io_KE,'(4e20.12)')  tim_t, E, Ek0, Em0
      
      if(E-Ek0>d_minE3d .or. tim_t<20d0) return
      print*, 'io_write_energy: Relaminarised!'
      open(99,file='RUNNING')
      close(99, status='delete')

   end subroutine io_write_energy


!--------------------------------------------------------------------------
!  (Total) Energy/E_lam, Input/I_lam, Dissip/D_lam
!--------------------------------------------------------------------------
   subroutine io_write_totEID()
      double precision :: E,E_,E__, Ek(0:i_K1), Em(0:i_M1)
      double precision :: Inp, Diss, Enrg, Lz

      call var_coll_copy(vel_uz, c3)
      if(mpi_rnk==0) c3%Re(:,0) = c3%Re(:,0) + vel_U
      call var_coll_norm(vel_ur, E,Ek,Em)
      call var_coll_norm(vel_ut, E_,Ek,Em)
      call var_coll_norm(c3,     E__,Ek,Em)
      Enrg = E + E_ + E__
      
      call var_coll_curl(vel_ur,vel_ut,c3, c1,c2,c3)
      call var_coll_norm(c1, E,Ek,Em)
      call var_coll_norm(c2, E_,Ek,Em)
      call var_coll_norm(c3, E__,Ek,Em)
      Diss = E + E_ + E__
      
      Lz = 2d0*d_PI/d_alpha
      Enrg = Enrg / (d_PI**2/(3d0*d_alpha))
      Diss = Diss * 2d0/d_Re
      Diss = Diss / abs(vel_Up(i_N)*d_PI*Lz/d_Re)
                                                   !   I/I_lam = 1+beta
      Inp = 1d0 + vel_Pr0
      
      if(mpi_rnk/=0) return
      write(io_ID,'(4e20.12)')  tim_t, Enrg, Inp, Diss

   end subroutine io_write_totEID


!--------------------------------------------------------------------------
!  write velocity at points.  r,t,z components at 3 points.
!--------------------------------------------------------------------------
   subroutine io_write_pointvel()
      integer :: n,rt(3),r(3)
      double precision :: d(9), d_
      if(_Ns/=1) stop 'io_write_pointvel: put _Ns=1'

      d = 1d0      
      do n = 1, i_N
         d_ = dabs(mes_D%r(n,1)-0.6d0) !0.6667d0)
         if(d_<d(1)) r(1) = n 
         if(d_<d(1)) d(1) = d_
         d_ = dabs(mes_D%r(n,1)-0.6d0) !0.7550d0)
         if(d_<d(2)) r(2) = n 
         if(d_<d(2)) d(2) = d_
         d_ = dabs(mes_D%r(n,1)-0.6d0) !0.9217d0)
         if(d_<d(3)) r(3) = n 
         if(d_<d(3)) d(3) = d_
      end do
      do n = 0, mpi_sze-1
         if(mes_D%pNi_(n)<=r(1)  &
            .and. mes_D%pNi_(n)+mes_D%pN_(n)-1>=r(1)) rt(1) = n
         if(mes_D%pNi_(n)<=r(2)  &
            .and. mes_D%pNi_(n)+mes_D%pN_(n)-1>=r(2)) rt(2) = n
         if(mes_D%pNi_(n)<=r(3)  &
            .and. mes_D%pNi_(n)+mes_D%pN_(n)-1>=r(3)) rt(3) = n
      end do

      if(rt(1)==mpi_rnk) then
         d(1) = vel_r%Re(0,0,r(1)-mes_D%pNi+1)
         d(2) = vel_t%Re(0,0,r(1)-mes_D%pNi+1)
         d(3) = vel_z%Re(0,0,r(1)-mes_D%pNi+1)
      end if      
      if(rt(2)==mpi_rnk) then
         d(4) = vel_r%Re(i_pZ/3,i_Th/3,r(2)-mes_D%pNi+1)
         d(5) = vel_t%Re(i_pZ/3,i_Th/3,r(2)-mes_D%pNi+1)
         d(6) = vel_z%Re(i_pZ/3,i_Th/3,r(2)-mes_D%pNi+1)
      end if
      if(rt(3)==mpi_rnk) then
         d(7) = vel_r%Re(2*i_pZ/3,2*i_Th/3,r(3)-mes_D%pNi+1)
         d(8) = vel_t%Re(2*i_pZ/3,2*i_Th/3,r(3)-mes_D%pNi+1)
         d(9) = vel_z%Re(2*i_pZ/3,2*i_Th/3,r(3)-mes_D%pNi+1)
      end if

#ifdef _MPI
      call mpi_bcast(d(1), 3, mpi_double_precision,  &
         rt(1), mpi_comm_world, mpi_er)
      call mpi_bcast(d(4), 3, mpi_double_precision,  &
         rt(2), mpi_comm_world, mpi_er)
      call mpi_bcast(d(7), 3, mpi_double_precision,  &
         rt(3), mpi_comm_world, mpi_er)
#endif
      if(mpi_rnk/=0) return
      if(tim_step==tim_steps)  write(io_pt,*) '# r=,', real(mes_D%r(r(1),1)),  &
         real(mes_D%r(r(2),1)), real(mes_D%r(r(3),1))
      write(io_pt,'(10e16.8)') tim_t, (d(n),n=1,9)

   end subroutine io_write_pointvel


!--------------------------------------------------------------------------
!  write:  1, time;  2,  bulk vel / excess pressure fraction if fixed flux;
!          3, mean uz at r=0;  4, friction velocity. 
!--------------------------------------------------------------------------
   subroutine io_write_friction()
      double precision :: Ub, Uc, Ufr

      Ub = 0.5d0 + 2d0*dot_product(vel_uz%Re(:,0),mes_D%intrdr)
      if(b_const_flux)  Ub = vel_Pr0
      Uc = 1d0 + dot_product(vel_uz%Re(1:1+i_KL,0),mes_D%dr0(:,0))
      Ufr = dot_product(vel_uz%Re(i_N-i_KL:i_N,0),mes_D%dr1(:,1))
      Ufr = dsqrt(dabs((Ufr-2d0)/d_Re))
      
      if(mpi_rnk/=0) return
      write(io_fr,'(4e20.12)') tim_t, Ub, Uc, Ufr
   
   end subroutine io_write_friction


!--------------------------------------------------------------------------
!  write to timestep file
!--------------------------------------------------------------------------
   subroutine io_write_timestep()
   10 format(3e18.10,I2)
   11 format(1e18.10,I8,3e17.9,I2)
      if(tim_new_dt) then
         if(mpi_rnk/=0) return
         write(io_dt,11) tim_t, tim_step,  &
            tim_dt, tim_corr_dt, tim_cfl_dt, tim_cfl_dir
      else
         call vel_maxtstep()
         if(mpi_rnk/=0) return
         write(io_dt,10) tim_t, tim_corr_dt, tim_cfl_dt, tim_cfl_dir
      end if
   end subroutine io_write_timestep


 
!--------------------------------------------------------------------------
!  subroutine to do statistics
!--------------------------------------------------------------------------
   subroutine io_write_velocity_statistics()
   use timestep
   implicit none
   double precision :: M_r(i_N), M_t(i_N),M_z(i_N),M_zt(i_N),M_or(i_N), M_ot(i_N),M_oz(i_N)
   double precision :: s_rr(i_N),s_rt(i_N),s_rz(i_N),s_tt(i_N),s_tz(i_N),s_zz(i_N)
   double precision :: s_rrr(i_N),s_ttt(i_N),s_zzz(i_N),s_r4(i_N),s_t4(i_N),s_z4(i_N)
   double precision :: s_r5(i_N),s_t5(i_N),s_z5(i_N),s_r6(i_N),s_t6(i_N),s_z6(i_N)
   double precision :: s_r7(i_N),s_t7(i_N),s_z7(i_N),s_r8(i_N),s_t8(i_N),s_z8(i_N)
   double precision :: s_r9(i_N),s_t9(i_N),s_z9(i_N),s_r10(i_N),s_t10(i_N),s_z10(i_N)
   double precision :: s_rtt(i_N),s_rzz(i_N),s_rrz(i_N),s_ttz(i_N)
   double precision :: s_orr(i_N),s_ort(i_N),s_orz(i_N),s_ott(i_N),s_otz(i_N),s_ozz(i_N)
   double precision :: s_rot(i_N),s_tor(i_N)
   

   double precision :: M_urr(i_N), M_utr(i_N),M_uzr(i_N),M_urt(i_N),M_utt(i_N),M_uzt(i_N),M_urz(i_N),M_utz(i_N),M_uzz(i_N)
   double precision :: S_urr(i_N), S_utr(i_N),S_uzr(i_N),S_urt(i_N),S_utt(i_N),S_uzt(i_N),S_urz(i_N),S_utz(i_N),S_uzz(i_N)
   double precision :: S_urruzr(i_N), S_urtuzt(i_N),S_urzuzz(i_N)
   double precision :: M_p(i_N), S_pp(i_N),S_pr(i_N),S_pt(i_N),S_pz(i_N),S_purr(i_N),S_putt(i_N),S_puzz(i_N),S_purzuzr(i_N)   


   double precision :: d1, d(i_N),dfact
   double precision :: Ub, Uc, Ufr

   double precision :: BCR(0:i_pH1), BCI(0:i_pH1)
   character(8) :: cnum
   integer :: e, f, kd,md,rd,rid,intrdrid,iKL1d,iKL2d,iKL3d,drid,dr0id,dr1id
   integer :: mrid,mtid,mzid,mztid,morid,motid,mozid,srrid,srtid,srzid,sttid,stzid,szzid
   integer :: srrrid,stttid,szzzid,srttid,srzzid,srrzid,sttzid,sr4id,st4id,sz4id
   integer :: sr5id,st5id,sz5id,sr6id,st6id,sz6id,sr7id,st7id,sz7id,sr8id,st8id,sz8id
   integer :: sr9id,st9id,sz9id,sr10id,st10id,sz10id
   integer :: sorrid,sortid,sorzid,sottid,sotzid,sozzid,srotid,storid
   integer :: murrid,mutrid,muzrid,murtid,muttid,muztid,murzid,mutzid,muzzid
   integer :: surrid,sutrid,suzrid,surtid,suttid,suztid,surzid,sutzid,suzzid, surruzrid, surtuztid, surzuzzid
   integer :: mpid,sppid,sprid,sptid,spzid,spurrid,sputtid,spuzzid,spurzuzrid

   integer :: n,n_, ff,fl,fi
   _loop_km_vars

   ! statistic : Bulk, center, and friction velocity
        
   Ub = 0.5d0 + 2d0*dot_product(vel_uz%Re(:,0),mes_D%intrdr)
   if(b_const_flux)  Ub = vel_Pr0
   Uc = 1d0 + dot_product(vel_uz%Re(1:1+i_KL,0),mes_D%dr0(:,0))
   Ufr = dot_product(vel_uz%Re(i_N-i_KL:i_N,0),mes_D%dr1(:,1))
   Ufr = dsqrt(dabs((Ufr-2d0)/d_Re))
      

      


   m_r = 0d0; m_t = 0d0; m_z = 0d0; m_or = 0d0;m_ot = 0d0; M_oz = 0d0
   s_rr= 0d0; s_rt= 0d0; s_rz= 0d0; s_tt= 0d0; s_tz= 0d0; s_zz= 0d0;
   s_rrr= 0d0; s_ttt= 0d0; s_zzz= 0d0; s_r4= 0d0; s_t4= 0d0; s_z4= 0d0;
   s_r5= 0d0; s_t5= 0d0; s_z5= 0d0;s_r6= 0d0; s_t6= 0d0; s_z6= 0d0;
   s_r7= 0d0; s_t7= 0d0; s_z7= 0d0;s_r8= 0d0; s_t8= 0d0; s_z8= 0d0;
   s_r9= 0d0; s_t9= 0d0; s_z9= 0d0;s_r10= 0d0; s_t10= 0d0; s_z10= 0d0;
   s_rtt= 0d0; s_rzz= 0d0; s_rrz= 0d0; s_ttz= 0d0;
   s_orr= 0d0; s_ort= 0d0; s_orz= 0d0; s_ott= 0d0; s_otz= 0d0; s_ozz= 0d0;
   s_rot= 0d0; s_tor= 0d0;

   M_urr=0d0;  M_utr=0d0;  M_uzr=0d0;  M_urt=0d0;  M_utt=0d0;  M_uzt=0d0;
   M_urz=0d0;  M_utz=0d0;  M_uzz=0d0;
   S_urr=0d0;  S_utr=0d0;  S_uzr=0d0;  S_urt=0d0;  S_utt=0d0;  S_uzt=0d0;
   S_urz=0d0;  S_utz=0d0;  S_uzz=0d0;  S_urruzr=0d0; S_urtuzt=0d0; S_urzuzz=0d0
   M_p=0; S_pp=0; S_pr=0; s_pt=0; s_pz=0; s_purr=0; s_putt=0; s_puzz=0; s_purzuzr=0

   
 
   write(cnum,'(I8.8)') tim_step
   if(mpi_rnk==0) then
         print*, ' saving statistics'//cnum//'  t=', tim_t
   end if

    ! statistic 1: Mean
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
         m_r(n_) = m_r(n_) + sum(vel_r%Re(:,:,n))
         m_t(n_) = m_t(n_) + sum(vel_t%Re(:,:,n))
         m_z(n_) = m_z(n_) + sum(vel_z%Re(:,:,n))
	 m_or(n_) = m_or(n_) + sum(vel_curlr%Re(:,:,n))
         m_ot(n_) = m_ot(n_) + sum(vel_curlt%Re(:,:,n))
	 M_oz(n_) = M_oz(n_) + sum(vel_curlz%Re(:,:,n))

         s_rr(n_) = s_rr(n_) + sum(vel_r%Re(:,:,n)**2)
         s_tt(n_) = s_tt(n_) + sum(vel_t%Re(:,:,n)**2)
	 s_zz(n_) = s_zz(n_) + sum(vel_z%Re(:,:,n)**2)
         s_rt(n_) = s_rt(n_) + sum(vel_r%Re(:,:,n)*vel_t%Re(:,:,n))
	 s_rz(n_) = s_rz(n_) + sum(vel_r%Re(:,:,n)*vel_z%Re(:,:,n))
         s_tz(n_) = s_tz(n_) + sum(vel_t%Re(:,:,n)*vel_z%Re(:,:,n))

	 s_rrr(n_) = s_rrr(n_) + sum(vel_r%Re(:,:,n)**3)
         s_ttt(n_) = s_ttt(n_) + sum(vel_t%Re(:,:,n)**3)
	 s_zzz(n_) = s_zzz(n_) + sum(vel_z%Re(:,:,n)**3)

	 s_rtt(n_) = s_rtt(n_) + sum(vel_r%Re(:,:,n)*vel_t%Re(:,:,n)**2)
         s_rzz(n_) = s_rzz(n_) + sum(vel_r%Re(:,:,n)*vel_z%Re(:,:,n)**2)
	 s_rrz(n_) = s_rrz(n_) + sum(vel_r%Re(:,:,n)**2*vel_z%Re(:,:,n))
	 s_ttz(n_) = s_ttz(n_) + sum(vel_t%Re(:,:,n)**2*vel_z%Re(:,:,n))

	 s_r4(n_) = s_r4(n_) + sum(vel_r%Re(:,:,n)**4)
         s_t4(n_) = s_t4(n_) + sum(vel_t%Re(:,:,n)**4)
	 s_z4(n_) = s_z4(n_) + sum(vel_z%Re(:,:,n)**4)
	 
	 s_r5(n_) = s_r5(n_) + sum(vel_r%Re(:,:,n)**5)
         s_t5(n_) = s_t5(n_) + sum(vel_t%Re(:,:,n)**5)
	 s_z5(n_) = s_z5(n_) + sum(vel_z%Re(:,:,n)**5)
	 
	 s_r6(n_) = s_r6(n_) + sum(vel_r%Re(:,:,n)**6)
         s_t6(n_) = s_t6(n_) + sum(vel_t%Re(:,:,n)**6)
	 s_z6(n_) = s_z6(n_) + sum(vel_z%Re(:,:,n)**6)
	 	 
	 s_r7(n_) = s_r7(n_) + sum(vel_r%Re(:,:,n)**7)
         s_t7(n_) = s_t7(n_) + sum(vel_t%Re(:,:,n)**7)
	 s_z7(n_) = s_z7(n_) + sum(vel_z%Re(:,:,n)**7)
	 	 
	 s_r8(n_) = s_r8(n_) + sum(vel_r%Re(:,:,n)**8)
         s_t8(n_) = s_t8(n_) + sum(vel_t%Re(:,:,n)**8)
	 s_z8(n_) = s_z8(n_) + sum(vel_z%Re(:,:,n)**8)
	 
	 s_r9(n_) = s_r9(n_) + sum(vel_r%Re(:,:,n)**9)
         s_t9(n_) = s_t9(n_) + sum(vel_t%Re(:,:,n)**9)
	 s_z9(n_) = s_z9(n_) + sum(vel_z%Re(:,:,n)**9)
	 
	 s_r10(n_) = s_r10(n_) + sum(vel_r%Re(:,:,n)**10)
         s_t10(n_) = s_t10(n_) + sum(vel_t%Re(:,:,n)**10)
	 s_z10(n_) = s_z10(n_) + sum(vel_z%Re(:,:,n)**10)
	 

	 s_orr(n_) = s_orr(n_) + sum(vel_curlr%Re(:,:,n)**2)
	 s_ott(n_) = s_ott(n_) + sum(vel_curlt%Re(:,:,n)**2)
	 s_ozz(n_) = s_ozz(n_) + sum(vel_curlz%Re(:,:,n)**2)
	 s_ort(n_) = s_ort(n_) + sum(vel_curlr%Re(:,:,n)*vel_curlt%Re(:,:,n))
	 s_orz(n_) = s_orz(n_) + sum(vel_curlr%Re(:,:,n)*vel_curlz%Re(:,:,n))
	 s_otz(n_) = s_otz(n_) + sum(vel_curlt%Re(:,:,n)*vel_curlz%Re(:,:,n))
 
         s_rot(n_) = s_rot(n_) + sum(vel_r%Re(:,:,n)*vel_curlt%Re(:,:,n))
	 s_tor(n_) = s_tor(n_) + sum(vel_t%Re(:,:,n)*vel_curlr%Re(:,:,n))
	
      end do

      
      !! Statistics 2: TKE
      !Laplacian of U
      call var_coll_curl(vel_ur,vel_ut,vel_uz, c1,c2,c3)
      call var_coll_curl(c1,c2,c3, c1,c2,c3)
      _loop_km_begin
        BCR(nh) = - c1%Re(i_N,nh)/d_Re
        BCI(nh) = - c1%Im(i_N,nh)/d_Re
      _loop_km_end
      ! r equation 
      !durdr
      call var_coll_meshmult(1,mes_D%dr(1),vel_ur,c1)
      call tra_coll2phys(c1,p1)
      pr%Re=vel_r%Re*p1%Re
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
	 m_urr(n_) = m_urr(n_) + sum(p1%Re(:,:,n))
	 s_urr(n_) = s_urr(n_) + sum(p1%Re(:,:,n)**2) 
      end do
      ! 1/r(durdt-ut)
      _loop_km_begin
            c1%Re(:,nh) = mes_D%r(:,-1)*(-vel_ur%Im(:,nh)*m*i_Mp-vel_ut%Re(:,nh))
            c1%Im(:,nh) = mes_D%r(:,-1)*( vel_ur%Re(:,nh)*m*i_Mp-vel_ut%Im(:,nh))   
      _loop_km_end
      call tra_coll2phys(c1,p1)
      pr%Re=pr%Re+vel_t%Re*p1%Re
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
	 m_urt(n_) = m_urt(n_) + sum(p1%Re(:,:,n))
	 s_urt(n_) = s_urt(n_) + sum(p1%Re(:,:,n)**2)
      end do
      ! durdz
      _loop_km_begin
            c1%Re(:,nh) = -vel_ur%Im(:,nh)*d_alpha*k
            c1%Im(:,nh) =  vel_ur%Re(:,nh)*d_alpha*k       
      _loop_km_end
      call tra_coll2phys(c1,p1)
      pr%Re=pr%Re+vel_z%Re*p1%Re
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
	 m_urz(n_) = m_urz(n_) + sum(p1%Re(:,:,n))
	 s_urz(n_) = s_urz(n_) + sum(p1%Re(:,:,n)**2)
	 pr%Re(:,:,n)=pr%Re(:,:,n)+p1%Re(:,:,n)*vel_U(n_)	
      end do
      call tra_phys2coll(pr,c1)

      ! theta equation
      ! dutdr
      call var_coll_meshmult(1,mes_D%dr(1),vel_ut,c2)
      call tra_coll2phys(c2,p1)
      pr%Re=vel_r%Re*p1%Re
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1	 
	 m_utr(n_) = m_utr(n_) + sum(p1%Re(:,:,n))	
	 s_utr(n_) = s_utr(n_) + sum(p1%Re(:,:,n)**2)	 
      end do
      ! 1/r(dutdt+ur)
      _loop_km_begin   
            c2%Re(:,nh) = mes_D%r(:,-1)*(-vel_ut%Im(:,nh)*m*i_Mp+vel_ur%Re(:,nh))
            c2%Im(:,nh) = mes_D%r(:,-1)*( vel_ut%Re(:,nh)*m*i_Mp+vel_ur%Im(:,nh))
      _loop_km_end
      call tra_coll2phys(c2,p1)
      pr%Re=pr%Re+vel_t%Re*p1%Re
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1	
	 m_utt(n_) = m_utt(n_) + sum(p1%Re(:,:,n))	 
	 s_utt(n_) = s_utt(n_) + sum(p1%Re(:,:,n)**2)
      end do
      ! dutdz
      _loop_km_begin          
            c2%Re(:,nh) = -vel_ut%Im(:,nh)*d_alpha*k
            c2%Im(:,nh) =  vel_ut%Re(:,nh)*d_alpha*k          
      _loop_km_end
      call tra_coll2phys(c2,p1)
      pr%Re=pr%Re+vel_z%Re*p1%Re
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1	
	 m_utz(n_) = m_utz(n_) + sum(p1%Re(:,:,n))	
	 s_utz(n_) = s_utz(n_) + sum(p1%Re(:,:,n)**2)
         pr%Re(:,:,n)=pr%Re(:,:,n)+p1%Re(:,:,n)*vel_U(n_)
      end do
      call tra_phys2coll(pr,c2)

      ! z equation
      ! duzdr 
      call var_coll_meshmult(0,mes_D%dr(1),vel_uz,c3)
      call tra_coll2phys(c3,p1)
      pr%Re=vel_r%Re*p1%Re
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1	
	 m_uzr(n_) = m_uzr(n_) + sum(p1%Re(:,:,n))	
	 s_uzr(n_) = s_uzr(n_) + sum(p1%Re(:,:,n)**2)	
      end do
      ! 1/r(duzdt)
      _loop_km_begin      
            c3%Re(:,nh) = mes_D%r(:,-1)*(-vel_uz%Im(:,nh)*m*i_Mp)
            c3%Im(:,nh) = mes_D%r(:,-1)*( vel_uz%Re(:,nh)*m*i_Mp)
      _loop_km_end
      call tra_coll2phys(c3,p1)
      pr%Re=pr%Re+vel_t%Re*p1%Re
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1	
	 m_uzt(n_) = m_uzt(n_) + sum(p1%Re(:,:,n))	
	 s_uzt(n_) = s_uzt(n_) + sum(p1%Re(:,:,n)**2)
      end do
      ! duzdz
      _loop_km_begin
            c3%Re(:,nh) = -vel_uz%Im(:,nh)*d_alpha*k
            c3%Im(:,nh) =  vel_uz%Re(:,nh)*d_alpha*k
      _loop_km_end
      call tra_coll2phys(c3,p1)
      pr%Re=pr%Re+vel_z%Re*p1%Re
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1	
	 m_uzz(n_) = m_uzz(n_) + sum(p1%Re(:,:,n))	
	 s_uzz(n_) = s_uzz(n_) + sum(p1%Re(:,:,n)**2)
         pr%Re(:,:,n)=pr%Re(:,:,n)+p1%Re(:,:,n)*vel_U(n_)+vel_r%Re(:,:,n)*vel_Up(n_)	
      end do
      call tra_phys2coll(pr,c3)

      _loop_km_begin
        BCR(nh) = BCR(nh) - c1%Re(i_N,nh)
        BCI(nh) = BCI(nh) - c1%Im(i_N,nh)
      _loop_km_end
      call var_coll_div(c1,c2,c3, c1)
      c1%Re=-c1%Re;
      c1%Im=-c1%Im;
      _loop_km_begin
        c1%Re(i_N,nh) = BCR(nh)
        c1%Im(i_N,nh) = BCI(nh)
      _loop_km_end

      !call tim_zerobc(c1)
      call tim_lumesh_inits( 0,1,0d0,1d0, LNps)
      call tim_lumesh_invert(0,LNps, c1)
      call tra_coll2phys(c1,pr)
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1      
	 M_p(n_)  = M_p(n_) + sum(pr%Re(:,:,n))
         s_pp(n_) = s_pp(n_) + sum(pr%Re(:,:,n)**2)
         s_pr(n_) = s_pr(n_) + sum(pr%Re(:,:,n)*vel_r%Re(:,:,n))
         s_pt(n_) = s_pt(n_) + sum(pr%Re(:,:,n)*vel_t%Re(:,:,n))
	 s_pz(n_) = s_pz(n_) + sum(pr%Re(:,:,n)*vel_z%Re(:,:,n))
	 
      end do

      !pdurdr
      call var_coll_meshmult(1,mes_D%dr(1),vel_ur,c1)
      call tra_coll2phys(c1,p1)  
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
	 s_purr(n_) = s_purr(n_) + sum(p1%Re(:,:,n)*pr%Re(:,:,n)) 
      end do
      ! p*1/r(dutdt+ur)
      _loop_km_begin   
            c2%Re(:,nh) = mes_D%r(:,-1)*(-vel_ut%Im(:,nh)*m*i_Mp+vel_ur%Re(:,nh))
            c2%Im(:,nh) = mes_D%r(:,-1)*( vel_ut%Re(:,nh)*m*i_Mp+vel_ur%Im(:,nh))
      _loop_km_end
      call tra_coll2phys(c2,p1)
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
	 s_putt(n_) = s_putt(n_) + sum(p1%Re(:,:,n)*pr%Re(:,:,n)) 
      end do
      ! pduzdz
      _loop_km_begin
            c3%Re(:,nh) = -vel_uz%Im(:,nh)*d_alpha*k
            c3%Im(:,nh) =  vel_uz%Re(:,nh)*d_alpha*k
      _loop_km_end
      call tra_coll2phys(c3,p1)
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
	 s_puzz(n_) = s_puzz(n_) + sum(p1%Re(:,:,n)*pr%Re(:,:,n)) 
      end do

     ! p(durdz+duzdr)
     ! durdz
      _loop_km_begin
            c1%Re(:,nh) = -vel_ur%Im(:,nh)*d_alpha*k
            c1%Im(:,nh) =  vel_ur%Re(:,nh)*d_alpha*k       
      _loop_km_end
      call tra_coll2phys(c1,p1)
      !duzdr
      call var_coll_meshmult(0,mes_D%dr(1),vel_uz,c3)
      call tra_coll2phys(c3,p2)

      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
	 s_purzuzr(n_) = s_purzuzr(n_) + &
		sum((p1%Re(:,:,n)+p2%Re(:,:,n))*pr%Re(:,:,n)) 
      end do

      ! viscous dissipation 
      !(durdr*duzdr)
      call var_coll_meshmult(1,mes_D%dr(1),vel_ur,c1)
      call tra_coll2phys(c1,p1)
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
	 S_urruzr(n_) = S_urruzr(n_) + sum(p1%Re(:,:,n)*p2%Re(:,:,n))
      end do

     !!(durdt*duzdt)
      _loop_km_begin
            c1%Re(:,nh) = mes_D%r(:,-1)*(-vel_ur%Im(:,nh)*m*i_Mp-vel_ut%Re(:,nh))
            c1%Im(:,nh) = mes_D%r(:,-1)*( vel_ur%Re(:,nh)*m*i_Mp-vel_ut%Im(:,nh))   
      _loop_km_end
      call tra_coll2phys(c1,p1)
      _loop_km_begin      
            c3%Re(:,nh) = mes_D%r(:,-1)*(-vel_uz%Im(:,nh)*m*i_Mp)
            c3%Im(:,nh) = mes_D%r(:,-1)*( vel_uz%Re(:,nh)*m*i_Mp)
      _loop_km_end
      call tra_coll2phys(c3,p2)
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
	 S_urtuzt(n_) = S_urtuzt(n_) + sum(p1%Re(:,:,n)*p2%Re(:,:,n))
      end do

     !(durdz*duzdz)
      _loop_km_begin
            c1%Re(:,nh) = -vel_ur%Im(:,nh)*d_alpha*k
            c1%Im(:,nh) =  vel_ur%Re(:,nh)*d_alpha*k       
      _loop_km_end
      call tra_coll2phys(c1,p1)
     ! duzdz
      _loop_km_begin
            c3%Re(:,nh) = -vel_uz%Im(:,nh)*d_alpha*k
            c3%Im(:,nh) =  vel_uz%Re(:,nh)*d_alpha*k
      _loop_km_end
      call tra_coll2phys(c3,p2)
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
	 S_urzuzz(n_) = S_urzuzz(n_) + sum(p1%Re(:,:,n)*p2%Re(:,:,n))
      end do

#ifdef _MPI
   call mpi_allreduce(m_r, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_r = d
   call mpi_allreduce(m_t, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_t = d
   call mpi_allreduce(m_z, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_z = d
   call mpi_allreduce(m_or, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_or = d
   call mpi_allreduce(m_ot, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_ot = d
   call mpi_allreduce(m_oz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_oz = d

   call mpi_allreduce(s_rr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_rr = d
   call mpi_allreduce(s_rt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_rt = d
   call mpi_allreduce(s_rz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_rz = d
   call mpi_allreduce(s_tt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_tt = d
   call mpi_allreduce(s_tz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_tz = d
   call mpi_allreduce(s_zz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_zz = d

   call mpi_allreduce(s_rrr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_rrr = d
   call mpi_allreduce(s_ttt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_ttt = d
   call mpi_allreduce(s_zzz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_zzz = d

   call mpi_allreduce(s_rtt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_rtt = d
   call mpi_allreduce(s_rzz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_rzz = d
   call mpi_allreduce(s_rrz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_rrz = d
   call mpi_allreduce(s_ttz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_ttz = d

   call mpi_allreduce(s_r4, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_r4 = d
   call mpi_allreduce(s_t4, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_t4 = d
   call mpi_allreduce(s_z4, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_z4 = d
   
   call mpi_allreduce(s_r5, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_r5 = d
   call mpi_allreduce(s_t5, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_t5 = d
   call mpi_allreduce(s_z5, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_z5 = d
   
   call mpi_allreduce(s_r6, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_r6 = d
   call mpi_allreduce(s_t6, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_t6 = d
   call mpi_allreduce(s_z6, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_z6 = d
   
   call mpi_allreduce(s_r7, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_r7 = d
   call mpi_allreduce(s_t7, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_t7 = d
   call mpi_allreduce(s_z7, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_z7 = d
   
   call mpi_allreduce(s_r8, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_r8 = d
   call mpi_allreduce(s_t8, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_t8 = d
   call mpi_allreduce(s_z8, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_z8 = d
   
   call mpi_allreduce(s_r9, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_r8 = d
   call mpi_allreduce(s_t9, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_t9 = d
   call mpi_allreduce(s_z9, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_z9 = d
   
   call mpi_allreduce(s_r10, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_r10 = d
   call mpi_allreduce(s_t10, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_t10 = d
   call mpi_allreduce(s_z10, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_z10 = d
   


   call mpi_allreduce(s_orr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_orr = d
   call mpi_allreduce(s_ort, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_ort = d
   call mpi_allreduce(s_orz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_orz = d
   call mpi_allreduce(s_ott, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_ott = d
   call mpi_allreduce(s_otz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_otz = d
   call mpi_allreduce(s_ozz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_ozz = d

   call mpi_allreduce(s_rot, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_rot = d
   call mpi_allreduce(s_tor, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_tor = d

   call mpi_allreduce(m_urr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_urr = d
   call mpi_allreduce(s_urr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_urr = d
   call mpi_allreduce(m_urt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_urt = d
   call mpi_allreduce(s_urt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_urt = d
   call mpi_allreduce(m_urz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_urz = d
   call mpi_allreduce(s_urz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_urz = d

   call mpi_allreduce(m_utr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_utr = d
   call mpi_allreduce(s_utr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_utr = d
   call mpi_allreduce(m_utt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_utt = d
   call mpi_allreduce(s_utt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_utt = d
   call mpi_allreduce(m_utz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_utz = d
   call mpi_allreduce(s_utz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_utz = d

   call mpi_allreduce(m_uzr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_uzr = d
   call mpi_allreduce(s_uzr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_uzr = d
   call mpi_allreduce(m_uzt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_uzt = d
   call mpi_allreduce(s_uzt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_uzt = d
   call mpi_allreduce(m_uzz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_uzz = d
   call mpi_allreduce(s_uzz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_uzz = d

 

   call mpi_allreduce(m_p, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_p = d
   call mpi_allreduce(s_pp, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_pp = d
   call mpi_allreduce(s_pr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_pr = d
   call mpi_allreduce(s_pt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_pt = d
   call mpi_allreduce(s_pz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_pz = d
   call mpi_allreduce(s_purr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_purr = d
   call mpi_allreduce(s_putt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_putt = d
   call mpi_allreduce(s_puzz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_puzz = d
   call mpi_allreduce(s_purzuzr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_purzuzr = d

   call mpi_allreduce(s_urruzr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_urruzr = d
   call mpi_allreduce(s_urtuzt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_urtuzt = d
   call mpi_allreduce(s_urzuzz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_urzuzz = d

   

#endif
    d1 = dble(i_Th*i_Z)
    m_r = m_r / d1
    m_t = m_t / d1
    m_z = m_z / d1
    m_or = m_or / d1
    m_ot = m_ot / d1
    m_oz = m_oz / d1
   
    s_rr = ( s_rr/d1 )
    s_rt = ( s_rt/d1 )
    s_rz = ( s_rz/d1 )
    s_tt = ( s_tt/d1 )
    s_tz = ( s_tz/d1 )
    s_zz = ( s_zz/d1 )

    s_rtt = ( s_rtt/d1 )
    s_rzz = ( s_rzz/d1 )
    s_rrz = ( s_rrz/d1 )
    s_ttz = ( s_ttz/d1 )

    s_rrr = ( s_rrr/d1 )
    s_ttt = ( s_ttt/d1 )
    s_zzz = ( s_zzz/d1 )

    s_r4 = ( s_r4/d1 )
    s_t4 = ( s_t4/d1 )
    s_z4 = ( s_z4/d1 )

    s_r5 = ( s_r5/d1 )
    s_t5 = ( s_t5/d1 )
    s_z5 = ( s_z5/d1 )
    
    s_r6 = ( s_r6/d1 )
    s_t6 = ( s_t6/d1 )
    s_z6 = ( s_z6/d1 )
    
    s_r7 = ( s_r7/d1 )
    s_t7 = ( s_t7/d1 )
    s_z7 = ( s_z7/d1 )
    
    s_r8 = ( s_r8/d1 )
    s_t8 = ( s_t8/d1 )
    s_z8 = ( s_z8/d1 )
    
    s_r9 = ( s_r9/d1 )
    s_t9 = ( s_t9/d1 )
    s_z9 = ( s_z9/d1 )
    
    s_r10 = ( s_r10/d1 )
    s_t10 = ( s_t10/d1 )
    s_z10 = ( s_z10/d1 )
    
    s_orr = ( s_orr/d1 )
    s_ort = ( s_ort/d1 )
    s_orz = ( s_orz/d1 )
    s_ott = ( s_ott/d1 )
    s_otz = ( s_otz/d1 )
    s_ozz = ( s_ozz/d1 )

    s_rot = ( s_rot/d1 )
    s_tor = ( s_tor/d1 )
    

    m_urr = m_urr / d1
    m_urt = m_urt / d1
    m_urz = m_urz / d1
    m_utr = m_utr / d1
    m_utt = m_utt / d1
    m_utz = m_utz / d1
    m_uzr = m_uzr / d1
    m_uzt = m_uzt / d1
    m_uzz = m_uzz / d1

    s_urr = (s_urr / d1  )
    s_urt = (s_urt / d1  )
    s_urz = (s_urz / d1  )
    s_utr = (s_utr / d1  )
    s_utt = (s_utt / d1  )
    s_utz = (s_utz / d1  )
    s_uzr = (s_uzr / d1  )
    s_uzt = (s_uzt / d1  )
    s_uzz = (s_uzz / d1  )
    
  
    m_p = m_p / d1
    s_pp = ( s_pp/d1 )
    s_pr = ( s_pr/d1 )
    s_pt = ( s_pt/d1 )
    s_pz = ( s_pz/d1 )
    s_purr = ( s_purr/d1 )
    s_putt = ( s_putt/d1 )
    s_puzz = ( s_puzz/d1 )
    s_purzuzr = ( s_purzuzr/d1 )

    s_urruzr = (s_urruzr / d1)
    s_urtuzt = (s_urtuzt / d1)
    s_urzuzz = (s_urzuzz / d1)
    
   
    m_zt = m_z + 1d0-mes_D%r(:,2)

   

     if(mpi_rnk==0) then
    	
        !print*, 'saving statistics'//cnum//'  t=', tim_t
    	e=nf90_create('Statist'//cnum//'.cdf.dat', NF90_64BIT_OFFSET, f)
    	e=nf90_put_att(f, nf90_global, 'Re', d_Re)
    	e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
    	e=nf90_put_att(f, nf90_global, 't', tim_t)
    	e=nf90_put_att(f, nf90_global, 'Ub', Ub)
    	e=nf90_put_att(f, nf90_global, 'Uc', Uc)
    	e=nf90_put_att(f, nf90_global, 'Ut', Ufr)    	
    	    	
    	e=nf90_def_dim(f, 'r', i_N, rd)
    	e=nf90_def_dim(f, 'K', 2*i_K-1, Kd)
    	e=nf90_def_dim(f, 'M', i_M, Md)
        e=nf90_def_dim(f, 'iKL1', 2*i_KL+1, iKL1d)
        e=nf90_def_dim(f, 'iKL2', i_N+i_KL, iKL2d)
        e=nf90_def_dim(f, 'iKL3', i_KL+1, iKL3d)
    
    	e=nf90_def_var(f, 'r', nf90_double, (/rd/), rid)
        e=nf90_def_var(f, 'intrdr', nf90_double, (/rd/), intrdrid)
        e=nf90_def_var(f, 'dr0', nf90_double, (/iKL3d,iKL3d/), dr0id)
        e=nf90_def_var(f, 'dr1', nf90_double, (/iKL3d,iKL3d/), dr1id)

    	e=nf90_def_var(f, 'Ur', nf90_double, (/rd/), mrid)
    	e=nf90_def_var(f, 'Ut', nf90_double, (/rd/), mtid)
        e=nf90_def_var(f, 'Uz', nf90_double, (/rd/), mztid)
        e=nf90_def_var(f, 'Uzp', nf90_double, (/rd/), mzid)
        e=nf90_def_var(f, 'Omegar', nf90_double, (/rd/), morid)
    	e=nf90_def_var(f, 'Omegat', nf90_double, (/rd/), motid)
        e=nf90_def_var(f, 'Omegaz', nf90_double, (/rd/), mozid)

        e=nf90_def_var(f, 'urr', nf90_double, (/rd/), srrid)
	e=nf90_def_var(f, 'urt', nf90_double, (/rd/), srtid)
	e=nf90_def_var(f, 'urz', nf90_double, (/rd/), srzid)
	e=nf90_def_var(f, 'utt', nf90_double, (/rd/), sttid)
	e=nf90_def_var(f, 'utz', nf90_double, (/rd/), stzid)
	e=nf90_def_var(f, 'uzz', nf90_double, (/rd/), szzid)

	e=nf90_def_var(f, 'urrr', nf90_double, (/rd/), srrrid)
	e=nf90_def_var(f, 'uttt', nf90_double, (/rd/), stttid)
	e=nf90_def_var(f, 'uzzz', nf90_double, (/rd/), szzzid)
	e=nf90_def_var(f, 'urtt', nf90_double, (/rd/), srttid)
	e=nf90_def_var(f, 'urzz', nf90_double, (/rd/), srzzid)
	e=nf90_def_var(f, 'urrz', nf90_double, (/rd/), srrzid)
	e=nf90_def_var(f, 'uttz', nf90_double, (/rd/), sttzid)
	e=nf90_def_var(f, 'ur4', nf90_double, (/rd/), sr4id)
	e=nf90_def_var(f, 'ut4', nf90_double, (/rd/), st4id)
	e=nf90_def_var(f, 'uz4', nf90_double, (/rd/), sz4id)
	e=nf90_def_var(f, 'ur5', nf90_double, (/rd/), sr5id)
	e=nf90_def_var(f, 'ut5', nf90_double, (/rd/), st5id)
	e=nf90_def_var(f, 'uz5', nf90_double, (/rd/), sz5id)
	e=nf90_def_var(f, 'ur6', nf90_double, (/rd/), sr6id)
	e=nf90_def_var(f, 'ut6', nf90_double, (/rd/), st6id)
	e=nf90_def_var(f, 'uz6', nf90_double, (/rd/), sz6id)
	e=nf90_def_var(f, 'ur7', nf90_double, (/rd/), sr7id)
	e=nf90_def_var(f, 'ut7', nf90_double, (/rd/), st7id)
	e=nf90_def_var(f, 'uz7', nf90_double, (/rd/), sz7id)
	e=nf90_def_var(f, 'ur8', nf90_double, (/rd/), sr8id)
	e=nf90_def_var(f, 'ut8', nf90_double, (/rd/), st8id)
	e=nf90_def_var(f, 'uz8', nf90_double, (/rd/), sz8id)
	e=nf90_def_var(f, 'ur9', nf90_double, (/rd/), sr9id)
	e=nf90_def_var(f, 'ut9', nf90_double, (/rd/), st9id)
	e=nf90_def_var(f, 'uz9', nf90_double, (/rd/), sz9id)
	e=nf90_def_var(f, 'ur10', nf90_double, (/rd/), sr10id)
	e=nf90_def_var(f, 'ut10', nf90_double, (/rd/), st10id)
	e=nf90_def_var(f, 'uz10', nf90_double, (/rd/), sz10id)

	e=nf90_def_var(f, 'orr', nf90_double, (/rd/), sorrid)
	e=nf90_def_var(f, 'ort', nf90_double, (/rd/), sortid)
	e=nf90_def_var(f, 'orz', nf90_double, (/rd/), sorzid)
	e=nf90_def_var(f, 'ott', nf90_double, (/rd/), sottid)
	e=nf90_def_var(f, 'otz', nf90_double, (/rd/), sotzid)
	e=nf90_def_var(f, 'ozz', nf90_double, (/rd/), sozzid)

	e=nf90_def_var(f, 'rot', nf90_double, (/rd/), srotid)
	e=nf90_def_var(f, 'tor', nf90_double, (/rd/), storid)

	e=nf90_def_var(f, 'dmurr', nf90_double, (/rd/), murrid)   
        e=nf90_def_var(f, 'dmurt', nf90_double, (/rd/), murtid)  
	e=nf90_def_var(f, 'dmurz', nf90_double, (/rd/), murzid)  
	e=nf90_def_var(f, 'dmutr', nf90_double, (/rd/), mutrid)   
        e=nf90_def_var(f, 'dmutt', nf90_double, (/rd/), muttid)  
	e=nf90_def_var(f, 'dmutz', nf90_double, (/rd/), mutzid)  
	e=nf90_def_var(f, 'dmuzr', nf90_double, (/rd/), muzrid)   
        e=nf90_def_var(f, 'dmuzt', nf90_double, (/rd/), muztid)  
	e=nf90_def_var(f, 'dmuzz', nf90_double, (/rd/), muzzid)  

	e=nf90_def_var(f, 'dsurr', nf90_double, (/rd/), surrid)   
        e=nf90_def_var(f, 'dsurt', nf90_double, (/rd/), surtid)  
	e=nf90_def_var(f, 'dsurz', nf90_double, (/rd/), surzid)  
	e=nf90_def_var(f, 'dsutr', nf90_double, (/rd/), sutrid)   
        e=nf90_def_var(f, 'dsutt', nf90_double, (/rd/), suttid)  
	e=nf90_def_var(f, 'dsutz', nf90_double, (/rd/), sutzid)  
	e=nf90_def_var(f, 'dsuzr', nf90_double, (/rd/), suzrid)   
        e=nf90_def_var(f, 'dsuzt', nf90_double, (/rd/), suztid)  
	e=nf90_def_var(f, 'dsuzz', nf90_double, (/rd/), suzzid)  


    	e=nf90_def_var(f, 'P', nf90_double, (/rd/), mpid)   
        e=nf90_def_var(f, 'pp', nf90_double, (/rd/), sppid)
	e=nf90_def_var(f, 'pr', nf90_double, (/rd/), sprid)
	e=nf90_def_var(f, 'pt', nf90_double, (/rd/), sptid)
	e=nf90_def_var(f, 'pz', nf90_double, (/rd/), spzid)
        e=nf90_def_var(f, 'purr', nf90_double, (/rd/), spurrid)
        e=nf90_def_var(f, 'putt', nf90_double, (/rd/), sputtid)
        e=nf90_def_var(f, 'puzz', nf90_double, (/rd/), spuzzid)
	e=nf90_def_var(f, 'purzuzr', nf90_double, (/rd/), spurzuzrid)

	e=nf90_def_var(f, 'dsurruzr', nf90_double, (/rd/), surruzrid)   
        e=nf90_def_var(f, 'dsurtuzt', nf90_double, (/rd/), surtuztid)  
	e=nf90_def_var(f, 'dsurzuzz', nf90_double, (/rd/), surzuzzid)  

    	e=nf90_enddef(f)


    	e=nf90_put_var(f, rid, mes_D%r(1:i_N,1))
        e=nf90_put_var(f, intrdrid, mes_D%intrdr)
        e=nf90_put_var(f, dr0id, mes_D%dr0, start=(/1,1/))
        e=nf90_put_var(f, dr1id, mes_D%dr1, start=(/1,1/))

        e=nf90_put_var(f, mrid, m_r)
	e=nf90_put_var(f, mtid, m_t)
	e=nf90_put_var(f, mzid, m_z)
	e=nf90_put_var(f, mztid, m_zt)
	e=nf90_put_var(f, morid, m_or)
	e=nf90_put_var(f, motid, m_ot)
	e=nf90_put_var(f, mozid, m_oz)
       
	e=nf90_put_var(f, srrid, s_rr)
	e=nf90_put_var(f, srtid, s_rt)
	e=nf90_put_var(f, srzid, s_rz)
	e=nf90_put_var(f, sttid, s_tt)
	e=nf90_put_var(f, stzid, s_tz)
	e=nf90_put_var(f, szzid, s_zz)

	e=nf90_put_var(f, srrrid, s_rrr)
	e=nf90_put_var(f, stttid, s_ttt)
	e=nf90_put_var(f, szzzid, s_zzz)
	e=nf90_put_var(f, srttid, s_rtt)
	e=nf90_put_var(f, srzzid, s_rzz)
	e=nf90_put_var(f, srrzid, s_rrz)
	e=nf90_put_var(f, sttzid, s_ttz)
	e=nf90_put_var(f, sr4id, s_r4)
	e=nf90_put_var(f, st4id, s_t4)
	e=nf90_put_var(f, sz4id, s_z4)
	e=nf90_put_var(f, sr5id, s_r5)
	e=nf90_put_var(f, st5id, s_t5)
	e=nf90_put_var(f, sz5id, s_z5)
	e=nf90_put_var(f, sr6id, s_r6)
	e=nf90_put_var(f, st6id, s_t6)
	e=nf90_put_var(f, sz6id, s_z6)
	e=nf90_put_var(f, sr7id, s_r7)
	e=nf90_put_var(f, st7id, s_t7)
	e=nf90_put_var(f, sz7id, s_z7)
	e=nf90_put_var(f, sr8id, s_r8)
	e=nf90_put_var(f, st8id, s_t8)
	e=nf90_put_var(f, sz8id, s_z8)
	e=nf90_put_var(f, sr9id, s_r9)
	e=nf90_put_var(f, st9id, s_t9)
	e=nf90_put_var(f, sz9id, s_z9)
	e=nf90_put_var(f, sr10id, s_r10)
	e=nf90_put_var(f, st10id, s_t10)
	e=nf90_put_var(f, sz10id, s_z10)




	e=nf90_put_var(f, sorrid, s_orr)
	e=nf90_put_var(f, sortid, s_ort)
	e=nf90_put_var(f, sorzid, s_orz)
	e=nf90_put_var(f, sottid, s_ott)
	e=nf90_put_var(f, sotzid, s_otz)
	e=nf90_put_var(f, sozzid, s_ozz)
	e=nf90_put_var(f, srotid, s_rot)
	e=nf90_put_var(f, storid, s_tor)

	e=nf90_put_var(f, murrid, m_urr) 
	e=nf90_put_var(f, murtid, m_urt) 
	e=nf90_put_var(f, murzid, m_urz) 
	e=nf90_put_var(f, mutrid, m_utr) 
	e=nf90_put_var(f, muttid, m_utt) 
	e=nf90_put_var(f, mutzid, m_utz) 
	e=nf90_put_var(f, muzrid, m_uzr) 
	e=nf90_put_var(f, muztid, m_uzt) 
	e=nf90_put_var(f, muzzid, m_uzz) 

	e=nf90_put_var(f, surrid, s_urr) 
	e=nf90_put_var(f, surtid, s_urt) 
	e=nf90_put_var(f, surzid, s_urz) 
	e=nf90_put_var(f, sutrid, s_utr) 
	e=nf90_put_var(f, suttid, s_utt) 
	e=nf90_put_var(f, sutzid, s_utz) 
	e=nf90_put_var(f, suzrid, s_uzr) 
	e=nf90_put_var(f, suztid, s_uzt) 
	e=nf90_put_var(f, suzzid, s_uzz) 

	e=nf90_put_var(f, mpid, m_p)    
	e=nf90_put_var(f, sppid, s_pp)
	e=nf90_put_var(f, sprid, s_pr)
	e=nf90_put_var(f, sptid, s_pt)
	e=nf90_put_var(f, spzid, s_pz)
        e=nf90_put_var(f, spurrid, s_purr)
	e=nf90_put_var(f, sputtid, s_putt)
	e=nf90_put_var(f, spuzzid, s_puzz)
	e=nf90_put_var(f, spurzuzrid, s_purzuzr)

	e=nf90_put_var(f, surruzrid, s_urruzr) 
	e=nf90_put_var(f, surtuztid, s_urtuzt) 
	e=nf90_put_var(f, surzuzzid, s_urzuzz) 

     

        e=nf90_close(f)
    end if
    
   
    s_rr = ( s_rr - m_r**2 )
    s_rt = ( s_rt - m_r*m_t)
    s_rz = ( s_rz - m_r*m_z)
    s_tt = ( s_tt - m_t**2 )
    s_tz = ( s_tz - m_t*m_z)
    s_zz = ( s_zz - m_z**2 )

    s_rtt = ( s_rtt- m_r*m_t**2)
    s_rzz = ( s_rzz- m_r*m_z**2-m_r*s_zz-2*m_z*s_rz)
    s_rrz = ( s_rrz - m_r**2*m_z-s_rr*m_z)
    s_ttz = ( s_ttz - m_t**2*m_z-s_tt*m_z)

    s_rrr = ( s_rrr - m_r**3)
    s_ttt = ( s_ttt - m_t**3)
    s_zzz = ( s_zzz -3*s_zz*m_z+2*m_z**3)

    s_r4 = ( s_r4 - m_r**4)
    s_t4 = ( s_t4- m_t**4)
    s_z4 = ( s_z4 - m_z**4-6*m_z**2*s_zz-4*m_z*s_zzz)
    
    s_r5 = ( s_r5 - m_r**5)
    s_t5 = ( s_t5- m_t**5)
    s_z5 = ( s_z5 -5*s_z4*m_z+10*s_zzz*m_z**2-10*s_zz*m_z**3+4*m_z**5)

    s_r6 = ( s_r6 - m_r**6)
    s_t6 = ( s_t6- m_t**6)
    s_z6 = ( s_z6 -6*s_z5*m_z+15*s_z4*m_z**2-20*s_zzz*m_z**3+15*s_zz*m_z**4-5*m_z**6)

    
  
    s_orr = ( s_orr - m_or**2 )
    s_ort = ( s_ort - m_or*m_ot)
    s_orz = ( s_orz- m_or*m_oz)
    s_ott = ( s_ott - m_ot**2 )
    s_otz = ( s_otz - m_ot*m_oz)
    s_ozz = ( s_ozz - m_oz**2 )

    s_rot = ( s_rot - m_ot*m_r)
    s_tor = ( s_tor - m_or*m_t)

    
    s_urr = (s_urr - m_urr**2 )
    s_urt = (s_urt - m_urt**2 )
    s_urz = (s_urz - m_urz**2 )
    s_utr = (s_utr - m_utr**2 )
    s_utt = (s_utt - m_utt**2 )
    s_utz = (s_utz - m_utz**2 )
    s_uzr = (s_uzr - m_uzr**2 )
    s_uzt = (s_uzt - m_uzt**2 )
    s_uzz = (s_uzz - m_uzz**2 )

   
    
    
    s_pp = ( s_pp - m_p**2 )
    s_pr = ( s_pr - m_p*m_r )
    s_pt = ( s_pt - m_p*m_t )
    s_pz = ( s_pz - m_p*m_z )
    s_purr = ( s_purr- m_p*m_urr )
    s_putt = ( s_putt - m_p*m_utt )
    s_puzz = ( s_puzz - m_p*m_uzz )
    s_purzuzr = ( s_purzuzr - m_p*(m_urz+m_uzr))

    s_urruzr = (s_urruzr- m_urr* m_uzr)
    s_urtuzt = (s_urtuzt- m_urt* m_uzt)
    s_urzuzz = (s_urzuzz- m_urz* m_uzz)
    
     if(mpi_rnk==0) then
    
    	e=nf90_create('Statistp'//cnum//'.cdf.dat', NF90_64BIT_OFFSET, f)
    	e=nf90_put_att(f, nf90_global, 'Re', d_Re)
    	e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
    	e=nf90_put_att(f, nf90_global, 't', tim_t)    	
    	e=nf90_def_dim(f, 'r', i_N, rd)
    	e=nf90_def_dim(f, 'K', 2*i_K-1, Kd)
    	e=nf90_def_dim(f, 'M', i_M, Md)
        e=nf90_def_dim(f, 'iKL1', 2*i_KL+1, iKL1d)
        e=nf90_def_dim(f, 'iKL2', i_N+i_KL, iKL2d)
        e=nf90_def_dim(f, 'iKL3', i_KL+1, iKL3d)
    
    	e=nf90_def_var(f, 'r', nf90_double, (/rd/), rid)
        e=nf90_def_var(f, 'intrdr', nf90_double, (/rd/), intrdrid)
        e=nf90_def_var(f, 'dr0', nf90_double, (/iKL3d,iKL3d/), dr0id)
        e=nf90_def_var(f, 'dr1', nf90_double, (/iKL3d,iKL3d/), dr1id)

    	e=nf90_def_var(f, 'Ur', nf90_double, (/rd/), mrid)
    	e=nf90_def_var(f, 'Ut', nf90_double, (/rd/), mtid)
        e=nf90_def_var(f, 'Uz', nf90_double, (/rd/), mztid)
        e=nf90_def_var(f, 'Uzp', nf90_double, (/rd/), mzid)
        e=nf90_def_var(f, 'Omegar', nf90_double, (/rd/), morid)
    	e=nf90_def_var(f, 'Omegat', nf90_double, (/rd/), motid)
        e=nf90_def_var(f, 'Omegaz', nf90_double, (/rd/), mozid)

        e=nf90_def_var(f, 'urr', nf90_double, (/rd/), srrid)
	e=nf90_def_var(f, 'urt', nf90_double, (/rd/), srtid)
	e=nf90_def_var(f, 'urz', nf90_double, (/rd/), srzid)
	e=nf90_def_var(f, 'utt', nf90_double, (/rd/), sttid)
	e=nf90_def_var(f, 'utz', nf90_double, (/rd/), stzid)
	e=nf90_def_var(f, 'uzz', nf90_double, (/rd/), szzid)

	e=nf90_def_var(f, 'urrr', nf90_double, (/rd/), srrrid)
	e=nf90_def_var(f, 'uttt', nf90_double, (/rd/), stttid)
	e=nf90_def_var(f, 'uzzz', nf90_double, (/rd/), szzzid)
	e=nf90_def_var(f, 'urtt', nf90_double, (/rd/), srttid)
	e=nf90_def_var(f, 'urzz', nf90_double, (/rd/), srzzid)
	e=nf90_def_var(f, 'urrz', nf90_double, (/rd/), srrzid)
	e=nf90_def_var(f, 'uttz', nf90_double, (/rd/), sttzid)
	e=nf90_def_var(f, 'ur4', nf90_double, (/rd/), sr4id)
	e=nf90_def_var(f, 'ut4', nf90_double, (/rd/), st4id)
	e=nf90_def_var(f, 'uz4', nf90_double, (/rd/), sz4id)
	e=nf90_def_var(f, 'ur5', nf90_double, (/rd/), sr5id)
	e=nf90_def_var(f, 'ut5', nf90_double, (/rd/), st5id)
	e=nf90_def_var(f, 'uz5', nf90_double, (/rd/), sz5id)
	e=nf90_def_var(f, 'ur6', nf90_double, (/rd/), sr6id)
	e=nf90_def_var(f, 'ut6', nf90_double, (/rd/), st6id)
	e=nf90_def_var(f, 'uz6', nf90_double, (/rd/), sz6id)
	e=nf90_def_var(f, 'ur7', nf90_double, (/rd/), sr7id)
	e=nf90_def_var(f, 'ut7', nf90_double, (/rd/), st7id)
	e=nf90_def_var(f, 'uz7', nf90_double, (/rd/), sz7id)
	e=nf90_def_var(f, 'ur8', nf90_double, (/rd/), sr8id)
	e=nf90_def_var(f, 'ut8', nf90_double, (/rd/), st8id)
	e=nf90_def_var(f, 'uz8', nf90_double, (/rd/), sz8id)
	e=nf90_def_var(f, 'ur9', nf90_double, (/rd/), sr9id)
	e=nf90_def_var(f, 'ut9', nf90_double, (/rd/), st9id)
	e=nf90_def_var(f, 'uz9', nf90_double, (/rd/), sz9id)
	e=nf90_def_var(f, 'ur10', nf90_double, (/rd/), sr10id)
	e=nf90_def_var(f, 'ut10', nf90_double, (/rd/), st10id)
	e=nf90_def_var(f, 'uz10', nf90_double, (/rd/), sz10id)

	e=nf90_def_var(f, 'orr', nf90_double, (/rd/), sorrid)
	e=nf90_def_var(f, 'ort', nf90_double, (/rd/), sortid)
	e=nf90_def_var(f, 'orz', nf90_double, (/rd/), sorzid)
	e=nf90_def_var(f, 'ott', nf90_double, (/rd/), sottid)
	e=nf90_def_var(f, 'otz', nf90_double, (/rd/), sotzid)
	e=nf90_def_var(f, 'ozz', nf90_double, (/rd/), sozzid)

	e=nf90_def_var(f, 'rot', nf90_double, (/rd/), srotid)
	e=nf90_def_var(f, 'tor', nf90_double, (/rd/), storid)

	e=nf90_def_var(f, 'dmurr', nf90_double, (/rd/), murrid)   
        e=nf90_def_var(f, 'dmurt', nf90_double, (/rd/), murtid)  
	e=nf90_def_var(f, 'dmurz', nf90_double, (/rd/), murzid)  
	e=nf90_def_var(f, 'dmutr', nf90_double, (/rd/), mutrid)   
        e=nf90_def_var(f, 'dmutt', nf90_double, (/rd/), muttid)  
	e=nf90_def_var(f, 'dmutz', nf90_double, (/rd/), mutzid)  
	e=nf90_def_var(f, 'dmuzr', nf90_double, (/rd/), muzrid)   
        e=nf90_def_var(f, 'dmuzt', nf90_double, (/rd/), muztid)  
	e=nf90_def_var(f, 'dmuzz', nf90_double, (/rd/), muzzid)  

	e=nf90_def_var(f, 'dsurr', nf90_double, (/rd/), surrid)   
        e=nf90_def_var(f, 'dsurt', nf90_double, (/rd/), surtid)  
	e=nf90_def_var(f, 'dsurz', nf90_double, (/rd/), surzid)  
	e=nf90_def_var(f, 'dsutr', nf90_double, (/rd/), sutrid)   
        e=nf90_def_var(f, 'dsutt', nf90_double, (/rd/), suttid)  
	e=nf90_def_var(f, 'dsutz', nf90_double, (/rd/), sutzid)  
	e=nf90_def_var(f, 'dsuzr', nf90_double, (/rd/), suzrid)   
        e=nf90_def_var(f, 'dsuzt', nf90_double, (/rd/), suztid)  
	e=nf90_def_var(f, 'dsuzz', nf90_double, (/rd/), suzzid)  


    	e=nf90_def_var(f, 'P', nf90_double, (/rd/), mpid)   
        e=nf90_def_var(f, 'pp', nf90_double, (/rd/), sppid)
	e=nf90_def_var(f, 'pr', nf90_double, (/rd/), sprid)
	e=nf90_def_var(f, 'pt', nf90_double, (/rd/), sptid)
	e=nf90_def_var(f, 'pz', nf90_double, (/rd/), spzid)
        e=nf90_def_var(f, 'purr', nf90_double, (/rd/), spurrid)
        e=nf90_def_var(f, 'putt', nf90_double, (/rd/), sputtid)
        e=nf90_def_var(f, 'puzz', nf90_double, (/rd/), spuzzid)
	e=nf90_def_var(f, 'purzuzr', nf90_double, (/rd/), spurzuzrid)

	e=nf90_def_var(f, 'dsurruzr', nf90_double, (/rd/), surruzrid)   
        e=nf90_def_var(f, 'dsurtuzt', nf90_double, (/rd/), surtuztid)  
	e=nf90_def_var(f, 'dsurzuzz', nf90_double, (/rd/), surzuzzid)  

    	e=nf90_enddef(f)


    	e=nf90_put_var(f, rid, mes_D%r(1:i_N,1))
        e=nf90_put_var(f, intrdrid, mes_D%intrdr)
        e=nf90_put_var(f, dr0id, mes_D%dr0, start=(/1,1/))
        e=nf90_put_var(f, dr1id, mes_D%dr1, start=(/1,1/))

        e=nf90_put_var(f, mrid, m_r)
	e=nf90_put_var(f, mtid, m_t)
	e=nf90_put_var(f, mzid, m_z)
	e=nf90_put_var(f, mztid, m_zt)
	e=nf90_put_var(f, morid, m_or)
	e=nf90_put_var(f, motid, m_ot)
	e=nf90_put_var(f, mozid, m_oz)
       
	e=nf90_put_var(f, srrid, s_rr)
	e=nf90_put_var(f, srtid, s_rt)
	e=nf90_put_var(f, srzid, s_rz)
	e=nf90_put_var(f, sttid, s_tt)
	e=nf90_put_var(f, stzid, s_tz)
	e=nf90_put_var(f, szzid, s_zz)

	e=nf90_put_var(f, srrrid, s_rrr)
	e=nf90_put_var(f, stttid, s_ttt)
	e=nf90_put_var(f, szzzid, s_zzz)
	e=nf90_put_var(f, srttid, s_rtt)
	e=nf90_put_var(f, srzzid, s_rzz)
	e=nf90_put_var(f, srrzid, s_rrz)
	e=nf90_put_var(f, sttzid, s_ttz)
	e=nf90_put_var(f, sr4id, s_r4)
	e=nf90_put_var(f, st4id, s_t4)
	e=nf90_put_var(f, sz4id, s_z4)
	e=nf90_put_var(f, sr5id, s_r5)
	e=nf90_put_var(f, st5id, s_t5)
	e=nf90_put_var(f, sz5id, s_z5)
	e=nf90_put_var(f, sr6id, s_r6)
	e=nf90_put_var(f, st6id, s_t6)
	e=nf90_put_var(f, sz6id, s_z6)
	e=nf90_put_var(f, sr7id, s_r7)
	e=nf90_put_var(f, st7id, s_t7)
	e=nf90_put_var(f, sz7id, s_z7)
	e=nf90_put_var(f, sr8id, s_r8)
	e=nf90_put_var(f, st8id, s_t8)
	e=nf90_put_var(f, sz8id, s_z8)
	e=nf90_put_var(f, sr9id, s_r9)
	e=nf90_put_var(f, st9id, s_t9)
	e=nf90_put_var(f, sz9id, s_z9)
	e=nf90_put_var(f, sr10id, s_r10)
	e=nf90_put_var(f, st10id, s_t10)
	e=nf90_put_var(f, sz10id, s_z10)




	e=nf90_put_var(f, sorrid, s_orr)
	e=nf90_put_var(f, sortid, s_ort)
	e=nf90_put_var(f, sorzid, s_orz)
	e=nf90_put_var(f, sottid, s_ott)
	e=nf90_put_var(f, sotzid, s_otz)
	e=nf90_put_var(f, sozzid, s_ozz)
	e=nf90_put_var(f, srotid, s_rot)
	e=nf90_put_var(f, storid, s_tor)

	e=nf90_put_var(f, murrid, m_urr) 
	e=nf90_put_var(f, murtid, m_urt) 
	e=nf90_put_var(f, murzid, m_urz) 
	e=nf90_put_var(f, mutrid, m_utr) 
	e=nf90_put_var(f, muttid, m_utt) 
	e=nf90_put_var(f, mutzid, m_utz) 
	e=nf90_put_var(f, muzrid, m_uzr) 
	e=nf90_put_var(f, muztid, m_uzt) 
	e=nf90_put_var(f, muzzid, m_uzz) 

	e=nf90_put_var(f, surrid, s_urr) 
	e=nf90_put_var(f, surtid, s_urt) 
	e=nf90_put_var(f, surzid, s_urz) 
	e=nf90_put_var(f, sutrid, s_utr) 
	e=nf90_put_var(f, suttid, s_utt) 
	e=nf90_put_var(f, sutzid, s_utz) 
	e=nf90_put_var(f, suzrid, s_uzr) 
	e=nf90_put_var(f, suztid, s_uzt) 
	e=nf90_put_var(f, suzzid, s_uzz) 

	e=nf90_put_var(f, mpid, m_p)    
	e=nf90_put_var(f, sppid, s_pp)
	e=nf90_put_var(f, sprid, s_pr)
	e=nf90_put_var(f, sptid, s_pt)
	e=nf90_put_var(f, spzid, s_pz)
        e=nf90_put_var(f, spurrid, s_purr)
	e=nf90_put_var(f, sputtid, s_putt)
	e=nf90_put_var(f, spuzzid, s_puzz)
	e=nf90_put_var(f, spurzuzrid, s_purzuzr)

	e=nf90_put_var(f, surruzrid, s_urruzr) 
	e=nf90_put_var(f, surtuztid, s_urtuzt) 
	e=nf90_put_var(f, surzuzzid, s_urzuzz) 

     

        e=nf90_close(f)
    end if
   
   end subroutine io_write_velocity_statistics
   
   
!--------------------------------------------------------------------------
!  subroutine to do energy spectrum statistics
!--------------------------------------------------------------------------
   subroutine io_write_energy_statistics()
   use timestep
   implicit none

   double precision :: Ekz(1:i_N,-i_K1:i_K1),Ekt(1:i_N,0:i_M1),Ekz_tmp(1:i_N,-i_K1:i_K1),Ekt_tmp(1:i_N,0:i_M1)
   double precision :: d1, d(i_N),dfact

   character(8) :: cnum
   integer :: e, f, kd,md,rd,rid,intrdrid,iKL1d,iKL2d,iKL3d,drid,dr0id,dr1id
   integer :: Ekzrr,Ektrr,Ekztt,Ekttt,Ekzzz,Ektzz,Ekzrz,Ektrz,Ekzorr,Ektorr,Ekzott,Ektott,Ekzozz,Ektozz
   integer :: n,n_, ff,fl,fi

   _loop_km_vars

   
   write(cnum,'(I8.8)') tim_step
   if(mpi_rnk==0) then
         print*, ' saving statistics'//cnum//'  t=', tim_t

    	e=nf90_create('Eenergy'//cnum//'.cdf.dat', NF90_64BIT_OFFSET, f)
    	e=nf90_put_att(f, nf90_global, 'Re', d_Re)
    	e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
    	e=nf90_put_att(f, nf90_global, 't', tim_t)
    	e=nf90_def_dim(f, 'r', i_N, rd)
    	e=nf90_def_dim(f, 'K', 2*i_K-1, Kd)
    	e=nf90_def_dim(f, 'M', i_M, Md)
        e=nf90_def_dim(f, 'iKL1', 2*i_KL+1, iKL1d)
        e=nf90_def_dim(f, 'iKL2', i_N+i_KL, iKL2d)
        e=nf90_def_dim(f, 'iKL3', i_KL+1, iKL3d)
    
    	e=nf90_def_var(f, 'r', nf90_double, (/rd/), rid)
        e=nf90_def_var(f, 'intrdr', nf90_double, (/rd/), intrdrid)
        e=nf90_def_var(f, 'dr0', nf90_double, (/iKL3d,iKL3d/), dr0id)
        e=nf90_def_var(f, 'dr1', nf90_double, (/iKL3d,iKL3d/), dr1id)


	
        e=nf90_def_var(f, 'Ekz_rr', nf90_double, (/rd,Kd/), Ekzrr)
    	e=nf90_def_var(f, 'Ekt_rr', nf90_double, (/rd,Md/), Ektrr)
	e=nf90_def_var(f, 'Ekz_tt', nf90_double, (/rd,Kd/), Ekztt)
    	e=nf90_def_var(f, 'Ekt_tt', nf90_double, (/rd,Md/), Ekttt)
	e=nf90_def_var(f, 'Ekz_zz', nf90_double, (/rd,Kd/), Ekzzz)
    	e=nf90_def_var(f, 'Ekt_zz', nf90_double, (/rd,Md/), Ektzz)
	e=nf90_def_var(f, 'Ekz_rz', nf90_double, (/rd,Kd/), Ekzrz)
    	e=nf90_def_var(f, 'Ekt_rz', nf90_double, (/rd,Md/), Ektrz)
	e=nf90_def_var(f, 'Ekz_orr', nf90_double, (/rd,Kd/), Ekzorr)
    	e=nf90_def_var(f, 'Ekt_orr', nf90_double, (/rd,Md/), Ektorr)
	e=nf90_def_var(f, 'Ekz_ott', nf90_double, (/rd,Kd/), Ekzott)
    	e=nf90_def_var(f, 'Ekt_ott', nf90_double, (/rd,Md/), Ektott)
	e=nf90_def_var(f, 'Ekz_ozz', nf90_double, (/rd,Kd/), Ekzozz)
    	e=nf90_def_var(f, 'Ekt_ozz', nf90_double, (/rd,Md/), Ektozz)

    	e=nf90_enddef(f)

    	e=nf90_put_var(f, rid, mes_D%r(1:i_N,1))
        e=nf90_put_var(f, intrdrid, mes_D%intrdr)
        e=nf90_put_var(f, dr0id, mes_D%dr0, start=(/1,1/))
        e=nf90_put_var(f, dr1id, mes_D%dr1, start=(/1,1/))
                
    end if
!! ur 
    Ekz_tmp=0d0; Ekt_tmp=0d0;
    _loop_km_begin
	dfact=1d0
	if (m.eq.0) dfact=2d0
	Ekz_tmp(:,k)=Ekz_tmp(:,k)+dfact*(vel_ur%Re(:,nh)*vel_ur%Re(:,nh)+&
		vel_ur%Im(:,nh)*vel_ur%Im(:,nh))
	Ekt_tmp(:,m)=Ekt_tmp(:,m)+dfact*(vel_ur%Re(:,nh)*vel_ur%Re(:,nh)+&
		vel_ur%Im(:,nh)*vel_ur%Im(:,nh))
    _loop_km_end
      
#ifdef _MPI
    call mpi_allreduce(Ekz_tmp, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    call mpi_allreduce(Ekt_tmp, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
#endif
    d1 = dble( i_M)
    Ekz= Ekz/d1
    d1 = dble( (2*i_K-1))
    Ekt= Ekt/d1
    if(mpi_rnk==0) then
    	e=nf90_put_var(f,Ekzrr,Ekz, start=(/1,1/))
        e=nf90_put_var(f,Ektrr,Ekt, start=(/1,1/))
    end if
    
!! ut 
    Ekz_tmp=0d0; Ekt_tmp=0d0;
    _loop_km_begin
	dfact=1d0
	if (m.eq.0) dfact=2d0
	Ekz_tmp(:,k)=Ekz_tmp(:,k)+dfact*(vel_ut%Re(:,nh)*vel_ut%Re(:,nh)+&
		vel_ut%Im(:,nh)*vel_ut%Im(:,nh))
	Ekt_tmp(:,m)=Ekt_tmp(:,m)+dfact*(vel_ut%Re(:,nh)*vel_ut%Re(:,nh)+&
		vel_ut%Im(:,nh)*vel_ut%Im(:,nh))
    _loop_km_end
      
#ifdef _MPI
    call mpi_allreduce(Ekz_tmp, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    call mpi_allreduce(Ekt_tmp, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
#endif
    d1 = dble( i_M)
    Ekz= Ekz/d1
    d1 = dble( (2*i_K-1))
    Ekt= Ekt/d1
    if(mpi_rnk==0) then
	e=nf90_put_var(f,Ekztt,Ekz, start=(/1,1/))
        e=nf90_put_var(f,Ekttt,Ekt, start=(/1,1/))
    end if
    
!! uz
    Ekz_tmp=0d0; Ekt_tmp=0d0;
    _loop_km_begin
	dfact=1d0
	if (m.eq.0) dfact=2d0
	Ekz_tmp(:,k)=Ekz_tmp(:,k)+dfact*(vel_uz%Re(:,nh)*vel_uz%Re(:,nh)+&
		vel_uz%Im(:,nh)*vel_uz%Im(:,nh))
	Ekt_tmp(:,m)=Ekt_tmp(:,m)+dfact*(vel_uz%Re(:,nh)*vel_uz%Re(:,nh)+&
		vel_uz%Im(:,nh)*vel_uz%Im(:,nh))
    _loop_km_end
      
#ifdef _MPI
    call mpi_allreduce(Ekz_tmp, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    call mpi_allreduce(Ekt_tmp, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
#endif
    d1 = dble( i_M)
    Ekz= Ekz/d1
    d1 = dble( (2*i_K-1))
    Ekt= Ekt/d1
    if(mpi_rnk==0) then
	e=nf90_put_var(f,Ekzzz,Ekz, start=(/1,1/))
        e=nf90_put_var(f,Ektzz,Ekt, start=(/1,1/))
    end if
    
!! uruz
    Ekz_tmp=0d0; Ekt_tmp=0d0;
    _loop_km_begin
	dfact=1d0
	if (m.eq.0) dfact=2d0
	Ekz_tmp(:,k)=Ekz_tmp(:,k)+dfact*(vel_ur%Re(:,nh)*vel_uz%Re(:,nh)+&
		vel_ur%Im(:,nh)*vel_uz%Im(:,nh))
	Ekt_tmp(:,m)=Ekt_tmp(:,m)+dfact*(vel_ur%Re(:,nh)*vel_uz%Re(:,nh)+&
		vel_ur%Im(:,nh)*vel_uz%Im(:,nh))
    _loop_km_end
      
#ifdef _MPI
    call mpi_allreduce(Ekz_tmp, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    call mpi_allreduce(Ekt_tmp, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
#endif
    d1 = dble( i_M)
    Ekz= Ekz/d1
    d1 = dble( (2*i_K-1))
    Ekt= Ekt/d1
    if(mpi_rnk==0) then
	e=nf90_put_var(f,Ekzrz,Ekz, start=(/1,1/))
        e=nf90_put_var(f,Ektrz,Ekt, start=(/1,1/))
    end if
    
         
!! Vorticity Statistis Spectrum
      call var_coll_curl(vel_ur,vel_ut,vel_uz, c1,c2,c3)
      	
!! or
    Ekz_tmp=0d0; Ekt_tmp=0d0;
    _loop_km_begin
	dfact=1d0
	if (m.eq.0) dfact=2d0
	Ekz_tmp(:,k)=Ekz_tmp(:,k)+dfact*(c1%Re(:,nh)*c1%Re(:,nh)+&
		c1%Im(:,nh)*c1%Im(:,nh))
	Ekt_tmp(:,m)=Ekt_tmp(:,m)+dfact*(c1%Re(:,nh)*c1%Re(:,nh)+&
		c1%Im(:,nh)*c1%Im(:,nh))
    _loop_km_end
      
#ifdef _MPI
    call mpi_allreduce(Ekz_tmp, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    call mpi_allreduce(Ekt_tmp, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
#endif
    d1 = dble( i_M)
    Ekz= Ekz/d1
    d1 = dble( (2*i_K-1))
    Ekt= Ekt/d1
    if(mpi_rnk==0) then
	e=nf90_put_var(f,Ekzorr,Ekz, start=(/1,1/))
        e=nf90_put_var(f,Ektorr,Ekt, start=(/1,1/))
    end if
    	

!! ot
    Ekz_tmp=0d0; Ekt_tmp=0d0;
    _loop_km_begin
	dfact=1d0
	if (m.eq.0) dfact=2d0
	Ekz_tmp(:,k)=Ekz_tmp(:,k)+dfact*(c2%Re(:,nh)*c2%Re(:,nh)+&
		c2%Im(:,nh)*c2%Im(:,nh))
	Ekt_tmp(:,m)=Ekt_tmp(:,m)+dfact*(c2%Re(:,nh)*c2%Re(:,nh)+&
		c2%Im(:,nh)*c2%Im(:,nh))
    _loop_km_end
      
#ifdef _MPI
    call mpi_allreduce(Ekz_tmp, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    call mpi_allreduce(Ekt_tmp, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
#endif
    d1 = dble( i_M)
    Ekz= Ekz/d1
    d1 = dble( (2*i_K-1))
    Ekt= Ekt/d1
    if(mpi_rnk==0) then
	e=nf90_put_var(f,Ekzott,Ekz, start=(/1,1/))
        e=nf90_put_var(f,Ektott,Ekt, start=(/1,1/))
    end if
    
!! oz
    Ekz_tmp=0d0; Ekt_tmp=0d0;
    _loop_km_begin
	dfact=1d0
	if (m.eq.0) dfact=2d0	
	Ekz_tmp(:,k)=Ekz_tmp(:,k)+dfact*(c3%Re(:,nh)*c3%Re(:,nh)+&
		c3%Im(:,nh)*c3%Im(:,nh))
	Ekt_tmp(:,m)=Ekt_tmp(:,m)+dfact*(c3%Re(:,nh)*c3%Re(:,nh)+&
		c3%Im(:,nh)*c3%Im(:,nh))
    _loop_km_end
      
#ifdef _MPI
    call mpi_allreduce(Ekz_tmp, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    call mpi_allreduce(Ekt_tmp, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
#endif
    d1 = dble( i_M)
    Ekz= Ekz/d1
    d1 = dble( (2*i_K-1))
    Ekt= Ekt/d1
    if(mpi_rnk==0) then
	e=nf90_put_var(f,Ekzozz,Ekz, start=(/1,1/))
        e=nf90_put_var(f,Ektozz,Ekt, start=(/1,1/))
        e=nf90_close(f)
    end if
    
     
   end subroutine io_write_energy_statistics   


   !--------------------------------------------------------------------------
!  subroutine to calculate pressure
!--------------------------------------------------------------------------
   subroutine io_pressure()
   use timestep
   implicit none

   double precision :: BCR(0:i_pH1), BCI(0:i_pH1)
   character(8) :: cnum


   integer :: n,n_, ff,fl,fi
   _loop_km_vars

   write(cnum,'(I8.8)') tim_step
   if(mpi_rnk==0) then
         print*, ' calculating pressure'//cnum//'  t=', tim_t
   end if



      !Laplacian of U
      call var_coll_curl(vel_ur,vel_ut,vel_uz, c1,c2,c3)
      call var_coll_curl(c1,c2,c3, c1,c2,c3)
      _loop_km_begin
        BCR(nh) = - c1%Re(i_N,nh)/d_Re
        BCI(nh) = - c1%Im(i_N,nh)/d_Re
      _loop_km_end
      ! r equation
      !durdr
      call var_coll_meshmult(1,mes_D%dr(1),vel_ur,c1)
      call tra_coll2phys(c1,p1)
      pr%Re=vel_r%Re*p1%Re

      ! 1/r(durdt-ut)
      _loop_km_begin
            c1%Re(:,nh) = mes_D%r(:,-1)*(-vel_ur%Im(:,nh)*m*i_Mp-vel_ut%Re(:,nh))
            c1%Im(:,nh) = mes_D%r(:,-1)*( vel_ur%Re(:,nh)*m*i_Mp-vel_ut%Im(:,nh))
      _loop_km_end
      call tra_coll2phys(c1,p1)
      pr%Re=pr%Re+vel_t%Re*p1%Re

      ! durdz
      _loop_km_begin
            c1%Re(:,nh) = -vel_ur%Im(:,nh)*d_alpha*k
            c1%Im(:,nh) =  vel_ur%Re(:,nh)*d_alpha*k
      _loop_km_end
      call tra_coll2phys(c1,p1)
      pr%Re=pr%Re+vel_z%Re*p1%Re
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
	 pr%Re(:,:,n)=pr%Re(:,:,n)+p1%Re(:,:,n)*vel_U(n_)
      end do
      call tra_phys2coll(pr,c1)

      ! theta equation
      ! dutdr
      call var_coll_meshmult(1,mes_D%dr(1),vel_ut,c2)
      call tra_coll2phys(c2,p1)
      pr%Re=vel_r%Re*p1%Re

      ! 1/r(dutdt+ur)
      _loop_km_begin
            c2%Re(:,nh) = mes_D%r(:,-1)*(-vel_ut%Im(:,nh)*m*i_Mp+vel_ur%Re(:,nh))
            c2%Im(:,nh) = mes_D%r(:,-1)*( vel_ut%Re(:,nh)*m*i_Mp+vel_ur%Im(:,nh))
      _loop_km_end
      call tra_coll2phys(c2,p1)
      pr%Re=pr%Re+vel_t%Re*p1%Re

      ! dutdz
      _loop_km_begin
            c2%Re(:,nh) = -vel_ut%Im(:,nh)*d_alpha*k
            c2%Im(:,nh) =  vel_ut%Re(:,nh)*d_alpha*k
      _loop_km_end
      call tra_coll2phys(c2,p1)
      pr%Re=pr%Re+vel_z%Re*p1%Re
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
         pr%Re(:,:,n)=pr%Re(:,:,n)+p1%Re(:,:,n)*vel_U(n_)
      end do
      call tra_phys2coll(pr,c2)

      ! z equation
      ! duzdr
      call var_coll_meshmult(0,mes_D%dr(1),vel_uz,c3)
      call tra_coll2phys(c3,p1)
      pr%Re=vel_r%Re*p1%Re

      ! 1/r(duzdt)
      _loop_km_begin
            c3%Re(:,nh) = mes_D%r(:,-1)*(-vel_uz%Im(:,nh)*m*i_Mp)
            c3%Im(:,nh) = mes_D%r(:,-1)*( vel_uz%Re(:,nh)*m*i_Mp)
      _loop_km_end
      call tra_coll2phys(c3,p1)
      pr%Re=pr%Re+vel_t%Re*p1%Re

      ! duzdz
      _loop_km_begin
            c3%Re(:,nh) = -vel_uz%Im(:,nh)*d_alpha*k
            c3%Im(:,nh) =  vel_uz%Re(:,nh)*d_alpha*k
      _loop_km_end
      call tra_coll2phys(c3,p1)
      pr%Re=pr%Re+vel_z%Re*p1%Re
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
         pr%Re(:,:,n)=pr%Re(:,:,n)+p1%Re(:,:,n)*vel_U(n_)+vel_r%Re(:,:,n)*vel_Up(n_)
      end do
      call tra_phys2coll(pr,c3)

      _loop_km_begin
        BCR(nh) = BCR(nh) - c1%Re(i_N,nh)
        BCI(nh) = BCI(nh) - c1%Im(i_N,nh)
      _loop_km_end
      call var_coll_div(c1,c2,c3, c1)
      c1%Re=-c1%Re;
      c1%Im=-c1%Im;
      _loop_km_begin
        c1%Re(i_N,nh) = BCR(nh)
        c1%Im(i_N,nh) = BCI(nh)
      _loop_km_end

      !call tim_zerobc(c1)
      call tim_lumesh_inits( 0,1,0d0,1d0, LNps)
      call tim_lumesh_invert(0,LNps, c1)
      call tra_coll2phys(c1,pr)

   end subroutine io_pressure
  
!**************************************************************************
 end module io
!**************************************************************************

