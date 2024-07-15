#include "../parallel.h"
!*************************************************************************
 PROGRAM LOADSTATE
!*************************************************************************
   use io
   implicit none
  
   double precision :: Ekz_rr(1:i_N,-i_K1:i_K1),Ekt_rr(1:i_N,0:i_M1),Ekz(1:i_N,-i_K1:i_K1),Ekt(1:i_N,0:i_M1)
   double precision :: Ekz_tt(1:i_N,-i_K1:i_K1),Ekt_tt(1:i_N,0:i_M1),Ekz_zz(1:i_N,-i_K1:i_K1),Ekt_zz(1:i_N,0:i_M1)
   double precision :: Ekz_rz(1:i_N,-i_K1:i_K1),Ekt_rz(1:i_N,0:i_M1),Ekz_orr(1:i_N,-i_K1:i_K1),Ekt_orr(1:i_N,0:i_M1)
   double precision :: Ekz_ott(1:i_N,-i_K1:i_K1),Ekt_ott(1:i_N,0:i_M1),Ekz_ozz(1:i_N,-i_K1:i_K1),Ekt_ozz(1:i_N,0:i_M1)

   double precision :: d1,dfact
   character(4) :: cnum
   integer :: e, f,r,kd,md,rd,Ekzrr,Ektrr,Ekztt,Ekttt,Ekzzz,Ektzz,Ekzrz,Ektrz,Ekzorr,Ektorr,Ekzott,Ektott,Ekzozz,Ektozz
   integer :: n,n_, ff,fl,fi
   _loop_km_vars
 9 format(A)

   print*, 'initialising...'
   call mpi_precompute()
   call par_precompute()
   call mes_precompute()
   call var_precompute()
   call tra_precompute()
   call tim_precompute()
   call vel_precompute()
   call  io_precompute()
   
   if(mpi_sze==1) then
      print*, 'Enter first and last file numbers'
      print*, ' (current dir; approx equally spaced in time):'
      read(*,*)  ff, fl
   else
      open(99,status='old',file='MEANFLUCT.IN')
      read(99,*) ff, fl
      close(99)
   end if

   

   Ekz_rr=0d0;  Ekt_rr=0d0; Ekz_tt=0d0;  Ekt_tt=0d0; Ekz_zz=0d0;  Ekt_zz=0d0;
   Ekz_rz=0d0;  Ekt_rz=0d0; Ekz_orr=0d0;  Ekt_orr=0d0; Ekz_ott=0d0;  Ekt_ott=0d0;
   Ekz_ozz=0d0;  Ekt_ozz=0d0

   do fi = ff, fl
      write(cnum,'(I4.4)') fi
      io_statefile = 'state'//cnum//'.cdf.dat'
      call io_load_state()
      call var_coll_curl(vel_ur,vel_ut,vel_uz, vel_Nr,vel_Nt,vel_Nz)
   
      _loop_km_begin
	dfact=1d0
	if (m.eq.0) dfact=2d0

	Ekz_rr(:,k)=Ekz_rr(:,k)+dfact*(vel_ur%Re(:,nh)*vel_ur%Re(:,nh)+&
		vel_ur%Im(:,nh)*vel_ur%Im(:,nh))
	Ekt_rr(:,m)=Ekt_rr(:,m)+dfact*(vel_ur%Re(:,nh)*vel_ur%Re(:,nh)+&
		vel_ur%Im(:,nh)*vel_ur%Im(:,nh))
	
	Ekz_tt(:,k)=Ekz_tt(:,k)+dfact*(vel_ut%Re(:,nh)*vel_ut%Re(:,nh)+&
		vel_ut%Im(:,nh)*vel_ut%Im(:,nh))
	Ekt_tt(:,m)=Ekt_tt(:,m)+dfact*(vel_ut%Re(:,nh)*vel_ut%Re(:,nh)+&
		vel_ut%Im(:,nh)*vel_ut%Im(:,nh))

	Ekz_zz(:,k)=Ekz_zz(:,k)+dfact*(vel_uz%Re(:,nh)*vel_uz%Re(:,nh)+&
		vel_uz%Im(:,nh)*vel_uz%Im(:,nh))
	Ekt_zz(:,m)=Ekt_zz(:,m)+dfact*(vel_uz%Re(:,nh)*vel_uz%Re(:,nh)+&
		vel_uz%Im(:,nh)*vel_uz%Im(:,nh))

	Ekz_rz(:,k)=Ekz_rz(:,k)+dfact*(vel_ur%Re(:,nh)*vel_uz%Re(:,nh)+&
		vel_ur%Im(:,nh)*vel_uz%Im(:,nh))
	Ekt_rz(:,m)=Ekt_rz(:,m)+dfact*(vel_ur%Re(:,nh)*vel_uz%Re(:,nh)+&
		vel_ur%Im(:,nh)*vel_uz%Im(:,nh))

	Ekz_orr(:,k)=Ekz_orr(:,k)+dfact*(vel_Nr%Re(:,nh)*vel_Nr%Re(:,nh)+&
		vel_Nr%Im(:,nh)*vel_Nr%Im(:,nh))
	Ekt_orr(:,m)=Ekt_orr(:,m)+dfact*(vel_Nr%Re(:,nh)*vel_Nr%Re(:,nh)+&
		vel_Nr%Im(:,nh)*vel_Nr%Im(:,nh))

	Ekz_ott(:,k)=Ekz_ott(:,k)+dfact*(vel_Nt%Re(:,nh)*vel_Nt%Re(:,nh)+&
		vel_Nt%Im(:,nh)*vel_Nt%Im(:,nh))
	Ekt_ott(:,m)=Ekt_ott(:,m)+dfact*(vel_Nt%Re(:,nh)*vel_Nt%Re(:,nh)+&
		vel_Nt%Im(:,nh)*vel_Nt%Im(:,nh))

	Ekz_ozz(:,k)=Ekz_ozz(:,k)+dfact*(vel_Nz%Re(:,nh)*vel_Nz%Re(:,nh)+&
		vel_Nz%Im(:,nh)*vel_Nz%Im(:,nh))
	Ekt_ozz(:,m)=Ekt_ozz(:,m)+dfact*(vel_Nz%Re(:,nh)*vel_Nz%Re(:,nh)+&
		vel_Nz%Im(:,nh)*vel_Nz%Im(:,nh))

      _loop_km_end


   end do

#ifdef _MPI
    call mpi_allreduce(Ekz_rr, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    Ekz_rr= Ekz
    call mpi_allreduce(Ekt_rr, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
    Ekt_rr= Ekt

    call mpi_allreduce(Ekz_tt, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    Ekz_tt= Ekz
    call mpi_allreduce(Ekt_tt, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
    Ekt_tt= Ekt

    call mpi_allreduce(Ekz_zz, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    Ekz_zz= Ekz
    call mpi_allreduce(Ekt_zz, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
    Ekt_zz= Ekt

    call mpi_allreduce(Ekz_rz, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    Ekz_rz= Ekz
    call mpi_allreduce(Ekt_rz, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
    Ekt_rz= Ekt

    call mpi_allreduce(Ekz_orr, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    Ekz_orr= Ekz
    call mpi_allreduce(Ekt_orr, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
    Ekt_orr= Ekt

    call mpi_allreduce(Ekz_ott, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    Ekz_ott= Ekz
    call mpi_allreduce(Ekt_ott, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
    Ekt_ott= Ekt

    call mpi_allreduce(Ekz_ozz, Ekz, i_N*(2*i_K-1), mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er) 
    Ekz_ozz= Ekz
    call mpi_allreduce(Ekt_ozz, Ekt, i_N*i_M, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
    Ekt_ozz= Ekt
#endif

    d1 = dble((fl-ff+1) * i_M)
    Ekz_rr= Ekz_rr/d1
    Ekz_tt= Ekz_tt/d1
    Ekz_zz= Ekz_zz/d1
    Ekz_rz= Ekz_rz/d1
    Ekz_orr= Ekz_orr/d1
    Ekz_ott= Ekz_ott/d1
    Ekz_ozz= Ekz_ozz/d1
   


    d1 = dble((fl-ff+1) * (2*i_K-1))
    Ekt_rr= Ekt_rr/d1
    Ekt_tt= Ekt_tt/d1
    Ekt_zz= Ekt_zz/d1
    Ekt_rz= Ekt_rz/d1
    Ekt_orr= Ekt_orr/d1
    Ekt_ott= Ekt_ott/d1
    Ekt_ozz= Ekt_ozz/d1

    if(mpi_rnk==0) then
    	print*, ' saving spectrum'
    
    	e=nf90_create('Energy.cdf.dat', nf90_clobber, f)
    	e=nf90_put_att(f, nf90_global, 'Re', d_Re)
    	e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
    	e=nf90_def_dim(f, 'r', i_N, rd)
    	e=nf90_def_dim(f, 'K', 2*i_K-1, Kd)
    	e=nf90_def_dim(f, 'M', i_M, Md)
    
    	e=nf90_def_var(f, 'r',     nf90_double, (/rd/), r)
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


    	e=nf90_put_var(f, r, mes_D%r(1:i_N,1))
        e=nf90_put_var(f,Ekzrr,Ekz_rr, start=(/1,1/))
        e=nf90_put_var(f,Ektrr,Ekt_rr, start=(/1,1/))
	e=nf90_put_var(f,Ekztt,Ekz_tt, start=(/1,1/))
        e=nf90_put_var(f,Ekttt,Ekt_tt, start=(/1,1/))
	e=nf90_put_var(f,Ekzzz,Ekz_zz, start=(/1,1/))
        e=nf90_put_var(f,Ektzz,Ekt_zz, start=(/1,1/))
	e=nf90_put_var(f,Ekzrz,Ekz_rz, start=(/1,1/))
        e=nf90_put_var(f,Ektrz,Ekt_rz, start=(/1,1/))
	e=nf90_put_var(f,Ekzorr,Ekz_orr, start=(/1,1/))
        e=nf90_put_var(f,Ektorr,Ekt_orr, start=(/1,1/))
	e=nf90_put_var(f,Ekzott,Ekz_ott, start=(/1,1/))
        e=nf90_put_var(f,Ektott,Ekt_ott, start=(/1,1/))
	e=nf90_put_var(f,Ekzozz,Ekz_ozz, start=(/1,1/))
        e=nf90_put_var(f,Ektozz,Ekt_ozz, start=(/1,1/))

        e=nf90_close(f)
    end if

#ifdef _MPI
   call mpi_barrier(mpi_comm_world, mpi_er)
   call mpi_finalize(mpi_er)
#endif

!*************************************************************************
 END PROGRAM LOADSTATE
!*************************************************************************

