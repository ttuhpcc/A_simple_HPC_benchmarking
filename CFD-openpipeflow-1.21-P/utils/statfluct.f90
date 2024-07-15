#include "../parallel.h"
!*************************************************************************
 PROGRAM LOADSTATE
!*************************************************************************
   use io
   implicit none
   double precision :: M_u(i_N), M_v(i_N),M_w(i_N),M_ox(i_N), M_oy(i_N),M_oz(i_N)
   double precision :: S_uu(i_N),S_uv(i_N),S_uw(i_N),S_vv(i_N),S_vw(i_N),S_ww(i_N)
   double precision :: S_uuu(i_N),S_vvv(i_N),S_www(i_N),S_u4(i_N),S_v4(i_N),S_w4(i_N)
   double precision :: S_oxx(i_N),S_oxy(i_N),S_oxz(i_N),S_oyy(i_N),S_oyz(i_N),S_ozz(i_N)
   double precision :: d1, d(i_N)
   double precision :: Ekz_uu(1:i_N,-i_K1:i_K1),Ekt_uu(1:i_N,0:i_M1),Ekz(1:i_N,-i_K1:i_K1),Ekt(1:i_N,0:i_M1)
   character(4) :: cnum
   integer :: e, f, kd,md,rd,rid,muid,mvid,mwid,moxid,moyid,mozid,suuid,suvid,suwid,svvid,svwid,swwid
   integer :: suuuid,svvvid,swwwid,su4id,sv4id,sw4id,soxxid,soxyid,soxzid,soyyid,soyzid,sozzid
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

   M_u = 0d0; M_v = 0d0; M_w = 0d0; M_ox = 0d0;M_oy = 0d0; M_oz = 0d0
   S_uu= 0d0; S_uv= 0d0; S_uw= 0d0; S_vv= 0d0; S_vw= 0d0; S_ww= 0d0;
   S_uuu= 0d0; S_vvv= 0d0; S_www= 0d0; S_u4= 0d0; S_v4= 0d0; S_w4= 0d0;
   S_oxx= 0d0; S_oxy= 0d0; S_oxz= 0d0; S_oyy= 0d0; S_oyz= 0d0; S_ozz= 0d0;

   Ekz_uu=0d0;  Ekt_uu=0d0;
   do fi = ff, fl
      write(cnum,'(I4.4)') fi
      io_statefile = 'state'//cnum//'.cdf.dat'
      call io_load_state()
      call vel_transform()
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
         M_u(n_) = M_u(n_) + sum(vel_r%Re(:,:,n))
         M_v(n_) = M_v(n_) + sum(vel_t%Re(:,:,n))
         M_w(n_) = M_w(n_) + sum(vel_z%Re(:,:,n))
	 M_ox(n_) = M_ox(n_) + sum(vel_curlr%Re(:,:,n))
         M_oy(n_) = M_oy(n_) + sum(vel_curlt%Re(:,:,n))
	 M_oz(n_) = M_oz(n_) + sum(vel_curlz%Re(:,:,n))

         s_uu(n_) = s_uu(n_) + sum(vel_r%Re(:,:,n)**2)
         s_vv(n_) = s_vv(n_) + sum(vel_t%Re(:,:,n)**2)
	 s_ww(n_) = s_ww(n_) + sum(vel_z%Re(:,:,n)**2)
         s_uv(n_) = s_uv(n_) + sum(vel_r%Re(:,:,n)*vel_t%Re(:,:,n))
	 s_uw(n_) = s_uw(n_) + sum(vel_r%Re(:,:,n)*vel_z%Re(:,:,n))
         s_vw(n_) = s_vw(n_) + sum(vel_t%Re(:,:,n)*vel_z%Re(:,:,n))

	 s_uuu(n_) = s_uuu(n_) + sum(vel_r%Re(:,:,n)**3)
         s_vvv(n_) = s_vvv(n_) + sum(vel_t%Re(:,:,n)**3)
	 s_www(n_) = s_www(n_) + sum(vel_z%Re(:,:,n)**3)

	 s_u4(n_) = s_u4(n_) + sum(vel_r%Re(:,:,n)**4)
         s_v4(n_) = s_v4(n_) + sum(vel_t%Re(:,:,n)**4)
	 s_w4(n_) = s_w4(n_) + sum(vel_z%Re(:,:,n)**4)

	 s_oxx(n_) = s_oxx(n_) + sum(vel_curlr%Re(:,:,n)**2)
	 s_oyy(n_) = s_oyy(n_) + sum(vel_curlt%Re(:,:,n)**2)
	 s_ozz(n_) = s_ozz(n_) + sum(vel_curlz%Re(:,:,n)**2)
	 s_oxy(n_) = s_oxy(n_) + sum(vel_curlr%Re(:,:,n)*vel_curlt%Re(:,:,n))
	 s_oxz(n_) = s_oxz(n_) + sum(vel_curlr%Re(:,:,n)*vel_curlz%Re(:,:,n))
	 s_oyz(n_) = s_oyz(n_) + sum(vel_curlt%Re(:,:,n)*vel_curlz%Re(:,:,n))
      end do


      
     

   end do

#ifdef _MPI
   call mpi_allreduce(m_u, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_u = d
   call mpi_allreduce(m_v, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_v = d
   call mpi_allreduce(m_w, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_w = d
   call mpi_allreduce(m_ox, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_ox = d
   call mpi_allreduce(m_oy, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_oy = d
   call mpi_allreduce(m_oz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   m_oz = d

   call mpi_allreduce(s_uu, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_uu = d
   call mpi_allreduce(s_uv, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_uv = d
   call mpi_allreduce(s_uw, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_uw = d
   call mpi_allreduce(s_vv, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_vv = d
   call mpi_allreduce(s_vw, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_vw = d
   call mpi_allreduce(s_ww, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_ww = d

   call mpi_allreduce(s_uuu, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_uuu = d
   call mpi_allreduce(s_vvv, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_vvv = d
   call mpi_allreduce(s_www, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_www = d
   call mpi_allreduce(s_u4, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_u4 = d
   call mpi_allreduce(s_v4, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_v4 = d
   call mpi_allreduce(s_w4, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_w4 = d


   call mpi_allreduce(s_oxx, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_oxx = d
   call mpi_allreduce(s_oxy, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_oxy = d
   call mpi_allreduce(s_oxz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_oxz = d
   call mpi_allreduce(s_oyy, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_oyy = d
   call mpi_allreduce(s_oyz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_oyz = d
   call mpi_allreduce(s_ozz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_ozz = d
#endif

   d1 = dble((fl-ff+1) * i_Th*i_Z)
   m_u = m_u / d1
   m_v = m_v / d1
   m_w = m_w / d1
   m_ox = m_ox / d1
   m_oy = m_oy / d1
   m_oz = m_oz / d1
   

   s_uu = ( s_uu/d1 - m_u**2 )
   s_uv = ( s_uv/d1 - m_u*m_v)
   s_uw = ( s_uw/d1 - m_u*m_w)
   s_vv = ( s_vv/d1 - m_v**2 )
   s_vw = ( s_vw/d1 - m_v*m_w)
   s_ww = ( s_ww/d1 - m_w**2 )

   s_uuu = ( s_uuu/d1 - m_u**3)
   s_vvv = ( s_vvv/d1 - m_v**3)
   s_www = ( s_www/d1 - m_w**3)
   s_u4 = ( s_u4/d1 - m_u**4)
   s_v4 = ( s_v4/d1 - m_v**4)
   s_w4 = ( s_w4/d1 - m_w**4)

  
   s_oxx = ( s_oxx/d1 - m_ox**2 )
   s_oxy = ( s_oxy/d1 - m_ox*m_oy)
   s_oxz = ( s_oxz/d1 - m_ox*m_oz)
   s_oyy = ( s_oyy/d1 - m_oy**2 )
   s_oyz = ( s_oyz/d1 - m_oy*m_oz)
   s_ozz = ( s_ozz/d1 - m_oz**2 )

   m_w = m_w + 1d0-mes_D%r(:,2)

   
   
   
   !if(mpi_rnk==0) then
   !   print*, 'Writing vel_meanstdv.dat ...'
   !   open(99, status='unknown', file='vel_meanstdv.dat')
   !   do n = 1, i_N
   !      write(99,'(19e16.8)') mes_D%r(n,1),  &
   !         m_u(n), m_v(n), m_w(n),m_ox(n), m_oy(n), m_oz(n),&
   !         s_uu(n),s_uv(n),s_uw(n),s_vv(n),s_vw(n),s_ww(n),&
   !	   s_oxx(n),s_oxy(n),s_oxz(n),s_oyy(n),s_oyz(n),s_ozz(n)
   !   end do
   !   close(99)
   !endif
   

   

    if(mpi_rnk==0) then
    	print*, ' saving statistics'
    
    	e=nf90_create('State.cdf.dat', nf90_clobber, f)
    	e=nf90_put_att(f, nf90_global, 'Re', d_Re)
    	e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
    	e=nf90_def_dim(f, 'r', i_N, rd)
    	e=nf90_def_dim(f, 'K', 2*i_K-1, Kd)
    	e=nf90_def_dim(f, 'M', i_M, Md)
    
    	e=nf90_def_var(f, 'r', nf90_double, (/rd/), rid)
    	e=nf90_def_var(f, 'U', nf90_double, (/rd/), muid)
    	e=nf90_def_var(f, 'V', nf90_double, (/rd/), mvid)
        e=nf90_def_var(f, 'W', nf90_double, (/rd/), mwid)
        e=nf90_def_var(f, 'Omegax', nf90_double, (/rd/), moxid)
    	e=nf90_def_var(f, 'Omegay', nf90_double, (/rd/), moyid)
        e=nf90_def_var(f, 'Omegaz', nf90_double, (/rd/), mozid)

        e=nf90_def_var(f, 'uu', nf90_double, (/rd/), suuid)
	e=nf90_def_var(f, 'uv', nf90_double, (/rd/), suvid)
	e=nf90_def_var(f, 'uw', nf90_double, (/rd/), suwid)
	e=nf90_def_var(f, 'vv', nf90_double, (/rd/), svvid)
	e=nf90_def_var(f, 'vw', nf90_double, (/rd/), svwid)
	e=nf90_def_var(f, 'ww', nf90_double, (/rd/), swwid)

	e=nf90_def_var(f, 'uuu', nf90_double, (/rd/), suuuid)
	e=nf90_def_var(f, 'vvv', nf90_double, (/rd/), svvvid)
	e=nf90_def_var(f, 'www', nf90_double, (/rd/), swwwid)
	e=nf90_def_var(f, 'u4', nf90_double, (/rd/), su4id)
	e=nf90_def_var(f, 'v4', nf90_double, (/rd/), sv4id)
	e=nf90_def_var(f, 'w4', nf90_double, (/rd/), sw4id)

	e=nf90_def_var(f, 'oxx', nf90_double, (/rd/), soxxid)
	e=nf90_def_var(f, 'oxy', nf90_double, (/rd/), soxyid)
	e=nf90_def_var(f, 'oxz', nf90_double, (/rd/), soxzid)
	e=nf90_def_var(f, 'oyy', nf90_double, (/rd/), soyyid)
	e=nf90_def_var(f, 'oyz', nf90_double, (/rd/), soyzid)
	e=nf90_def_var(f, 'ozz', nf90_double, (/rd/), sozzid)


    	e=nf90_enddef(f)


    	e=nf90_put_var(f, rid, mes_D%r(1:i_N,1))
        e=nf90_put_var(f, muid, m_u)
	e=nf90_put_var(f, mvid, m_v)
	e=nf90_put_var(f, mwid, m_w)
	e=nf90_put_var(f, moxid, m_ox)
	e=nf90_put_var(f, moyid, m_oy)
	e=nf90_put_var(f, mozid, m_oz)
       
	e=nf90_put_var(f, suuid, s_uu)
	e=nf90_put_var(f, suvid, s_uv)
	e=nf90_put_var(f, suwid, s_uw)
	e=nf90_put_var(f, svvid, s_vv)
	e=nf90_put_var(f, svwid, s_vw)
	e=nf90_put_var(f, swwid, s_ww)

	e=nf90_put_var(f, suuuid, s_uuu)
	e=nf90_put_var(f, svvvid, s_vvv)
	e=nf90_put_var(f, swwwid, s_www)
	e=nf90_put_var(f, su4id, s_u4)
	e=nf90_put_var(f, sv4id, s_v4)
	e=nf90_put_var(f, sw4id, s_w4)

	e=nf90_put_var(f, soxxid, s_oxx)
	e=nf90_put_var(f, soxyid, s_oxy)
	e=nf90_put_var(f, soxzid, s_oxz)
	e=nf90_put_var(f, soyyid, s_oyy)
	e=nf90_put_var(f, soyzid, s_oyz)
	e=nf90_put_var(f, sozzid, s_ozz)

        e=nf90_close(f)
    end if

#ifdef _MPI
   call mpi_barrier(mpi_comm_world, mpi_er)
   call mpi_finalize(mpi_er)
#endif

!*************************************************************************
 END PROGRAM LOADSTATE
!*************************************************************************

