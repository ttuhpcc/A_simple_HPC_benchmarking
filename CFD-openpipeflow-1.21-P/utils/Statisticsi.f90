#include "../parallel.h"
!*************************************************************************
 PROGRAM LOADSTATE
!*************************************************************************
   use io
   implicit none
   double precision :: M_r(i_N), M_t(i_N),M_z(i_N),M_zt(i_N),M_or(i_N), M_ot(i_N),M_oz(i_N)
   double precision :: s_rr(i_N),s_rt(i_N),s_rz(i_N),s_tt(i_N),s_tz(i_N),s_zz(i_N)
   double precision :: s_rrr(i_N),s_ttt(i_N),s_zzz(i_N),s_r4(i_N),s_t4(i_N),s_z4(i_N)
   double precision :: s_rtt(i_N),s_rzz(i_N),s_rrz(i_N)
   double precision :: s_orr(i_N),s_ort(i_N),s_orz(i_N),s_ott(i_N),s_otz(i_N),s_ozz(i_N)
   double precision :: s_rot(i_N),s_tor(i_N)
   double precision :: M_urr(i_N), M_utr(i_N),M_uzr(i_N),M_urt(i_N),M_utt(i_N),M_uzt(i_N),M_urz(i_N),M_utz(i_N),M_uzz(i_N)
   double precision :: S_urr(i_N), S_utr(i_N),S_uzr(i_N),S_urt(i_N),S_utt(i_N),S_uzt(i_N),S_urz(i_N),S_utz(i_N),S_uzz(i_N)
   double precision :: S_urruzr(i_N), S_urtuzt(i_N),S_urzuzz(i_N)
   double precision :: M_p(i_N), S_pp(i_N),S_pr(i_N),S_pt(i_N),S_pz(i_N),S_purr(i_N),S_putt(i_N),S_puzz(i_N),S_purzuzr(i_N) 
   double precision :: d1, d(i_N),a(i_N)
   character(8) :: cnum
   integer :: e, f, kd,md,rd,rid,intrdrid,iKL1d,iKL2d,iKL3d,drid,dr0id,dr1id
   integer :: mrid,mtid,mzid,mztid,morid,motid,mozid,srrid,srtid,srzid,sttid,stzid,szzid
   integer :: srrrid,stttid,szzzid,srttid,srzzid,srrzid,sr4id,st4id,sz4id,sorrid,sortid,sorzid,sottid,sotzid,sozzid
   integer :: srotid,storid
   integer :: murrid,mutrid,muzrid,murtid,muttid,muztid,murzid,mutzid,muzzid
   integer :: surrid,sutrid,suzrid,surtid,suttid,suztid,surzid,sutzid,suzzid, surruzrid, surtuztid, surzuzzid
   integer :: mpid,sppid,sprid,sptid,spzid,spurrid,sputtid,spuzzid,spurzuzrid
   integer :: n,n_, ff,fl,fi,fs,fn
   type (coll) :: c1,c2,c3,c4,c5,c6
   type (phys) :: pr,pt,pz,p1,p2,p3
   type (lumesh) ::  LNp(0:i_pH1)
   double precision :: BCR(0:i_pH1), BCI(0:i_pH1)
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
      read(*,*)  ff, fl, fn
   else
      open(99,status='old',file='MEANFLUCT.IN')
      read(99,*) ff, fl,fn
      close(99)
   end if

   fs=0 
   do fi = ff, fl,fn
      	fs=fs+1;
   	m_r = 0d0; m_t = 0d0; m_z = 0d0; m_zt = 0d0; m_or = 0d0;m_ot = 0d0; M_oz = 0d0
   	s_rr= 0d0; s_rt= 0d0; s_rz= 0d0; s_tt= 0d0; s_tz= 0d0; s_zz= 0d0;
   	s_rrr= 0d0; s_ttt= 0d0; s_zzz= 0d0; s_r4= 0d0; s_t4= 0d0; s_z4= 0d0;
   	s_rtt= 0d0; s_rzz= 0d0; s_rrz= 0d0;
   	s_orr= 0d0; s_ort= 0d0; s_orz= 0d0; s_ott= 0d0; s_otz= 0d0; s_ozz= 0d0;
   	s_rot= 0d0; s_tor= 0d0;
      
   	M_urr=0d0;  M_utr=0d0;  M_uzr=0d0;  M_urt=0d0;  M_utt=0d0;  M_uzt=0d0;
   	M_urz=0d0;  M_utz=0d0;  M_uzz=0d0;
   	S_urr=0d0;  S_utr=0d0;  S_uzr=0d0;  S_urt=0d0;  S_utt=0d0;  S_uzt=0d0;
   	S_urz=0d0;  S_utz=0d0;  S_uzz=0d0; S_urruzr=0d0;S_urtuzt=0d0;S_urzuzz=0d0
   	M_p=0; S_pp=0; S_pr=0; s_pt=0; s_pz=0; s_purr=0; s_putt=0; s_puzz=0;
   	s_purzuzr=0

   	write(cnum,'(I8.8)') fi
      	io_statefile = 'state'//cnum//'.cdf.dat'
      	call io_load_state()
      	call vel_transform()
      	! statistic 1:
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

	 	s_r4(n_) = s_r4(n_) + sum(vel_r%Re(:,:,n)**4)
         	s_t4(n_) = s_t4(n_) + sum(vel_t%Re(:,:,n)**4)
	 	s_z4(n_) = s_z4(n_) + sum(vel_z%Re(:,:,n)**4)

	 	s_orr(n_) = s_orr(n_) + sum(vel_curlr%Re(:,:,n)**2)
	 	s_ott(n_) = s_ott(n_) + sum(vel_curlt%Re(:,:,n)**2)
	 	s_ozz(n_) = s_ozz(n_) + sum(vel_curlz%Re(:,:,n)**2)
	 	s_ort(n_) = s_ort(n_) + sum(vel_curlr%Re(:,:,n)*vel_curlt%Re(:,:,n))
	 	s_orz(n_) = s_orz(n_) + sum(vel_curlr%Re(:,:,n)*vel_curlz%Re(:,:,n))
	 	s_otz(n_) = s_otz(n_) + sum(vel_curlt%Re(:,:,n)*vel_curlz%Re(:,:,n))
 
         	s_rot(n_) = s_rot(n_) + sum(vel_r%Re(:,:,n)*vel_curlt%Re(:,:,n))
	 	s_tor(n_) = s_tor(n_) + sum(vel_t%Re(:,:,n)*vel_curlr%Re(:,:,n))
	
      end do

 !! Statistics 3: TKE

      call var_coll_copy(vel_ur, c1)
      call var_coll_copy(vel_ut, c2)
      call var_coll_copy(vel_uz, c3)
      call var_coll_curl(vel_ur,vel_ut,vel_uz, c4,c5,c6)
      call var_coll_curl(c4,c5,c6, c4,c5,c6)  

      tim_dt = 0.0001d0
      call vel_matrices()
      call vel_transform()
      call vel_nonlinear()
      call vel_predictor()
      tim_it = 1

      do while(tim_it/=0)
         call vel_transform()
         call vel_nonlinear()
         call var_null(2)
         call vel_corrector()
         call tim_check_cgce()
      end do
      c1%Re = vel_Nr%Re - (vel_ur%Re-c1%Re)/tim_dt - c4%Re/d_Re
      c1%Im = vel_Nr%Im - (vel_ur%Im-c1%Im)/tim_dt - c4%Im/d_Re
      c2%Re = vel_Nt%Re - (vel_ut%Re-c2%Re)/tim_dt - c5%Re/d_Re
      c2%Im = vel_Nt%Im - (vel_ut%Im-c2%Im)/tim_dt - c5%Im/d_Re
      c3%Re = vel_Nz%Re - (vel_uz%Re-c3%Re)/tim_dt - c6%Re/d_Re
      c3%Im = vel_Nz%Im - (vel_uz%Im-c3%Im)/tim_dt - c6%Im/d_Re
         
      _loop_km_begin
            if(m/=0) then
               c1%Re(:,nh) =  c2%Im(:,nh)*mes_D%r(:,1)/dble(i_Mp*m)
               c1%Im(:,nh) = -c2%Re(:,nh)*mes_D%r(:,1)/dble(i_Mp*m)
            else if(k/=0) then
               c1%Re(:,nh) =  c3%Im(:,nh)/(d_alpha*k)
               c1%Im(:,nh) = -c3%Re(:,nh)/(d_alpha*k)
            else
               a(1) = 0d0
               do n = 2, i_N
                  a(n) = a(n-1) + c1%Re(n,0)*mes_D%intrdr(n)/mes_D%r(n,1)
               end do
               c1%Re(:,0) = a(:) - a(i_N)
               c1%Im(:,0) = 0d0
            end if
      _loop_km_end
      call tra_coll2phys(c1, pr)
      pr%Re = pr%Re - 0.5d0*  &
            (vel_r%Re*vel_r%Re + vel_t%Re*vel_t%Re + vel_z%Re*vel_z%Re)         
         
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1        
	 M_p(n_) = M_p(n_) + sum(pr%Re(:,:,n))
         s_pp(n_) = s_pp(n_) + sum(pr%Re(:,:,n)**2)        	 
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

   call mpi_allreduce(s_r4, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_r4 = d
   call mpi_allreduce(s_t4, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_t4 = d
   call mpi_allreduce(s_z4, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_z4 = d


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

   call mpi_allreduce(s_urruzr, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_urruzr = d
   call mpi_allreduce(s_urtuzt, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_urtuzt = d
   call mpi_allreduce(s_urzuzz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   s_urzuzz = d


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

   !d1 = dble((fl-ff+1) * i_Th*i_Z)
   d1 = dble(1 * i_Th*i_Z)
   m_r = m_r / d1
   m_t = m_t / d1
   m_z = m_z / d1
   m_or = m_or / d1
   m_ot = m_ot / d1
   m_oz = m_oz / d1
   

   s_rr = ( s_rr/d1 - m_r**2 )
   s_rt = ( s_rt/d1 - m_r*m_t)
   s_rz = ( s_rz/d1 - m_r*m_z)
   s_tt = ( s_tt/d1 - m_t**2 )
   s_tz = ( s_tz/d1 - m_t*m_z)
   s_zz = ( s_zz/d1 - m_z**2 )

   s_rtt = ( s_rtt/d1 - m_r*m_t**2)
   s_rzz = ( s_rzz/d1 - m_r*m_z**2-m_r*s_zz-2*m_z*s_rz)
   s_rrz = ( s_rrz/d1 - m_r**2*m_z-s_rr*m_z)

   s_rrr = ( s_rrr/d1 - m_r**3)
   s_ttt = ( s_ttt/d1 - m_t**3)
   s_zzz = ( s_zzz/d1 - m_z**3-3*m_z*s_zz)

   s_r4 = ( s_r4/d1 - m_r**4)
   s_t4 = ( s_t4/d1 - m_t**4)
   s_z4 = ( s_z4/d1 - m_z**4-6*m_z**2*s_zz-4*m_z*s_zzz)

  
   s_orr = ( s_orr/d1 - m_or**2 )
   s_ort = ( s_ort/d1 - m_or*m_ot)
   s_orz = ( s_orz/d1 - m_or*m_oz)
   s_ott = ( s_ott/d1 - m_ot**2 )
   s_otz = ( s_otz/d1 - m_ot*m_oz)
   s_ozz = ( s_ozz/d1 - m_oz**2 )

   s_rot = ( s_rot/d1 - m_ot*m_r)
   s_tor = ( s_tor/d1 - m_or*m_t)

   m_urr = m_urr / d1
   m_urt = m_urt / d1
   m_urz = m_urz / d1
   m_utr = m_utr / d1
   m_utt = m_utt / d1
   m_utz = m_utz / d1
   m_uzr = m_uzr / d1
   m_uzt = m_uzt / d1
   m_uzz = m_uzz / d1

   s_urr = (s_urr / d1 - m_urr**2 )
   s_urt = (s_urt / d1 - m_urt**2 )
   s_urz = (s_urz / d1 - m_urz**2 )
   s_utr = (s_utr / d1 - m_utr**2 )
   s_utt = (s_utt / d1 - m_utt**2 )
   s_utz = (s_utz / d1 - m_utz**2 )
   s_uzr = (s_uzr / d1 - m_uzr**2 )
   s_uzt = (s_uzt / d1 - m_uzt**2 )
   s_uzz = (s_uzz / d1 - m_uzz**2 )

   s_urruzr = (s_urruzr / d1- m_urr* m_uzr)
   s_urtuzt = (s_urtuzt / d1- m_urt* m_uzt)
   s_urzuzz = (s_urzuzz / d1- m_urz* m_uzz)

   m_p = m_p / d1
   s_pp = ( s_pp/d1 - m_p**2 )
   s_pr = ( s_pr/d1 - m_p*m_r )
   s_pt = ( s_pt/d1 - m_p*m_t )
   s_pz = ( s_pz/d1 - m_p*m_z )
 
    s_purr = ( s_purr/d1 - m_p*m_urr )
    s_putt = ( s_putt/d1 - m_p*m_utt )
    s_puzz = ( s_puzz/d1 - m_p*m_uzz )
    s_purzuzr = ( s_purzuzr/d1 - m_p*(m_urz+m_uzr))

    m_zt = m_z + 1d0-mes_D%r(:,2)



    

    if(mpi_rnk==0) then
    	print*, ' saving statistics'
        e=nf90_create('Statist'//cnum//'.cdf.dat', NF90_64BIT_OFFSET, f) 
    	!e=nf90_create('State.cdf.dat', nf90_clobber, f)
    	e=nf90_put_att(f, nf90_global, 'Re', d_Re)
    	e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
    	e=nf90_def_dim(f, 'r', i_N, rd)
    	e=nf90_def_dim(f, 'K', 2*i_K-1, Kd)
    	e=nf90_def_dim(f, 'M', i_M, Md)
        e=nf90_def_dim(f, 'iKL1', 2*i_KL+1, iKL1d)
        e=nf90_def_dim(f, 'iKL2', i_N+i_KL, iKL2d)
        e=nf90_def_dim(f, 'iKL3', i_KL+1, iKL3d)
    
    	e=nf90_def_var(f, 'r', nf90_double, (/rd/), rid)
        e=nf90_def_var(f, 'intrdr', nf90_double, (/rd/), intrdrid)
        e=nf90_def_var(f, 'dr', nf90_double, (/iKL1d,iKL2d/), drid)
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
	e=nf90_def_var(f, 'ur4', nf90_double, (/rd/), sr4id)
	e=nf90_def_var(f, 'ut4', nf90_double, (/rd/), st4id)
	e=nf90_def_var(f, 'uz4', nf90_double, (/rd/), sz4id)

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

	e=nf90_def_var(f, 'dsurruzr', nf90_double, (/rd/), surruzrid)   
        e=nf90_def_var(f, 'dsurtuzt', nf90_double, (/rd/), surtuztid)  
	e=nf90_def_var(f, 'dsurzuzz', nf90_double, (/rd/), surzuzzid)  

    	e=nf90_def_var(f, 'P', nf90_double, (/rd/), mpid)   
        e=nf90_def_var(f, 'pp', nf90_double, (/rd/), sppid)
	e=nf90_def_var(f, 'pr', nf90_double, (/rd/), sprid)
	e=nf90_def_var(f, 'pt', nf90_double, (/rd/), sptid)
	e=nf90_def_var(f, 'pz', nf90_double, (/rd/), spzid)
        e=nf90_def_var(f, 'purr', nf90_double, (/rd/), spurrid)
        e=nf90_def_var(f, 'putt', nf90_double, (/rd/), sputtid)
        e=nf90_def_var(f, 'puzz', nf90_double, (/rd/), spuzzid)
	e=nf90_def_var(f, 'purzuzr', nf90_double, (/rd/), spurzuzrid)
	


    	e=nf90_enddef(f)


    	e=nf90_put_var(f, rid, mes_D%r(1:i_N,1))
        e=nf90_put_var(f, intrdrid, mes_D%intrdr)
	!e=nf90_put_var(f, drid, mes_D%dr, start=(/1,1/))
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
	e=nf90_put_var(f, sr4id, s_r4)
	e=nf90_put_var(f, st4id, s_t4)
	e=nf90_put_var(f, sz4id, s_z4)

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

	e=nf90_put_var(f, surruzrid, s_urruzr) 
	e=nf90_put_var(f, surtuztid, s_urtuzt) 
	e=nf90_put_var(f, surzuzzid, s_urzuzz) 

	e=nf90_put_var(f, mpid, m_p)    
	e=nf90_put_var(f, sppid, s_pp)
	e=nf90_put_var(f, sprid, s_pr)
	e=nf90_put_var(f, sptid, s_pt)
	e=nf90_put_var(f, spzid, s_pz)
        e=nf90_put_var(f, spurrid, s_purr)
	e=nf90_put_var(f, sputtid, s_putt)
	e=nf90_put_var(f, spuzzid, s_puzz)
	e=nf90_put_var(f, spurzuzrid, s_purzuzr)
       

        e=nf90_close(f)
    end if
    end do

#ifdef _MPI
   call mpi_barrier(mpi_comm_world, mpi_er)
   call mpi_finalize(mpi_er)
#endif

!*************************************************************************
 END PROGRAM LOADSTATE
!*************************************************************************

