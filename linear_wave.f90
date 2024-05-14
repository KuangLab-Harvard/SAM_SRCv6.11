subroutine linear_wave

! Original version from Zhiming Kuang, Nov. 27, 2006
! Editted by Qiyu Song (2024)
! SAM specific routine parameterized large-scale linear wave dynamics

use vars
use params
use module_linearwave
implicit none
integer k,i,j,ilag

! local variables
real dwwavedt(nzm) ! tendency of large-scale wave vertical velocity
real wn            ! horizontal wavenumber
real coef, coef1, buffer(nzm,2), buffer1(nzm,2), buffer2(nzm*2,nzm*2), tmp(nzm)

dwwavedt = 0.

! 0 < nstep <= nstartlinearwave:
!    run to equilibrium
! nstartlinearwave < nstep < nstartlinearwave + nsteplinearwavebg:
!    compute background profiles
! nstep = nstartlinearwave + nsteplinearwavebg + n*nsteplinearwave:
!    wave coupling
if(nstep.gt.nstartlinearwave.and.nstep.le.nstartlinearwave+nsteplinearwavebg) then
   t_wavebg_dble=t_wavebg_dble+dble(tabs0)
   q_wavebg_dble=q_wavebg_dble+dble(qv0)
   rad_wavebg_dble=rad_wavebg_dble+dble(radqrlw+radqrsw)/dble(nx*ny)
endif

if (nstep.eq.nstartlinearwave+nsteplinearwavebg) then
   t_wavebg=t_wavebg_dble/dble(nsteplinearwavebg)
   q_wavebg=q_wavebg_dble/dble(nsteplinearwavebg)
   rad_wavebg=rad_wavebg_dble/dble(nsteplinearwavebg)
   !get the averages of all ensembles (i.e. subdomains)
   coef1 = 1./dble(nsubdomains)
   do k=1,nzm
      buffer(k,1) = t_wavebg(k)
      buffer(k,2) = q_wavebg(k)
   end do
   call task_sum_real8(buffer,buffer1,nzm*2)
   do k=1,nzm
      t_wavebg(k)=buffer1(k,1)*coef1
      q_wavebg(k)=buffer1(k,2)*coef1
   end do
   
   ! also initialize t_wave and q_wave here
   t_wave_local=tabs0
   q_wave_local=qv0
   !average T, q over nsteplinearwave
   t_wave=t_wave_local
   q_wave=q_wave_local
   coef1 = 1./dble(nsubdomains)
   do k=1,nzm
      buffer(k,1) = t_wave(k)
      buffer(k,2) = q_wave(k)
   end do
   call task_sum_real8(buffer,buffer1,nzm*2)
   do k=1,nzm
      t_wave(k)=buffer1(k,1)*coef1
      q_wave(k)=buffer1(k,2)*coef1
   end do


end if

if(nstep.gt.nstartlinearwave+nsteplinearwavebg) then
   !update wwave every nsteplinearwave steps
   if(mod(nstep-nstartlinearwave-nsteplinearwavebg,nsteplinearwave).eq.0) then
      t_wave_local=tabs0
      q_wave_local=qv0
      !average T, q over nsteplinearwave
      t_wave=t_wave_local!/float(nsteplinearwave)
      q_wave=q_wave_local!/float(nsteplinearwave)
      coef1 = 1./dble(nsubdomains)
      do k=1,nzm
         buffer(k,1) = t_wave(k)
         buffer(k,2) = q_wave(k)
      end do
      call task_sum_real8(buffer,buffer1,nzm*2)
      do k=1,nzm
         t_wave(k)=buffer1(k,1)*coef1
         q_wave(k)=buffer1(k,2)*coef1
      end do

      !before updating the anomalies, compute the lagged correlation for the noise
      !This part is not coded yet

      !compute virtual temperature
      tv_wavebg=t_wavebg*(1+0.61*q_wavebg)
      tv_wave=t_wave*(1+0.61*q_wave)
      if(dointernalnoise) then
         tv_wave=t_wave_local*(1+0.61*q_wave_local)
      end if

      !horizontal wavenumber
      wn=wavenumber_factor*6.283e-6 !2pi/1000km
      !radiation upper boundary condition
      !exclude the top layer because of doupperbound
      !exclude the top two layers as doupperbound affect the top two layers
      dwwavedt=0.
      call calc_wtend(           wn,  w_wave(1:nzm-2), &
                 tv_wave(1:nzm-2),                           &
                 tv_wavebg(1:nzm-2),                       rho(1:nzm-2),    &
                 z(1:nzm-2),                          zi(1:nzm-1),        &
                 dwwavedt(1:nzm-2),                           nzm-2)
      w_wave=w_wave*(1.-float(nsteplinearwave)*dt/86400./wavedampingtime)+dwwavedt*dt*nsteplinearwave
      !Use only one ensemble member (rank=0) for coupling with the wave
      !so that the wave will be noisy. The wave forcing is then
      !broadcasted to every ensemble member so that all ensemble members
      !are driven the same way. And the convective tendendies and the states
      !are computed using the ensembel mean
      if(dointernalnoise) then
         if(masterproc)   tmp=w_wave
         call task_bcast_float4(0,tmp,nzm)
         w_wave=tmp
      end if

   end if
end if

end subroutine linear_wave
