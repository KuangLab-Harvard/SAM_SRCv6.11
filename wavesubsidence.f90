!Add vertical advection of q,w due to the large-scale wave and some largescale damping of T, q
subroutine wavesubsidence()

use vars
use params
use microphysics
use module_linearwave
implicit none

integer i,j,k
real ttend_center(nzm), qtend_center(nzm) !cell center tendencies
real rrr,ranf_,persistence,amp,fixTQdampingtime !for random forcing tendencies,
integer number_of_steps
double precision coef, coef1, buffer(nzm), buffer1(nzm)

if(nstep.gt.nstartlinearwave+nsteplinearwavebg.and.icycle.eq.ncycle) then
   ttend_center=0.
   qtend_center=0.

   !Obtain the tendency from advecting the background/reference T,q fields
   if(doadvectbg) then
      call linearwavesubsidence(           w_wave(1:nzm-1),      &
              t_wavebg(1:nzm-1),                       q_wavebg(1:nzm-1),    &
              z(1:nzm-1),                          ttend_center(1:nzm-1),        &
              qtend_center(1:nzm-1),                           nzm-1 )
   else
      !advect the full T,q fields
      call linearwavesubsidence(           w_wave(1:nzm-1),      &
              t_wave(1:nzm-1),                       q_wave(1:nzm-1),    &
              z(1:nzm-1),                          ttend_center(1:nzm-1),        &
              qtend_center(1:nzm-1),                           nzm-1 )
   endif

   number_of_steps=1
   do k=1,nzm
      water_fill(k)=0.
      do j=1,ny
         do i=1,nx
            t(i,j,k) = t(i,j,k) + number_of_steps*dt * (ttend_center(k)-(t_wave(k)-t_wavebg(k))/86400./wavetqdampingtime)
            !In SAM1MOM, micro_field(:,:,:,index_water_vapor) actually includes nonprecipitating water as well
            micro_field(i,j,k,index_water_vapor)=micro_field(i,j,k,index_water_vapor) + number_of_steps*dt *&
                  (qtend_center(k)-(q_wave(k)-q_wavebg(k))/86400./wavetqdampingtime)
            if (micro_field(i,j,k,index_water_vapor).lt.0) then
               water_fill(k)=water_fill(k)-micro_field(i,j,k,index_water_vapor)
               micro_field(i,j,k,index_water_vapor) = 0.
            end if
         end do
      end do
      water_fill(k)=water_fill(k)/dble(nx*ny)
   end do

   if(mod(nstep-nstartlinearwave-nsteplinearwavebg,nsteplinearwave).eq.0) then
      !compute T and q profiles after large-scale advection and before convection
      !because of the removal of negative qv above, this may incur a small error
      number_of_steps=nsteplinearwave
      coef=1./dble(nx*ny)
      do k=1,nzm
         bc_tabs0_local(k)=t_wave_local(k)+number_of_steps*dt*ttend_center(k)
         !there may be slight inaccuracies due to zeroing of negative water field as above
         bc_qv0_local(k)=q_wave_local(k)+water_fill(k)+number_of_steps*dt*qtend_center(k)
      end do
      coef1 = 1./dble(nsubdomains)
      buffer=bc_qv0_local
      call task_sum_real8(buffer,buffer1,nzm)
      bc_qv0=buffer1*coef1
      buffer=bc_tabs0_local
      call task_sum_real8(buffer,buffer1,nzm)
      bc_tabs0=buffer1*coef1
   
      !record the tendencies for output statistics
      ttend_wave = ttend_center
      qtend_wave = qtend_center
   end if
end if

end subroutine wavesubsidence

