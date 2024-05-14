subroutine perturbtq
!Editted by Qiyu Song based on earlier version from Zhiming Kuang
!Add perturbation to one layer of T,q field
  use vars
  use params
  use microphysics, only: micro_field, index_water_vapor

  implicit none
  
  !iperturb1 varies from 1 to nT+nQ
  if(iperturb1.le.nT) then
    !perturb temperature layer
    iperturbt=iperturb1
    iperturbq=0
  else
    !perturb moisture layer
    iperturbt=0
    iperturbq=iperturb1-nT
  end if

  if(nstep.gt.nstartperturb.and.douniformintime) then
    !spread the random forcing uniformly over the nstat step
    if(iperturbt.gt.0) then
      t(:,:,iperturbt)=t(:,:,iperturbt)+ttend_random(iperturbt)/nperturbstep
    end if
    if(iperturbq.gt.0) then
      micro_field(:,:,iperturbq,index_water_vapor)=micro_field(:,:,iperturbq,index_water_vapor)+qtend_random(iperturbq)/nperturbstep
    end if

  elseif(nstep.gt.nstartperturb.and.mod(nstep-nstartperturb,nperturbstep).eq.0) then
    !only do the forcing every nstat steps
    if(iperturbt.gt.0) then
      t(:,:,iperturbt)=t(:,:,iperturbt)+ttend_random(iperturbt)
    end if
    if(iperturbq.gt.0) then
      micro_field(:,:,iperturbq,index_water_vapor)=micro_field(:,:,iperturbq,index_water_vapor)+qtend_random(iperturbq)
    end if
  end if

end subroutine perturbtq
