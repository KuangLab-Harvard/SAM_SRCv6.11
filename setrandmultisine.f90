subroutine setrandmultisine
  !compute the random phases for the random phase multisine at 1st step
  !then compute the step perturbation to field every nperturbstep steps
  !output at the same timestep, but show effect starting next output
  use vars
  use params
  use microphysics, only: micro_field, index_water_vapor

  implicit none
  integer i,j,k,f
  real pii,factor
  real delt_perturbt_current, delt_perturbq_current

  !-----note-----------------------------------------------------------------
  !unnormalized ttend_random and qtend_random are in K or kg/kg per 900s
  !not per day,
  !then transform into step perturbation (every nperturbstep steps)
  !--------------------------------------------------------------------------
  
  pii = atan2(0.d0,-1.d0)
    !update the forcing every nperturbstep time steps
    if(nstep.gt.nstartperturb.and.mod(nstep-nstartperturb,nperturbstep).eq.0) then
       !index for the nperturbstep*dt chunks
       !we update the forcing at the end of timesteps so the forcing for the very first chunk of the first period is zero.
       !But since the first period is spin up, this is okay. The subsequent periods don't have this problem.
       i=(nstep-nstartperturb)/nperturbstep
       i=mod(i,lenP) !reuse the same set of random numbers after 200 days
       ttend_random=0.
       qtend_random=0.
       !include all frequency forcings
       do k=1,lenP/freqfactor
          if(dorandmultisineoddonly) then
             f = 2*k-1
          else
             f = k
          end if
          !this is to reduce ampltude of lower frequencies by this factor:
          !new formula: 3.5-1.25*max(log10(1.5*freq),0)
          if(dowhitenoiseforcing) then
             factor=3.5D0
          else
             factor=max(1.D0,3.5-1.25*max(log10(1.5D0*real(f)/200.D0),0.D0))
          end if
          do j=1,nT
             ttend_random(j)=ttend_random(j)+cos((f*real(i)/real(lenP)+randcoef(j,k))*2.*pii)/factor
          end do
          do j=1,nQ
             qtend_random(j)=qtend_random(j)+cos((f*real(i)/real(lenP)+randcoef(j+nT,k))*2.*pii)/factor
          end do
       end do


       ! linearly increase perturbation amplitude (Qiyu 2024)
       delt_perturbt_current = delt_perturbt
       delt_perturbq_current = delt_perturbq
       if(increase_delt) then
          if(nstep.gt.nstep_increase_start.and.nstep.lt.nstep_increase_end) then
             delt_perturbt_current = delt_perturbt_current + (delt_perturbt_end-delt_perturbt) * &
                     (nstep-nstep_increase_start)/(nstep_increase_end-nstep_increase_start)
             delt_perturbq_current = delt_perturbq_current + (delt_perturbq_end-delt_perturbq) * &
                     (nstep-nstep_increase_start)/(nstep_increase_end-nstep_increase_start)
          else if (nstep.gt.nstep_increase_end) then
             delt_perturbt_current = delt_perturbt_end
             delt_perturbq_current = delt_perturbq_end
          end if
       end if


       !transform into step perturbation (every nperturbstep steps)
       !"unitary" change is 0.01K or qg0/960.kg/kg (in 900s)
       !sqrt(nperturbstep*dt/900.): normalized by 900 seconds
       !delt_perturbt: determine sign (+1 or -1), seems useless
       !1/sqrt(...): keep variance same as a single sine signal
       ttend_random=ttend_random*0.01*sqrt(nperturbstep*dt/900.)*delt_perturbt_current/sqrt(real(lenP/freqfactor))
       do j=1,nQ
          qtend_random(j)=qtend_random(j)*qg0(j)/960.*sqrt(nperturbstep*dt/900.)*delt_perturbq_current/sqrt(real(lenP/freqfactor))
       end do
    end if

end subroutine setrandmultisine
