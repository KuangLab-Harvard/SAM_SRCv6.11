subroutine setrandmultisinephases
  !compute the random phases for the random phase multisine at 1st step
  use vars
  use params
  
  implicit none
  real ranf_
  
  !integer,dimension(1)::seed=(/3/)
  integer j,k

  !call random_seed(put=seed+iensemble) !reset the random number generator
  !call ranset_(1000*(iensemble+1)) !run the random number generator for this many times (same across processors) to initialize
  
  lenP = floor(maxperturbperiod/(dt*nstat)+1.e-5)
  if(dorandmultisineoddonly) then
    freqfactor = 4
  else
    freqfactor = 2
  end if
  
  allocate(randcoef(nT+nQ, lenP/freqfactor))
  
  do j=1,nT+nQ
    do k=1,lenP/freqfactor
      randcoef(j,k)=ranf_()
    end do
  end do
  !print*,'randcoef(1,1):', randcoef(1,1), ', rank:', rank

end subroutine setrandmultisinephases
