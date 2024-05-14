
subroutine setperturb

!  Random noise

use vars
use params
use microphysics, only: micro_field, index_water_vapor
use sgs, only: setperturb_sgs

implicit none

! Modified so that there is an option to produced identical perturbations regardless of 
!  the parallel decomposition. All random numbers are computed on masterproc
!  and scattered to the other processes, by replacing scalar rrr with
!  rrr_loc and rrr_glob
! Only an option because we don't want a global size array for very large problems.
! Default will be false, so user will have to go in here and reset it true and recompile.
! Don Dazlich, dazlich@atmos.colostate.edu 2 Feb 2017

integer i,j,k,ptype,it,jt, n, sender
real rrr,ranf_
real xxx,yyy,zzz
logical :: do_decomp_independent_perturb = .false.
real rrr_loc(nx,ny,nzm)
real, allocatable, dimension(:,:,:) :: rrr_gl

call ranset_(3*rank+icopy*10000)

ptype = perturb_type

call setperturb_sgs(ptype)  ! set sgs fields

if(do_decomp_independent_perturb) then
   if(masterproc) then
      allocate(rrr_gl(nx_gl,ny_gl,nzm))
      do k=1,nzm
       do j=1,ny_gl
        do i=1,nx_gl
           rrr_gl(i,j,k) = ranf_()
        end do
       end do
      end do
   endif
   
   if(masterproc) then
     do n = 1, nsubdomains-1
! get the array section for rank n and send it
        call task_rank_to_index(n,it,jt)
        do k=1,nzm
         do j=1,ny
          do i=1,nx
             rrr_loc(i,j,k) = rrr_gl(i+it,j+jt,k)
          end do
         end do
        end do
 
        call task_bsend_float(n,rrr_loc,nx*ny*nzm,n)
     enddo

! get the masterproc array section
     do k=1,nzm
      do j=1,ny
       do i=1,nx
          rrr_loc(i,j,k) = rrr_gl(i,j,k)
       end do
      end do
     end do
    else

! receive array section from masterproc
       call task_breceive_float(rrr_loc,nx*ny*nzm,sender,rank)
    endif
     
   call task_barrier()
   if(masterproc) deallocate(rrr_gl)
else
     do k=1,nzm
      do j=1,ny
       do i=1,nx
          rrr_loc(i,j,k) = ranf_()
       end do
      end do
     end do
endif

select case (ptype)

  case(-1)  ! no perturbation is set

  case(0)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*rrr_loc(i,j,k)
         if(k.le.5) then
            t(i,j,k)=t(i,j,k)+0.02*rrr*(6-k)
         endif
       end do
      end do
     end do

  case(1)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*rrr_loc(i,j,k)
         if(q0(k).gt.6.e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
         endif
       end do
      end do
     end do

  case(2) ! warm bubble

     if(masterproc) then
       print*, 'initialize with warm bubble:'
       print*, 'bubble_x0=',bubble_x0
       print*, 'bubble_y0=',bubble_y0
       print*, 'bubble_z0=',bubble_z0
       print*, 'bubble_radius_hor=',bubble_radius_hor
       print*, 'bubble_radius_ver=',bubble_radius_ver
       print*, 'bubble_dtemp=',bubble_dtemp
       print*, 'bubble_dq=',bubble_dq
     end if

     call task_rank_to_index(rank,it,jt)
     do k=1,nzm
       zzz = z(k)
       do j=1,ny
         yyy = dy*(j+jt)
         do i=1,nx
          xxx = dx*(i+it)
           if((xxx-bubble_x0)**2+YES3D*(yyy-bubble_y0)**2.lt.bubble_radius_hor**2 &
            .and.(zzz-bubble_z0)**2.lt.bubble_radius_ver**2) then
              rrr = cos(pi/2.*(xxx-bubble_x0)/bubble_radius_hor)**2 &
               *cos(pi/2.*(yyy-bubble_y0)/bubble_radius_hor)**2 &
               *cos(pi/2.*(zzz-bubble_z0)/bubble_radius_ver)**2
              t(i,j,k) = t(i,j,k) + bubble_dtemp*rrr
              micro_field(i,j,k,index_water_vapor) = &
                  micro_field(i,j,k,index_water_vapor) + bubble_dq*rrr
           end if
         end do
       end do
     end do

  case(3)   ! gcss wg1 smoke-cloud case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*rrr_loc(i,j,k)
         if(q0(k).gt.0.5e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
         endif
       end do
      end do
     end do

  case(4)  ! gcss wg1 arm case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*rrr_loc(i,j,k)
         if(z(k).le.200.) then
            t(i,j,k)=t(i,j,k)+0.1*rrr*(1.-z(k)/200.)
         endif
       end do
      end do
     end do

  case(5)  ! gcss wg1 BOMEX case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*rrr_loc(i,j,k)
         if(z(k).le.1600.) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)+0.025e-3*rrr
         endif
       end do
      end do
     end do

  case(6)  ! GCSS Lagragngian ASTEX


     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*rrr_loc(i,j,k)
         if(q0(k).gt.6.e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)+2.5e-5*rrr
         endif
       end do
      end do
     end do


  case default

       if(masterproc) print*,'perturb_type is not defined in setperturb(). Exitting...'
       call task_abort()

end select


end

