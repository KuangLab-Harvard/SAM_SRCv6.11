! Non-blocking receives before blocking sends

subroutine pressure_mpiensemble
	
!       Original pressure solver based on horizontal slabs
!       (C) 1998, 2002 Marat Khairoutdinov
!       Works only when the number of slabs is equal to the number of processors.
!       Therefore, the number of processors shouldn't exceed the number of levels nzm
!       Also, used for a 2D version 
!       For more processors for the given number of levels and 3D, use pressure_big
!
!       MPI Ensemble run:
!       Modified by Song Qiyu (2022)
!       Renamed into a separate file/subroutine by Nathanael Wong (2022)

use vars
use params, only: dowallx, dowally, docolumn
implicit none
	
! Kuang Ensemble run: replace all 'nsubdomains' to 1, and all 'rank' to 0 (Song Qiyu, 2022)
! Kuang Ensemble run: replace 'nx_gl' to 'nx', and 'ny_gl' to 'ny' (Song Qiyu, 2022)
integer, parameter :: npressureslabs = 1
!solve in each processor separately
!each processor contains an ensemble member
integer, parameter :: nzslab = max(1,nzm / npressureslabs)
integer, parameter :: nx2=nx+2, ny2=ny+2*YES3D
integer, parameter :: n3i=3*nx/2+1,n3j=3*ny/2+1

real(8) f(nx2,ny2,nzslab) ! global rhs and array for FTP coefficeients
real ff(nx+1,ny+2*YES3D,nzm)	! local (subdomain's) version of f
real buff_slabs(nxp1,nyp2,nzslab,npressureslabs)
real buff_subs(nxp1,nyp2,nzslab,1) 
real bufp_slabs(0:nx,1-YES3D:ny,nzslab,npressureslabs)  
real bufp_subs(0:nx,1-YES3D:ny,nzslab,1)  
common/tmpstack/f,ff,buff_slabs,buff_subs
equivalence (buff_slabs,bufp_slabs)
equivalence (buff_subs,bufp_subs)

real(8) work(nx2,ny2),trigxi(n3i),trigxj(n3j) ! FFT stuff
integer ifaxj(100),ifaxi(100)

real(8) a(nzm),b,c(nzm),e,fff(nzm)	
real(8) xi,xj,xnx,xny,ddx2,ddy2,pii,factx,facty,eign
real(8) alfa(nzm-1),beta(nzm-1)

integer reqs_in(1)
integer i, j, k, id, jd, m, n, it, jt, ii, jj, tag, rf
integer nyp22, n_in, count
integer iii(0:nx),jjj(0:ny)
logical flag(1)
integer iwall,jwall

! check if the grid size allows the computation:

if(1.gt.nzm) then
  if(masterproc) print*,'pressure_orig: nzm < 1. STOP'
  call task_abort
endif

if(mod(nzm,npressureslabs).ne.0) then
  if(masterproc) print*,'pressure_orig: nzm/npressureslabs is not round number. STOP'
  call task_abort
endif

! for debug
!print*,'********************** rank = ***********************'
!print*,rank


!-----------------------------------------------------------------

if(docolumn) return

if(dowallx) then
  iwall=1
else
  iwall=0
end if
if(RUN2D) then  
  nyp22=1
  jwall=0
else
  nyp22=nyp2
  if(dowally) then
    jwall=2
  else
    jwall=0
  end if
endif
	
!-----------------------------------------------------------------
!  Compute the r.h.s. of the Poisson equation for pressure

call press_rhs()


!-----------------------------------------------------------------	 
!   Form the horizontal slabs of right-hand-sides of Poisson equation 
!   for the global domain. Request sending and receiving tasks.

! iNon-blocking receive first:

n_in = 0
do m = 0,1-1

  if(0.lt.npressureslabs.and.m.ne.1-1) then

    n_in = n_in + 1
    call task_receive_float(bufp_subs(0,1-YES3D,1,n_in), &
                           nzslab*nxp1*nyp1,reqs_in(n_in))
    flag(n_in) = .false.
 
  endif

  if(0.lt.npressureslabs.and.m.eq.1-1) then

    call task_rank_to_index(0,it,jt)	  
    n = 0*nzslab
    do k = 1,nzslab
     do j = 1,ny
       do i = 1,nx
         f(i+it,j+jt,k) = p(i,j,k+n)
       end do
     end do
    end do
  endif

end do ! m


! Blocking send now:


do m = 0,1-1

  if(m.lt.npressureslabs.and.m.ne.0) then

    n = m*nzslab + 1
    call task_bsend_float(m,p(0,1-YES3D,n),nzslab*nxp1*nyp1, 33)
  endif

end do ! m


! Fill slabs when receive buffers are full:

count = n_in
do while (count .gt. 0)
  do m = 1,n_in
   if(.not.flag(m)) then
	call task_test(reqs_in(m), flag(m), rf, tag)
        if(flag(m)) then 
	   count=count-1
           call task_rank_to_index(rf,it,jt)	  
           do k = 1,nzslab
            do j = 1,ny
             do i = 1,nx
               f(i+it,j+jt,k) = bufp_subs(i,j,k,m)
             end do
            end do
           end do
	endif   
   endif
  end do
end do


!-------------------------------------------------
! Perform Fourier transformation for a slab:

if(0.lt.npressureslabs) then

 call fftfax_crm(nx,ifaxi,trigxi)
 if(RUN3D) call fftfax_crm(ny,ifaxj,trigxj)

 do k=1,nzslab

  if(dowallx) then
   call cosft_crm(f(1,1,k),work,trigxi,ifaxi,1,nx2,nx,ny,-1)
  else
   call fft991_crm(f(1,1,k),work,trigxi,ifaxi,1,nx2,nx,ny,-1)
  end if

  if(RUN3D) then
   if(dowally) then
     call cosft_crm(f(1,1,k),work,trigxj,ifaxj,nx2,1,ny,nx+1,-1)
   else
     call fft991_crm(f(1,1,k),work,trigxj,ifaxj,nx2,1,ny,nx+1,-1)
   end if
  end if

 end do 

endif


! Synchronize all slabs:

call task_barrier()

!-------------------------------------------------
!   Send Fourier coeffiecients back to subdomains:

! Non-blocking receive first:

n_in = 0
do m = 0, 1-1
		
   call task_rank_to_index(m,it,jt)

   if(0.lt.npressureslabs.and.m.eq.0) then

     n = 0*nzslab
     do k = 1,nzslab
      do j = 1,nyp22-jwall
        do i = 1,nxp1-iwall
          ff(i,j,k+n) = f(i+it,j+jt,k) 
        end do
      end do
     end do 

   end if

   if(m.lt.npressureslabs-1.or.m.eq.npressureslabs-1 &
                            .and.0.ge.npressureslabs) then

     n_in = n_in + 1
     call task_receive_float(buff_slabs(1,1,1,n_in), &
                                nzslab*nxp1*nyp22,reqs_in(n_in))
     flag(n_in) = .false.	    
   endif

end do ! m

! Blocking send now:

do m = 0, 1-1

   call task_rank_to_index(m,it,jt)

   if(0.lt.npressureslabs.and.m.ne.0) then

     do k = 1,nzslab
      do j = 1,nyp22
       do i = 1,nxp1
         buff_subs(i,j,k,1) = f(i+it,j+jt,k)
       end do
      end do
     end do

     call task_bsend_float(m, buff_subs(1,1,1,1),nzslab*nxp1*nyp22,44)

   endif

end do ! m



! Fill slabs when receive buffers are complete:


count = n_in
do while (count .gt. 0)
  do m = 1,n_in
   if(.not.flag(m)) then
	call task_test(reqs_in(m), flag(m), rf, tag)
        if(flag(m)) then 
	   count=count-1
           n = rf*nzslab           
           do k = 1,nzslab
             do j=1,nyp22
               do i=1,nxp1
                 ff(i,j,k+n) = buff_slabs(i,j,k,m)
               end do
             end do
           end do
	endif   
   endif
  end do
end do

!-------------------------------------------------
!   Solve the tri-diagonal system for Fourier coeffiecients 
!   in the vertical for each subdomain:

do k=1,nzm
    a(k)=rhow(k)/rho(k)/(adz(k)*adzw(k)*dz*dz)
    c(k)=rhow(k+1)/rho(k)/(adz(k)*adzw(k+1)*dz*dz)	 
end do 

call task_rank_to_index(0,it,jt)
	
ddx2=1._8/(dx*dx)
ddy2=1._8/(dy*dy)
pii = acos(-1._8)
xnx=pii/nx
xny=pii/ny 	 
do j=1,nyp22-jwall
   if(dowally) then
      jd=j+jt-1
      facty = 1.d0
   else
      jd=(j+jt-0.1)/2.
      facty = 2.d0
   end if
   xj=jd
   do i=1,nxp1-iwall
      if(dowallx) then
        id=i+it-1
        factx = 1.d0
      else
        id=(i+it-0.1)/2.
        factx = 2.d0
      end if
      fff(1:nzm) = ff(i,j,1:nzm)
      xi=id
      eign=(2._8*cos(factx*xnx*xi)-2._8)*ddx2+ & 
            (2._8*cos(facty*xny*xj)-2._8)*ddy2
      if(id+jd.eq.0) then               
         b=1._8/(eign-a(1)-c(1))
         alfa(1)=-c(1)*b
         beta(1)=fff(1)*b
      else
         b=1._8/(eign-c(1))
         alfa(1)=-c(1)*b
         beta(1)=fff(1)*b
      end if
      do k=2,nzm-1
        e=1._8/(eign-a(k)-c(k)+a(k)*alfa(k-1))
        alfa(k)=-c(k)*e
        beta(k)=(fff(k)-a(k)*beta(k-1))*e
      end do

      fff(nzm)=(fff(nzm)-a(nzm)*beta(nzm-1))/ &
	        (eign-a(nzm)+a(nzm)*alfa(nzm-1))
	  
      do k=nzm-1,1,-1
       fff(k)=alfa(k)*fff(k+1)+beta(k)
      end do
      ff(i,j,1:nzm) = fff(1:nzm)

   end do  
end do 

call task_barrier()

!-----------------------------------------------------------------	 
!   Send the Fourier coefficient to the tasks performing
!   the inverse Fourier transformation:

! Non-blocking receive first:

n_in = 0
do m = 0,1-1

  if(0.lt.npressureslabs.and.m.ne.1-1) then
    n_in = n_in + 1
    call task_receive_float(buff_subs(1,1,1,n_in), &
                              nzslab*nxp1*nyp22, reqs_in(n_in))
    flag(n_in) = .false.	    
  endif

  if(0.lt.npressureslabs.and.m.eq.1-1) then

    call task_rank_to_index(0,it,jt)	  
    n = 0*nzslab
    do k = 1,nzslab
     do j = 1,nyp22-jwall
       do i = 1,nxp1-iwall
         f(i+it,j+jt,k) = ff(i,j,k+n)
       end do
     end do
    end do

  endif

end do ! m

! Blocking send now:

do m = 0,1-1

  if(m.lt.npressureslabs.and.m.ne.0) then
    n = m*nzslab+1
    call task_bsend_float(m,ff(1,1,n),nzslab*nxp1*nyp22, 33)
  endif

end do ! m


! Fill slabs when receive buffers are full:


count = n_in
do while (count .gt. 0)
  do m = 1,n_in
   if(.not.flag(m)) then
	call task_test(reqs_in(m), flag(m), rf, tag)
        if(flag(m)) then 
	   count=count-1
           call task_rank_to_index(rf,it,jt)	  
           do k = 1,nzslab
            do j = 1,nyp22-jwall
             do i = 1,nxp1-iwall
                f(i+it,j+jt,k) = buff_subs(i,j,k,m)
             end do
            end do
           end do
	endif   
   endif
  end do
end do

!-------------------------------------------------
!   Perform inverse Fourier transformation:

if(0.lt.npressureslabs) then

 do k=1,nzslab

  if(RUN3D) then
   if(dowally) then
     call cosft_crm(f(1,1,k),work,trigxj,ifaxj,nx2,1,ny,nx+1,+1)
   else
     call fft991_crm(f(1,1,k),work,trigxj,ifaxj,nx2,1,ny,nx+1,+1)
   end if
  end if
	 
  if(dowallx) then
   call cosft_crm(f(1,1,k),work,trigxi,ifaxi,1,nx2,nx,ny,+1)
  else
   call fft991_crm(f(1,1,k),work,trigxi,ifaxi,1,nx2,nx,ny,+1)
  end if

 end do 

endif

call task_barrier()

!-----------------------------------------------------------------	 
!   Fill the pressure field for each subdomain: 

do i=1,nx
 iii(i)=i
end do
iii(0)=nx
do j=1,ny
 jjj(j)=j
end do
jjj(0)=ny

! Non-blocking receive first:

n_in = 0
do m = 0, 1-1
		
   call task_rank_to_index(m,it,jt)

   if(m.lt.npressureslabs-1.or.  &
		m.eq.npressureslabs-1.and.0.ge.npressureslabs) then

     n_in = n_in + 1
     call task_receive_float(bufp_slabs(0,1-YES3D,1,n_in), &
                                nzslab*nxp1*nyp1, reqs_in(n_in))
     flag(n_in) = .false.    

   endif

   if(0.lt.npressureslabs.and.m.eq.0) then

     n = 0*nzslab
     do k = 1,nzslab
      do j = 1-YES3D,ny
       jj=jjj(j+jt)
        do i = 0,nx
	 ii=iii(i+it)
          p(i,j,k+n) = f(ii,jj,k) 
        end do
      end do
     end do 

   end if

end do ! m


! Blocking send now:

do m = 0, 1-1

   call task_rank_to_index(m,it,jt)

   if(0.lt.npressureslabs.and.m.ne.0) then

     do k = 1,nzslab
      do j = 1-YES3D,ny
       jj=jjj(j+jt)
       do i = 0,nx
         ii=iii(i+it)
         bufp_subs(i,j,k,1) = f(ii,jj,k)
       end do
      end do
     end do

     call task_bsend_float(m, bufp_subs(0,1-YES3D,1,1), nzslab*nxp1*nyp1,44)

   endif

end do ! m

! Fill the receive buffers:

count = n_in
do while (count .gt. 0)
  do m = 1,n_in
   if(.not.flag(m)) then
	call task_test(reqs_in(m), flag(m), rf, tag)
        if(flag(m)) then 
	   count=count-1
           n = rf*nzslab           
           do k = 1,nzslab
            do j=1-YES3D,ny
             do i=0,nx
               p(i,j,k+n) = bufp_slabs(i,j,k,m)
             end do
            end do
           end do
        endif   
   endif
  end do
end do


call task_barrier()

!  Add pressure gradient term to the rhs of the momentum equation:

call press_grad()

end subroutine pressure_mpiensemble



