! This subroutine solves for perturbation omega based on WTG approximation as
! described in equation (6) of Kuang (2008):
!
!   http://dx.doi.org/10.1175/2007JAS2399.1

subroutine wtg_jas2008(nzm, dtn, z, zi, rho, tabs_ref, qv_ref, tabs_model, &
     qv_model, qcond_model, lambda_wtg, am_wtg, w_ls, dwdt_ls)

  use grid, only: icycle, ncycle, nsubdomains
  use params, only: dowtg_timedependence, dompiensemble, cp, ggr, pi
  !use vars
  implicit none

  ! ======= inputs =======
  integer, intent(in) :: nzm ! number of model levels
  real, intent(in) :: dtn ! current dynamical time step
  real, intent(in) :: z(nzm) ! model cell center height
  real, intent(in) :: zi(nzm+1) ! model interface height
  real, intent(in) :: rho(nzm) ! model cell center density
  real, intent(in) :: tabs_ref(nzm) ! reference temperature sounding in K
  real, intent(in) :: qv_ref(nzm)   ! reference water vapor sounding in kg/kg
  real, intent(in) :: tabs_model(nzm) ! model temperature sounding in K (domain-mean for LES)
  real, intent(in) :: qv_model(nzm)   ! model water vapor sounding in kg/kg (domain-mean for LES)
  real, intent(in) :: qcond_model(nzm) ! model condensate (liq+ice) sounding in kg/kg (domain-mean for LES)

  real, intent(in) :: lambda_wtg ! WTG length scale in m (JAMES2009 value = 650.e6 m)

  ! WTG momentum damping rate
  real, intent(in) :: am_wtg     ! WTG momentum damping rate in 1/s (default = 1./86400. /s)

  ! ======= input+output =======
  real, intent(inout) :: w_ls(nzm) ! WTG large-scale vertical velocity in m/s on model levels.

  ! ======= output =======
  real, intent(out) :: dwdt_ls(nzm) ! WTG large-scale vertical velocity tendency on model levels.

  ! ======= local variables =======
   integer :: k
   integer :: ktrop
   integer :: ktrop1
   real :: ztrop
   real :: ztrop1
   real :: min_temp
   real :: coef, buffer(nzm), buffer1(nzm)
   real :: tv_lsbg(nzm)
   real :: tv_wtg(nzm)

   if (dowtg_timedependence) then

      ! ! ===== find index of cold point tropopause in vertical. =====
      ! ! reverse pressure coordinate, and find index
      ! !   of cold point tropopause in the vertical.
      ! ktrop = nzm+1 ! default is top of model/atmosphere (counting from surface)
      ! min_temp = tg0(nzm)
      ! do k = 1,nzm
      !    if(tg0(k).lt.min_temp) then
      !       ktrop = k
      !       min_temp = tg0(k)
      !    end if
      ! end do

      ! ! ===== find index of previous cold point tropopause in vertical. =====
      ! ! reverse pressure coordinate, and find index
      ! !   of cold point tropopause in the vertical.
      ! ktrop1 = nzm !
      ! do k = nzm,1,-1
      !    if(w_ls(k).EQ.0) then
      !       ktrop1 = k
      !    end if
      ! end do

      ! if ((ktrop1.NE.1).AND.(ktrop1.NE.ktrop))

      !    ztrop = z(ktrop)
      !    ztrop1 = z(ktrop1)
      !    wwtgi = 0.

      !    do k = 2 : (ktrop-1)
      !       do kk = 1 : nzm
      !          wwtgi(k)
      !       end do
      !    end do

      !    do k = 1,nzm
      !       w_ls(k) = wwtgi(k)
      !    end do

      ! end if

      ! ! At tropopause and above, quash any pre-existing vertical velocity from previous timestep to zero.
      ! do k = ktrop,nzm
      !    w_ls(k) = 0
      ! end do

      ktrop = nzm-2
   
   else

      ! ===== find index of cold point tropopause in vertical. =====
      ! reverse pressure coordinate, and find index
      !   of cold point tropopause in the vertical.
      ktrop = nzm+1 ! default is top of model/atmosphere (counting from surface)
      min_temp = tabs_model(nzm)
      do k = 1,nzm
         if(tabs_model(k).lt.min_temp) then
            ktrop = k
            min_temp = tabs_model(k)
         end if
      end do

   end if

   tv_lsbg = tabs_ref   * (1. + 0.61*qv_ref)
   tv_wtg = tabs_model * (1. + 0.61*qv_model - qcond_model)

   dwdt_ls = 0.
   
   call calc_wtend(0.5*pi/lambda_wtg, w_ls(1:ktrop), dwdt_ls(1:ktrop), &
                     tv_wtg(1:ktrop), tv_lsbg(1:ktrop), &
                     rho(1:ktrop), z(1:ktrop), zi(1:ktrop+1), ktrop)

   if (dowtg_timedependence) then

      w_ls(1:nzm) = (w_ls(1:nzm) + dwdt_ls * dtn) / (1. + dtn * am_wtg)

   else

      w_ls = 0.
      w_ls(1:nzm) = dwdt_ls / am_wtg

   end if

   contains

   subroutine calc_wtend(wn, w_curr, wtend, &
                           tv_curr, tv_fullbg, &
                           rho_full, z_full, z_half, nz)

   !     ------------------------------ input arguments ------------------------------

      integer, intent(in) :: nz  ! number of midpoint levels, the number of interface levels is nz+1

      real, dimension(nz), intent(in) ::                 &
      w_curr,            &       ! Cell center wave vertical velocity
      tv_curr,           &       ! Cell center virtual temperature
      tv_fullbg,         &       ! Cell center background virtual temperature
      z_full,            &       ! Cell center height
      rho_full                   ! Cell center density
      real, dimension(nz+1), intent(in) ::                 &
      z_half                     ! Interface-level height, including the surface, which is at height zero

      real, intent (in) :: wn    !the specified horizontal wavenumber

!     ------------------------------ output arguments -----------------------------

      real, dimension(nz), intent(out) ::                 &
      wtend                      ! Cell center w tendency

!     ------------------------------ local variables ------------------------------

      real  N2top                ! Top-level buoyancy frequency squared
      real, dimension(1:nz+1) ::                                             &
      dz                         ! grid spacing between midpoint levels
      real, dimension(1:nz) ::                                             &
      rhs,               &       ! The right hand side for the Gauss Elimination
      aa,                &       ! aa, bb, cc are the three columns of the tridiagonal matrix
      bb,                &
      cc

      integer :: k               ! Vertical loop counter

!     ------------------------------ executable code ------------------------------

      N2top=ggr/tv_fullbg(nz)*((tv_fullbg(nz)-tv_fullbg(nz-1))/(z_full(nz)-z_full(nz-1))+ggr/cp)

!     compute grid spacing between midpoint levels
      do k=2,nz
         dz(k)=z_full(k)-z_full(k-1)
      end do
      dz(1)=2*z_full(1)
      dz(nz+1)=2*(z_half(nz+1)-z_full(nz))

!     Gaussian Elimination
      rhs=0.
      do k=1,nz
         rhs(k)=-rho_full(k)*ggr*wn*wn*(tv_curr(k)-tv_fullbg(k))/tv_fullbg(k)*dz(k)*dz(k+1)/2.
      end do

!     set up the tridiagonal matrix
      do k=1,nz-1
         aa(k)=dz(k+1)/(dz(k)+dz(k+1))
         bb(k)=-1.
         cc(k)=dz(k)/(dz(k)+dz(k+1))
      end do

      ! Symmetric boundary conditions?
      aa(1)=0.
      bb(1)=-(2*dz(2)+dz(1))/(dz(1)+dz(2))

      if (dowtg_timedependence) then

         !radiating upper BC
         aa(nz)=dz(nz+1)/dz(nz)
         bb(nz)=-1.*dz(nz+1)/dz(nz)
         rhs(nz)=rhs(nz)+rho_full(nz)*sqrt(N2top)*wn*w_curr(nz)*dz(nz+1)

      else

         ! dirichlet boundary conditions at tropopause (w' = 0)
         aa(nz)=0.
         bb(nz)=1.
         cc(nz)=0.
         rhs(nz)=0

      end if

      !Gaussian Elimination with no pivoting
      do k=1,nz-1
         bb(k+1)=bb(k+1)-aa(k+1)/bb(k)*cc(k)
         rhs(k+1)=rhs(k+1)-aa(k+1)/bb(k)*rhs(k)
      end do

      !backward substitution
      rhs(nz)=rhs(nz)/bb(nz)
      do k= nz-1,1,-1
         rhs(k)=(rhs(k)-cc(k)*rhs(k+1))/bb(k)
      end do

      do k=1,nz
         wtend(k)=rhs(k)/rho_full(k)
      end do

   end subroutine calc_wtend

end subroutine wtg_jas2008
