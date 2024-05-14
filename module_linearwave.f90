!---------------------------------------------
!Two subroutines that performs the coupling between a single column model
!(or a limited domain cloud resolving model)
!and a 2D gravity wave with a single horizontal wave number
!Reference: Kuang ,2008, Journal of Atmospheric Sciences, pg 576-591
!Zhiming Kuang 2008
!---------------------------------------------
!There are two subroutines
!1. calc_wtend computes the tendency for vertical velocity, given the current virtual temperature profile.
!2. linearwavesubsidence computes the temperature and moisture tendencies given the current vertical velocity profile

! Editted by Qiyu Song (2024)

module module_linearwave

implicit none
private
real, parameter:: ggr=9.8  !gravity m/s2
real, parameter:: Cp=1004. !specific heat J/kg/K
public:: calc_wtend, linearwavesubsidence

contains

!------------------------------------------------------------------------------
!Compare the current sounding with the reference sounding, from which infer the vertical velocity tendency

subroutine calc_wtend(           wn,  w_curr, &
      tv_curr,        tv_fullbg,       rho_full ,   &
      z_full,                          z_half,      &
      wtend,                           nz)
!     ------------------------------ input arguments ------------------------------

      integer, intent(in) :: nz  ! number of midpoint levels, the number of interface levels is nz+1

      real, dimension(nz), intent(in) ::                 &
      tv_curr,           &       ! Cell center virtual temperature
      w_curr,            &       ! Cell center wave vertical velocity
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
     
      N2top=ggr/tv_fullbg(nz)*((tv_fullbg(nz)-tv_fullbg(nz-1))/(z_full(nz)-z_full(nz-1))+ggr/Cp)

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
      do k=1,nz
         aa(k)=dz(k+1)/(dz(k)+dz(k+1))
         bb(k)=-1.
         cc(k)=dz(k)/(dz(k)+dz(k+1))
      end do

      !symmetric lower BC
      aa(1)=0.
      bb(1)=-(2*dz(2)+dz(1))/(dz(1)+dz(2))

      !radiating upper BC
      aa(nz)=dz(nz+1)/dz(nz)
      bb(nz)=-1.*dz(nz+1)/dz(nz)
      rhs(nz)=rhs(nz)+rho_full(nz)*sqrt(N2top)*wn*w_curr(nz)*dz(nz+1)

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

!------------------------------------------------------------------------------

subroutine linearwavesubsidence(           w_curr,      &
      t_fullbg,                       q_fullbg,    &
      z_full,                          ttend,        &
      qtend,                           nz )

!     ------------------------------ input arguments ------------------------------

      integer, intent(in) :: nz ! number of midpoint levels, the number of interface levels is nz+1

      real, dimension(nz), intent(in) ::                 &
      t_fullbg,         &       ! Cell center background temperature
      q_fullbg,         &       ! Cell center background moisture
      w_curr,           &       ! Cell center wave vertical velocity
      z_full                    ! Cell center height


!     ------------------------------ output arguments -----------------------------

      real, dimension(nz), intent(out) ::                 &
      ttend,            &       ! Cell center T tendency
      qtend                     ! Cell center q tendency

!     ------------------------------ local variables ------------------------------
      real ggr_over_Cp          ! gravity divided specific heat

      integer :: k              ! Vertical loop counter

!     ------------------------------ executable code ------------------------------

      ggr_over_Cp=ggr/Cp

      do k=2,nz-1
         ttend(k) =  - w_curr(k) * ((t_fullbg(k+1)-t_fullbg(k-1))/(z_full(k+1)-z_full(k-1))+ggr_over_Cp)
         qtend(k) =  - w_curr(k) * (q_fullbg(k+1)-q_fullbg(k-1))/(z_full(k+1)-z_full(k-1))
      end do

      ! Use 1st order formular at upper and lower boundaries
      ttend(1)=-w_curr(1)*((t_fullbg(2)-t_fullbg(1))/(z_full(2)-z_full(1))+ggr_over_Cp)
      qtend(1)=-w_curr(1)*(q_fullbg(2)-q_fullbg(1))/(z_full(2)-z_full(1))
      ttend(nz)=-w_curr(nz)*((t_fullbg(nz)-t_fullbg(nz-1))/(z_full(nz)-z_full(nz-1))+ggr_over_Cp)
      qtend(nz)=-w_curr(nz)*(q_fullbg(nz)-q_fullbg(nz-1))/(z_full(nz)-z_full(nz-1))

end subroutine linearwavesubsidence

end module module_linearwave

