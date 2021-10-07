!
! SPDX-License-Identifier: MIT
!
module mod_funcs
  !
  use mod_param, only: pi
  use mod_types, only: rp
  !
  implicit none
  !
  private
  public ssign,heaviside,dirac,efun,interp_d,interp_g
  !
  contains
  !
  function ssign(delta,phi)
    !
    ! smooth sign function
    !
    implicit none
    !
    real(rp), intent(in) :: delta,phi
    real(rp)             :: ssign
    !
    ssign = sign(1._rp,phi)
    if(abs(phi).le.delta) then
      ssign = phi/(sqrt(phi**2+delta**2))
    endif
    !
    return
  end function ssign
  !
  function heaviside(r,eps)
    !
    ! smooth step function based
    ! on cosine/sine
    !
    implicit none
    !
    real(rp), intent(in) :: r,eps
    real(rp)             :: heaviside
    !
    if(r.lt.-eps) then
      heaviside = 0._rp
    elseif(r.lt.eps) then
      heaviside = 0.5_rp + 0.5_rp*r/eps + 0.5_rp/pi*sin(pi*r/eps)
    else
      heaviside = 1._rp
    endif
    !
    return
  end function heaviside
  !
  function dirac(r,eps)
    !
    ! smooth impulse function
    ! (note: the derivative of heaviside function w.r.t. eps)
    !
    implicit none
    !
    real(rp), intent(in) :: r,eps
    real(rp)             :: dirac
    !
    if(abs(r).ge.eps) then
      dirac = 0._rp
    else
      dirac = 0.5_rp/eps + 0.5_rp/eps*cos(pi*r/eps)
    endif
    !
    return
  end function dirac
  !
  function efun(r,eps)
    !
    ! smooth step function based 
    ! on the error function
    !
    implicit none
    !
    real(rp), intent(in) :: r,eps
    real(rp)             :: efun
    !
    efun = 0.5_rp*( 1._rp+erf(r/eps) ) 
    !
    return
  end function efun
  !
  function dir_efun(r,eps)
    !
    ! smooth impulse function based 
    ! on the error function
    ! (note: the derivative of error function w.r.t. eps)
    !
    implicit none
    !
    real(rp), intent(in) :: r,eps
    real(rp)             :: dir_efun
    !
    dir_efun = (1._rp/(eps*pi))*exp(-(r/eps)**2._rp)
    !
    return
  end function dir_efun
  !
  function interp_d(vec,vel,q) result(flux) ! divergence form
    !
    ! QUICK discretization (by Leonard) in divergence form 
    ! (see, e.g.: Numerical simulation of incompressible flows, pag. 112)
    !
    implicit none
    !
    real(rp), intent(in), dimension(-2:2) :: vec
    real(rp), intent(in) :: vel ! the transporting velocity (in d(uu)/dx is u, in d(uv)/dy is v)
    integer , intent(in) :: q   ! (q=-1 if f_i, q=0 if f_i+1 with f is the flux)
    !
    real(rp) :: flux
    !
    !flux = 0.5_rp*vel*(vec(0+q)+vec(1+q)) + 0.5_rp*abs(vel)*(+vec(0+q)-vec(1+q))
    !
    flux =      0.0625_rp*(vel)*(-vec(-1+q)+9.0_rp*vec(0+q)+9.0_rp*vec(1+q)-vec(2+q)) &
           + 0.0625_rp*abs(vel)*(-vec(-1+q)+3.0_rp*vec(0+q)-3.0_rp*vec(1+q)+vec(2+q)) 
    !
    return
  end function interp_d
  !
  function interp_g(vec,vel,dli) result(vel_grad) ! gradient form
    !
    ! II order (Kawamura and Kuwahara equation) discretization in gradient form 
    ! (see, e.g.: Numerical simulation of incompressible flows, pag. 113)
    !
    implicit none
    !
    real(rp), intent(in), dimension(-2:2) :: vec
    real(rp), intent(in)                  :: vel ! the transporting velocity (in ud(f)/dx is u, in vd(f)/dy is v)
    real(rp), intent(in)                  :: dli  
    !
    real(rp) :: vel_grad
    !
    !vel_grad = 0.5_rp*(vel+abs(vel))*(vec(0)-vec(-1))*dli + 0.5_rp*(vel-abs(vel))*(vec(1)-vec(0))*dli
    !
    vel_grad =             (vel)*(vec(-2)-8.0_rp*vec(-1)               +8.0_rp*vec(+1)-vec(+2))*(1.0_rp/12.0_rp)*dli & 
               + 3.0_rp*abs(vel)*(vec(-2)-4.0_rp*vec(-1)+6.0_rp*vec(+0)-4.0_rp*vec(+1)+vec(+2))*(1.0_rp/12.0_rp)*dli
    !
    return
  end function interp_g
  !
end module mod_funcs 
