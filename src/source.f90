!
! SPDX-License-Identifier: MIT
!
module mod_source
  !
  use mod_types
  use mod_common_mpi, only: myid,ierr,ijk_start
  use mod_param     , only: sigma,gacc,dpdl
  !
  implicit none
  !
  private
#if defined(_USE_VOF)
  public :: surft_src
#endif
  public :: grav_src,pres_src
#if defined(_TURB_FORCING)
  public :: forc_src
#endif
  !
  contains
  !
#if defined(_USE_VOF)
  subroutine surft_src(nx,ny,nz,dxi,dyi,dzi,kappa,psi,rho,dudt,dvdt,dwdt)
    !
    implicit none 
    !
    integer , intent(in   )                      :: nx,ny,nz
    real(rp), intent(in   )                      :: dxi,dyi,dzi
    real(rp), intent(in   ), dimension(0:,0:,0:) :: kappa,psi,rho
    real(rp), intent(inout), dimension(1:,1:,1:) :: dudt,dvdt,dwdt
    !
    real(rp) :: kappax,kappay,kappaz
    real(rp) :: rhox,rhoy,rhoz
    integer  :: i,j,k,ip,jp,kp
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: kappa, psi, rho, dudt, dvdt, dwdt
#endif
    !
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1   
          jp = j+1   
          kp = k+1  
          !
          rhox   = 0.5_rp*(rho(ip,j,k)+rho(i,j,k))
          rhoy   = 0.5_rp*(rho(i,jp,k)+rho(i,j,k))
          rhoz   = 0.5_rp*(rho(i,j,kp)+rho(i,j,k))
          kappax = 0.5_rp*(kappa(ip,j,k)+kappa(i,j,k))
          kappay = 0.5_rp*(kappa(i,jp,k)+kappa(i,j,k))
          kappaz = 0.5_rp*(kappa(i,j,kp)+kappa(i,j,k))
          !
          dudt(i,j,k) = dudt(i,j,k) + dxi*sigma*kappax*(psi(ip,j,k)-psi(i,j,k))/rhox
          dvdt(i,j,k) = dvdt(i,j,k) + dyi*sigma*kappay*(psi(i,jp,k)-psi(i,j,k))/rhoy
          dwdt(i,j,k) = dwdt(i,j,k) + dzi*sigma*kappaz*(psi(i,j,kp)-psi(i,j,k))/rhoz
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine surft_src
#endif
  !
  subroutine grav_src(nx,ny,nz,rho,rho_av, &
#if defined(_BOUSSINESQ)
                      nh_t,tmp, & 
#endif
                      dudt,dvdt,dwdt)  
    !
#if defined(_BOUSSINESQ)
    use mod_param     , only: tmp0,beta_th
#endif
    !
    implicit none
    ! 
    integer , intent(in   )                                     :: nx,ny,nz
    real(rp), intent(in   ), dimension(     0:,     0:,     0:) :: rho
    real(rp), intent(in   )                                     :: rho_av
#if defined(_BOUSSINESQ)
    integer , intent(in   )                                     :: nh_t
    real(rp), intent(in   ), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: tmp
#endif
    real(rp), intent(inout), dimension(     1:,     1:,     1:) :: dudt,dvdt,dwdt
    !
    real(rp) :: rhox,rhoy,rhoz, gacc1, gacc2, gacc3 
    real(rp) :: termx,termy,termz,tmppx,tmppy,tmppz
    integer  :: i,j,k,ip,jp,kp
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: rho, dudt, dvdt, dwdt
#if defined(_BOUSSINESQ)
    attributes(managed) :: tmp
#endif
#endif
    !
    ! in case it is periodic, subtract mean gravitational force per unit mass
    !
    gacc1 = gacc(1)
    gacc2 = gacc(2)
    gacc3 = gacc(3)
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1
          jp = j+1
          kp = k+1
          !
          rhox = 0.5_rp*(rho(ip,j,k)+rho(i,j,k))
          rhoy = 0.5_rp*(rho(i,jp,k)+rho(i,j,k))
          rhoz = 0.5_rp*(rho(i,j,kp)+rho(i,j,k))
          !
#if defined(_BOUSSINESQ)
          tmppx = 0.5_rp*(tmp(ip,j,k)+tmp(i,j,k))
          tmppy = 0.5_rp*(tmp(i,jp,k)+tmp(i,j,k))
          tmppz = 0.5_rp*(tmp(i,j,kp)+tmp(i,j,k))
          termx = gacc1*((1.0_rp-beta_th*(tmppx-tmp0))-rho_av/rhox)
          termy = gacc2*((1.0_rp-beta_th*(tmppy-tmp0))-rho_av/rhoy)
          termz = gacc3*((1.0_rp-beta_th*(tmppz-tmp0))-rho_av/rhoz)
#else
          termx = gacc1*(rhox-rho_av)/rhox
          termy = gacc2*(rhoy-rho_av)/rhoy
          termz = gacc3*(rhoz-rho_av)/rhoz
#endif
          dudt(i,j,k) = dudt(i,j,k) + termx
          dvdt(i,j,k) = dvdt(i,j,k) + termy
          dwdt(i,j,k) = dwdt(i,j,k) + termz
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine grav_src
  !
  subroutine pres_src(nx,ny,nz,dxi,dyi,dzi,rho0i,rho,p,pp,dudt,dvdt,dwdt)
    !
    implicit none 
    !
    integer , intent(in )                      :: nx,ny,nz
    real(rp), intent(in )                      :: dxi,dyi,dzi
    real(rp), intent(in )                      :: rho0i
    real(rp), intent(in ), dimension(0:,0:,0:) :: p,pp,rho
    real(rp), intent(out), dimension(1:,1:,1:) :: dudt,dvdt,dwdt
    !
    real(rp) :: rhox,rhoy,rhoz,rhoxi,rhoyi,rhozi
    integer  :: i,j,k,ip,jp,kp
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: p, pp, rho, dudt,dvdt, dwdt
#endif
    !
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1
          jp = j+1
          kp = k+1
          !
          rhox  = 0.5_rp*(rho(ip,j,k)+rho(i,j,k))
          rhoxi = 1.0_rp / rhox
          !
          rhoy  = 0.5_rp*(rho(i,jp,k)+rho(i,j,k))
          rhoyi = 1.0_rp / rhoy
          !
          rhoz  = 0.5_rp*(rho(i,j,kp)+rho(i,j,k))
          rhozi = 1.0_rp / rhoz
          !
          dudt(i,j,k) = ( - ( p(ip,j,k)-p(i,j,k) )*dxi - dpdl(1))*rho0i &
                          - (rhoxi  - rho0i)*(pp(ip,j,k)-pp(i,j,k))*dxi
          dvdt(i,j,k) = ( - ( p(i,jp,k)-p(i,j,k) )*dyi - dpdl(2))*rho0i &
                          - (rhoyi  - rho0i)*(pp(i,jp,k)-pp(i,j,k))*dyi
          dwdt(i,j,k) = ( - ( p(i,j,kp)-p(i,j,k) )*dzi - dpdl(3))*rho0i &
                          - (rhozi  - rho0i)*(pp(i,j,kp)-pp(i,j,k))*dzi
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine pres_src
  !
#if defined(_TURB_FORCING)
  subroutine forc_src(nx,ny,nz,dx,dy,dz,dudt,dvdt,dwdt)
    !
    use mod_param, only: f0_t,k0_t,turbType,abc,pi,lx,ly,lz
    !
    implicit none
    ! 
    integer , intent(in   )                      :: nx,ny,nz
    real(rp), intent(in   )                      :: dx,dy,dz
    real(rp), intent(inout), dimension(1:,1:,1:) :: dudt,dvdt,dwdt
    !
    real(rp) :: xc,yc,zc,xf,yf
    real(rp) :: abc1,abc2,abc3
    integer  :: i,j,k
#if defined(_OPENACC)
    attributes(managed) :: dudt,dvdt,dwdt
    integer  :: istat
#endif
    !
    abc1 = abc(1)
    abc2 = abc(2)
    abc3 = abc(3)
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          zc = (k+ijk_start(3)-0.5_rp)*dz/lz*2._rp*pi
          yc = (j+ijk_start(2)-0.5_rp)*dy/ly*2._rp*pi
          yf = (j+ijk_start(2)-0.0_rp)*dy/ly*2._rp*pi
          xc = (i+ijk_start(1)-0.5_rp)*dx/lx*2._rp*pi
          xf = (i+ijk_start(1)-0.0_rp)*dx/lx*2._rp*pi
          !
          if(    turbType.eq.'tgv') then
            dudt(i,j,k) = dudt(i,j,k) + f0_t*sin(k0_t*xf)*cos(k0_t*yc)*cos(k0_t*zc)
            dvdt(i,j,k) = dvdt(i,j,k) - f0_t*cos(k0_t*xc)*sin(k0_t*yf)*cos(k0_t*zc)
          elseif(turbType.eq.'abc') then
            dudt(i,j,k) = dudt(i,j,k) + f0_t*(abc1*sin(k0_t*zc) + abc3*cos(k0_t*yc))
            dvdt(i,j,k) = dvdt(i,j,k) + f0_t*(abc2*sin(k0_t*xc) + abc1*cos(k0_t*zc))
            dwdt(i,j,k) = dwdt(i,j,k) + f0_t*(abc3*sin(k0_t*yc) + abc2*cos(k0_t*xc))
          endif
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine forc_src
#endif
  !
end module mod_source
