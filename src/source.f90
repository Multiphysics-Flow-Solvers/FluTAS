!
! SPDX-License-Identifier: MIT
!
module mod_source
  !
  use mod_types
  use mod_common_mpi, only: myid,ierr,ijk_start
  use mod_sanity    , only: flutas_error
  !
  implicit none
  !
  private
  public :: bulk_forcing_src
#if defined(_USE_VOF)
  public :: surft_src
  public :: grav_tw_src,pres_tw_src
#else
  public :: grav_sp_src,pres_sp_src
#endif
#if defined(_TURB_FORCING)
  public :: forc_src
#endif
  !
  contains
  !
#if defined(_USE_VOF)
  subroutine surft_src(nx,ny,nz,nh_d,nh_u,f_t12,dxi,dyi,dzi,dzci,kappa,psi,rho,u,v,w)
    !
    use mod_param, only: sigma
    !
    implicit none 
    !
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   )                                     :: nh_d,nh_u
    real(rp), intent(in   )                                     :: f_t12
    real(rp), intent(in   )                                     :: dxi,dyi,dzi
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzci
    real(rp), intent(in   ), dimension(0:,0:,0:)                :: kappa,psi,rho
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp) :: kappax,kappay,kappaz
    real(rp) :: rhox,rhoy,rhoz
    integer  :: i,j,k,ip,jp,kp
    !@cuf attributes(managed) :: kappa, psi, rho, dzci, u, v, w
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
          u(i,j,k) = u(i,j,k) + f_t12*(sigma*kappax*(psi(ip,j,k)-psi(i,j,k))/rhox)*dxi
          v(i,j,k) = v(i,j,k) + f_t12*(sigma*kappay*(psi(i,jp,k)-psi(i,j,k))/rhoy)*dyi
          w(i,j,k) = w(i,j,k) + f_t12*(sigma*kappaz*(psi(i,j,kp)-psi(i,j,k))/rhoz)*dzci(k)
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
#if defined(_USE_VOF)
  subroutine grav_tw_src(nx,ny,nz,rho,f_t12,cbcpre,dxi,dyi,dzi,nh_d,dzfi, &
#if defined(_BOUSSINESQ)
                         vof,nh_t,tmp, & 
#endif
                         nh_u,u,v,w)  
    !
    use mod_param     , only: lx,ly,lz,gacc_x,gacc_y,gacc_z
#if defined(_BOUSSINESQ)
    use mod_param     , only: tmp0,beta1_th,beta2_th,rho1,rho2
#endif
    !
    implicit none
    ! 
    integer         , intent(in   )                                     :: nx,ny,nz
    real(rp)        , intent(in   ), dimension(     0:,     0:,     0:) :: rho
    real(rp)        , intent(in   )                                     :: f_t12
    character(len=1), intent(in   ), dimension(0:1,3)                   :: cbcpre
    real(rp)        , intent(in   )                                     :: dxi,dyi,dzi
    integer         , intent(in   )                                     :: nh_d
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzfi
#if defined(_BOUSSINESQ)
    real(rp)        , intent(in   ), dimension(     0:,     0:,     0:) :: vof
    integer         , intent(in   )                                     :: nh_t
    real(rp)        , intent(in   ), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: tmp
    !@cuf attributes(managed) :: vof,tmp
#endif
    integer         , intent(in   )                                     :: nh_u
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp) :: rho_av,rhox,rhoy,rhoz,vofpx,vofpy,vofpz
    real(rp) :: termx,termy,termz
    real(rp) :: tmppx,tmppy,tmppz
    integer  :: i,j,k,ip,jp,kp
    !@cuf attributes(managed) :: rho, dzfi, u, v, w
    !
    ! in case the gravity force acts along a periodic direction
    ! calculate average density (to subtract net gravitational force per unit mass)
    ! 
    rho_av = 0._rp
    if(((cbcpre(0,1)//cbcpre(1,1).eq.'PP').and.gacc_x.ne.0._rp).or. &
       ((cbcpre(0,2)//cbcpre(1,2).eq.'PP').and.gacc_y.ne.0._rp).or. &
       ((cbcpre(0,3)//cbcpre(1,3).eq.'PP').and.gacc_z.ne.0._rp)     ) then
#if defined(_OPENACC)
      !$acc parallel loop reduction(+:rho_av)
#endif
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rho_av = rho_av + rho(i,j,k)/(dxi*dyi*dzfi(k))
          enddo
        enddo
      enddo
#if defined(_OPENACC)
      !$acc end parallel loop 
#endif
      call MPI_ALLREDUCE(MPI_IN_PLACE,rho_av,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      rho_av = rho_av/(lx*ly*lz)
    endif
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
          rhox = 0.5_rp*(rho(ip,j,k)+rho(i,j,k))
          rhoy = 0.5_rp*(rho(i,jp,k)+rho(i,j,k))
          rhoz = 0.5_rp*(rho(i,j,kp)+rho(i,j,k))
          !
#if defined(_BOUSSINESQ)
          tmppx = 0.5_rp*(tmp(ip,j,k)+tmp(i,j,k))
          tmppy = 0.5_rp*(tmp(i,jp,k)+tmp(i,j,k))
          tmppz = 0.5_rp*(tmp(i,j,kp)+tmp(i,j,k))
          !
          vofpx = 0.5_rp*(vof(ip,j,k)+vof(i,j,k))
          vofpy = 0.5_rp*(vof(i,jp,k)+vof(i,j,k))
          vofpz = 0.5_rp*(vof(i,j,kp)+vof(i,j,k))
          !
          termx = rhox - (rho1*beta1_th*vofpx+rho2*beta2_th*(1._rp-vofpx))*(tmppx-tmp0)
          termy = rhoy - (rho1*beta1_th*vofpy+rho2*beta2_th*(1._rp-vofpy))*(tmppy-tmp0)
          termz = rhoz - (rho1*beta1_th*vofpz+rho2*beta2_th*(1._rp-vofpz))*(tmppz-tmp0)
#else
          termx = rhox
          termy = rhoy
          termz = rhoz
#endif
          !
          u(i,j,k) = u(i,j,k) + f_t12*gacc_x*(termx-rho_av)/rhox
          v(i,j,k) = v(i,j,k) + f_t12*gacc_y*(termy-rho_av)/rhoy
          w(i,j,k) = w(i,j,k) + f_t12*gacc_z*(termz-rho_av)/rhoz
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine grav_tw_src
#else
  subroutine grav_sp_src(nx,ny,nz,f_t12,cbcpre,dxi,dyi,dzi,nh_d,dzfi, &
#if defined(_BOUSSINESQ)
                         nh_t,tmp, & 
#endif
                         nh_u,u,v,w)  
    !
    use mod_param     , only: lx,ly,lz,gacc_x,gacc_y,gacc_z
#if defined(_BOUSSINESQ)
    use mod_param     , only: tmp0,beta_sp_th
#endif
    !
    implicit none
    ! 
    integer         , intent(in   )                                     :: nx,ny,nz
    real(rp)        , intent(in   )                                     :: f_t12
    character(len=1), intent(in   ), dimension(0:1,3)                   :: cbcpre
    real(rp)        , intent(in   )                                     :: dxi,dyi,dzi
    integer         , intent(in   )                                     :: nh_d
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzfi
#if defined(_BOUSSINESQ)
    integer         , intent(in   )                                     :: nh_t
    real(rp)        , intent(in   ), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: tmp
    !@cuf attributes(managed) :: tmp
#endif
    integer         , intent(in   )                                     :: nh_u
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp) :: termx,termy,termz
    real(rp) :: tmppx,tmppy,tmppz
    integer  :: i,j,k,ip,jp,kp
    !@cuf attributes(managed) :: dzfi, u, v, w
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
#if defined(_BOUSSINESQ)
          tmppx = 0.5_rp*(tmp(ip,j,k)+tmp(i,j,k))
          tmppy = 0.5_rp*(tmp(i,jp,k)+tmp(i,j,k))
          tmppz = 0.5_rp*(tmp(i,j,kp)+tmp(i,j,k))
          !
          termx = 1._rp - beta_sp_th*(tmppx-tmp0)
          termy = 1._rp - beta_sp_th*(tmppy-tmp0)
          termz = 1._rp - beta_sp_th*(tmppz-tmp0)
#else
          termx = 1._rp
          termy = 1._rp
          termz = 1._rp
#endif
          !
          u(i,j,k) = u(i,j,k) + f_t12*gacc_x*termx
          v(i,j,k) = v(i,j,k) + f_t12*gacc_y*termy
          w(i,j,k) = w(i,j,k) + f_t12*gacc_z*termz
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine grav_sp_src
#endif
  !
#if defined(_USE_VOF)
  subroutine pres_tw_src(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,rho0i,f_t12,f_t12_o,p,pold,rho,u,v,w)
    !
    implicit none 
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzci
    real(rp), intent(in )                                     :: rho0i,f_t12,f_t12_o
    real(rp), intent(in ), dimension(0:,0:,0:)                :: p,pold
    real(rp), intent(in ), dimension(0:,0:,0:)                :: rho
    real(rp), intent(out), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp) :: rhox,rhoy,rhoz,rhoxi,rhoyi,rhozi
#if defined(_CONSTANT_COEFFS_POISSON)
    real(rp) :: f1,f2
#endif
    integer  :: i,j,k,ip,jp,kp
    !@cuf attributes(managed) :: dzci, p, pold, rho, u, v, w
    !
#if defined(_CONSTANT_COEFFS_POISSON)
    f1 = 1._rp+(f_t12/f_t12_o)
    f2 =       (f_t12/f_t12_o)
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
#if defined(_CONSTANT_COEFFS_POISSON)
          u(i,j,k) = u(i,j,k) + f_t12*(( - ( p(ip,j,k)-p(i,j,k) )*dxi     )*rho0i &
                              - (rhoxi  - rho0i)*( (f1*p(ip,j,k)-f2*pold(ip,j,k))-(f1*p(i,j,k)-f2*pold(i,j,k)) )*dxi)
          v(i,j,k) = v(i,j,k) + f_t12*(( - ( p(i,jp,k)-p(i,j,k) )*dyi     )*rho0i &
                              - (rhoyi  - rho0i)*( (f1*p(i,jp,k)-f2*pold(i,jp,k))-(f1*p(i,j,k)-f2*pold(i,j,k)) )*dyi)
          w(i,j,k) = w(i,j,k) + f_t12*(( - ( p(i,j,kp)-p(i,j,k) )*dzci(k) )*rho0i &
                              - (rhozi  - rho0i)*( (f1*p(i,j,kp)-f2*pold(i,j,kp))-(f1*p(i,j,k)-f2*pold(i,j,k)) )*dzci(k))
#else
          u(i,j,k) = u(i,j,k) + f_t12*(( - ( p(ip,j,k)-p(i,j,k) )*dxi     )*rhoxi
          v(i,j,k) = v(i,j,k) + f_t12*(( - ( p(i,jp,k)-p(i,j,k) )*dyi     )*rhoyi
          w(i,j,k) = w(i,j,k) + f_t12*(( - ( p(i,j,kp)-p(i,j,k) )*dzci(k) )*rhozi
#endif
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine pres_tw_src
#else
  subroutine pres_sp_src(nx,ny,nz,f_t12,dxi,dyi,dzi,nh_d,nh_u,dzci,rho0i,pold,u,v,w)
    !
    implicit none 
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: f_t12
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzci
    real(rp), intent(in )                                     :: rho0i
    real(rp), intent(in ), dimension(0:,0:,0:)                :: pold
    real(rp), intent(out), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    integer :: i,j,k,ip,jp,kp
    !@cuf attributes(managed) :: dzci, pold, u, v, w
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
          u(i,j,k) = u(i,j,k) + f_t12*( - ( pold(ip,j,k)-pold(i,j,k) )*dxi     )*rho0i
          v(i,j,k) = v(i,j,k) + f_t12*( - ( pold(i,jp,k)-pold(i,j,k) )*dyi     )*rho0i
          w(i,j,k) = w(i,j,k) + f_t12*( - ( pold(i,j,kp)-pold(i,j,k) )*dzci(k) )*rho0i
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine pres_sp_src
#endif
  !
#if defined(_TURB_FORCING)
  subroutine forc_src(nx,ny,nz,nh_d,nh_u,f_t12,dx,dy,dz,zc,u,v,w)
    !
    use mod_param, only: abc_x,abc_y,abc_z,lx,ly,lz,f0_t,k0_t,pi,turb_type
    !
    implicit none
    ! 
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   )                                     :: nh_d,nh_u
    real(rp), intent(in   )                                     :: f_t12
    real(rp), intent(in   )                                     :: dx,dy,dz
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: zc
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp) :: xcc,ycc,zcc,xff,yff
    integer  :: i,j,k,ijk_start_x,ijk_start_y,ijk_start_z
    !@cuf attributes(managed) :: u, v, w, zc
    !
    ijk_start_x = ijk_start(1)
    ijk_start_y = ijk_start(2)
    ijk_start_z = ijk_start(3)
    !
    select case(turb_type)
    case('tgv')
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            !
            !zcc = zc(k)/lz*2._rp*pi
            zcc = (k+ijk_start_z-0.5_rp)*dz/lz*2._rp*pi ! in triperiodic we typically use unstretched grid
            ycc = (j+ijk_start_y-0.5_rp)*dy/ly*2._rp*pi
            yff = (j+ijk_start_y-0.0_rp)*dy/ly*2._rp*pi
            xcc = (i+ijk_start_x-0.5_rp)*dx/lx*2._rp*pi
            xff = (i+ijk_start_x-0.0_rp)*dx/lx*2._rp*pi
            !
            u(i,j,k) = u(i,j,k) + f_t12*( + f0_t*sin(k0_t*xff)*cos(k0_t*ycc)*cos(k0_t*zcc) )
            v(i,j,k) = v(i,j,k) + f_t12*( - f0_t*cos(k0_t*xcc)*sin(k0_t*yff)*cos(k0_t*zcc) )
            !
          enddo
        enddo
      enddo
      !$acc end kernels
    case('abc')
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            !
            !zcc = zc(k)/lz*2._rp*pi
            zcc = (k+ijk_start_z-0.5_rp)*dz/lz*2._rp*pi ! in triperiodic we typically use unstretched grid
            ycc = (j+ijk_start_y-0.5_rp)*dy/ly*2._rp*pi
            yff = (j+ijk_start_y-0.0_rp)*dy/ly*2._rp*pi
            xcc = (i+ijk_start_x-0.5_rp)*dx/lx*2._rp*pi
            xff = (i+ijk_start_x-0.0_rp)*dx/lx*2._rp*pi
            !
            u(i,j,k) = u(i,j,k) + f_t12*( + f0_t*(abc_x*sin(k0_t*zcc) + abc_z*cos(k0_t*ycc)) )
            v(i,j,k) = v(i,j,k) + f_t12*( + f0_t*(abc_y*sin(k0_t*xcc) + abc_x*cos(k0_t*zcc)) )
            w(i,j,k) = w(i,j,k) + f_t12*( + f0_t*(abc_z*sin(k0_t*ycc) + abc_y*cos(k0_t*xcc)) )
            !
          enddo
        enddo
      enddo
      !$acc end kernels
    case default
      call flutas_error('Invalid forcing for HIT - check forcing.in')
    end select
    !
    return
  end subroutine forc_src
#endif
  !
  subroutine bulk_forcing_src(forcing_type,is_forced,nx,ny,nz,dx,dy,dz,f_t12, &
                              nh_d,nh_u,dzc,dzf,rho,u,v,w,f)
    !
    ! to add bulk velocity forcing following three control strategies:
    !   --> constant flow rate (CFR);
    !   --> constant pressure gradient (CPR);
    !   --> constant power input (CPI) - not implemented yet.
    !
    use mod_debug, only: cmpt_mean
    use mod_param, only: lx,ly,lz,bvel_x,bvel_y,bvel_z,dpdl_x,dpdl_y,dpdl_z
#if !defined(_USE_VOF)
    use mod_param, only: rho_sp
#endif
    !
    implicit none
    !
    character(len=3), intent(in   )                                     :: forcing_type
    logical         , intent(in   ), dimension(3)                       :: is_forced
    integer         , intent(in   )                                     :: nx,ny,nz
    real(rp)        , intent(in   )                                     :: dx,dy,dz
    real(rp)        , intent(in   )                                     :: f_t12
    integer         , intent(in   )                                     :: nh_d,nh_u
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp)        , intent(in   ), dimension(     0:,     0:,     0:) :: rho
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp)        , intent(out  ), dimension(3)                       :: f
    !
    real(rp) :: meanvel,rhof
    integer  :: i,j,k
    !@cuf attributes(managed) :: dzc, dzf, u, v, w, rho
    !
#if !defined(_USE_VOF)
    rhof = rho_sp
#endif
    !
    f(1:3) = 0._rp
    select case(forcing_type)
    case('cfr')
      !
      ! ensure the flow rate is exactly constant
      !  note: to be generalized when the density field is not uniform
      !  (e.g., two-phase, low-mach flows)
      !
      if(is_forced(1)) then
        call cmpt_mean(nx,ny,nz,nh_d,nh_u,dx,dy,dzf,lx,ly,lz,u,meanvel)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              u(i,j,k) = u(i,j,k) + (bvel_x-meanvel)
            enddo
          enddo
        enddo
        !$acc end kernels
        f(1) = bvel_x-meanvel
      endif
      if(is_forced(2)) then
        call cmpt_mean(nx,ny,nz,nh_d,nh_u,dx,dy,dzf,lx,ly,lz,v,meanvel)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              v(i,j,k) = v(i,j,k) + (bvel_y-meanvel)
            enddo
          enddo
        enddo
        !$acc end kernels
        f(2) = bvel_y-meanvel
      endif
      if(is_forced(3)) then
        call cmpt_mean(nx,ny,nz,nh_d,nh_u,dx,dy,dzc,lx,ly,lz,w,meanvel)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              w(i,j,k) = w(i,j,k) + (bvel_z-meanvel)
            enddo
          enddo
        enddo
        !$acc end kernels
        f(3) = bvel_z-meanvel
      endif
    case('cpr')
      !
      ! ensure the pressure gradient is exactly constant
      !
      if(is_forced(1)) then
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
#if defined(_USE_VOF)
              rhof     = 0.5_rp*(rho(i+1,j,k)+rho(i,j,k))
#endif
              u(i,j,k) = u(i,j,k) + f_t12*dpdl_x/rhof
            enddo
          enddo
        enddo
        !$acc end kernels
      endif
      if(is_forced(2)) then
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
#if defined(_USE_VOF)
              rhof     = 0.5_rp*(rho(i,j+1,k)+rho(i,j,k))
#endif
              v(i,j,k) = v(i,j,k) + f_t12*dpdl_y/rhof
            enddo
          enddo
        enddo
        !$acc end kernels
      endif
      if(is_forced(3)) then
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
#if defined(_USE_VOF)
              rhof     = 0.5_rp*(rho(i,j,k+1)+rho(i,j,k))
#endif
              w(i,j,k) = w(i,j,k) + f_t12*dpdl_z/rhof
            enddo
          enddo
        enddo
        !$acc end kernels
      endif
    case('cpi')
      !
      ! ensure the power input is exactly constant
      !
      if(is_forced(1).or.is_forced(2).or.is_forced(3)) then
        call flutas_error('CPI control strategy not implemented yet')
      endif
    case default
      call flutas_error('Invalid selected control strategy - check dns.in')
    end select
    !
    return
  end subroutine bulk_forcing_src
  !
end module mod_source
