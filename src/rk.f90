!
! SPDX-License-Identifier: MIT
!
module mod_rk
  !
  use mpi
#if defined(_USE_VOF)
  use mod_mom   , only: momad_xyz_tw_cen,momad_xyz_tw_fll
#else
  use mod_mom   , only: momad_xyz_sp
#endif
  use mod_sanity, only: flutas_error 
  use mod_types
  !@cuf use cudafor
  !
  implicit none
  !
  real(rp), allocatable, dimension(:,:,:) :: dudtrk, dvdtrk, dwdtrk
  !@cuf attributes(managed) :: dudtrk, dvdtrk, dwdtrk
  !
  private
  public  :: rk,cmpt_time_factors
  !
  contains
  !
  subroutine rk(space_scheme_mom,f_t1,f_t2,nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w, &
                mu,rho,dudtrko,dvdtrko,dwdtrko)
    !
    ! subroutine to compute the velocity prediction
    !
    ! Note: --> suitable for both 2nd order Adams-Bashforth and 3rd order low-storage Runge-Kutta;
    !       --> source terms are included in the main program;
    !
    implicit none
    !
    character(len=3), intent(in   )                                     :: space_scheme_mom
    real(rp)        , intent(in   )                                     :: f_t1,f_t2
    integer         , intent(in   )                                     :: nx,ny,nz
    real(rp)        , intent(in   )                                     :: dxi,dyi,dzi
    integer         , intent(in   )                                     :: nh_d,nh_u
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzci,dzfi
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp)        , intent(in   ), dimension(     0:,     0:,     0:) :: mu,rho
    real(rp)        , intent(inout), dimension(      :,      :,      :) :: dudtrko,dvdtrko,dwdtrko
    !
    integer :: i,j,k
    !
    !@cuf attributes(managed) :: dzci, dzfi, u, v, w, dudtrko, dvdtrko, dwdtrko
    !@cuf attributes(managed) :: mu, rho
    !
    ! Allocate once
    if( (.not.allocated(dudtrk)) .or. (.not.allocated(dvdtrk)) .or. (.not.allocated(dwdtrk) ) ) then
      allocate (dudtrk(nx,ny,nz),dvdtrk(nx,ny,nz),dwdtrk(nx,ny,nz))
    endif
    !
    ! compute the r.h.s. of the momentum equation 
    !
#if defined(_USE_VOF)
    select case(space_scheme_mom)
    case('cen') ! ii order central scheme in divergence form
      call momad_xyz_tw_cen(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,mu,rho,dudtrk,dvdtrk,dwdtrk)  
    case('fll') ! ii order flux limiter in gradient form
      call momad_xyz_tw_fll(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,mu,rho,dudtrk,dvdtrk,dwdtrk)  
    case default
      call flutas_error('Invalid space discretization for momentum equation - check dns.in')
    end select
#else
    call momad_xyz_sp(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dudtrk,dvdtrk,dwdtrk)  
#endif
    !
    ! compute the velocity prediction 
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(nx,ny,nz,f_t1,f_t2,u,v,w,dudtrk,dvdtrk,dwdtrk,dudtrko,dvdtrko,dwdtrko)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          u(i,j,k) = u(i,j,k) + f_t1*dudtrk(i,j,k) + f_t2*dudtrko(i,j,k)
          v(i,j,k) = v(i,j,k) + f_t1*dvdtrk(i,j,k) + f_t2*dvdtrko(i,j,k)
          w(i,j,k) = w(i,j,k) + f_t1*dwdtrk(i,j,k) + f_t2*dwdtrko(i,j,k)
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
          !
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end kernels 
#else
    !$OMP END PARALLEL DO
#endif
    !
    return
  end subroutine rk
  !
  subroutine cmpt_time_factors(time_scheme,restart,istep,irk_ss,rkcoeff,dt,dto,f_t1,f_t2,f_t12)
    !
    ! compute f_t1 and f_t2 according to the input time_scheme 
    !
    implicit none
    !
    character(len=3), intent(in )                 :: time_scheme
    logical         , intent(in )                 :: restart
    integer         , intent(in )                 :: istep,irk_ss
    real(rp)        , intent(in ), dimension(2,3) :: rkcoeff
    real(rp)        , intent(in )                 :: dt,dto
    real(rp)        , intent(out)                 :: f_t1,f_t2,f_t12
    !
    select case(time_scheme)
    case('ab2')
      if(istep.eq.1) then ! we use 1st order Euler in the first time-step 
        f_t1 = 1._rp*dt
        f_t2 = 0._rp*dt
      else
        f_t1 = (1._rp+0.5_rp*(dt/dto))*dt
        f_t2 = (     -0.5_rp*(dt/dto))*dt
      endif
    case('rk3')
      f_t1 = rkcoeff(1,irk_ss)*dt
      f_t2 = rkcoeff(2,irk_ss)*dt
    case default
      call flutas_error('Invalid time discretization: only AB2 and RK3 available for now - check dns.in')
    end select
    f_t12 = f_t1+f_t2
    !
    return
  end subroutine cmpt_time_factors
  !
end module mod_rk
