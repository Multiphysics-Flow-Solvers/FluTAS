!
! SPDX-License-Identifier: MIT
!
module mod_cmpt_divth
  !
  ! note: these subroutines will be used for future extension to low Mach flows
  !
  use mod_param, only: gam_g
  use mod_types
  use mod_common_mpi
  !@cuf use cudafor
  !
  implicit none
  !
  public  :: cmpt_divth,cmpt_density,cmpt_e_div
  private 
  !
  contains
  !
  subroutine cmpt_divth(nx,ny,nz,dxi,dyi,dzi,nh_t,p0,dp0dt,vof,kappa,s,div_th)
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_t
    real(rp), intent(in )                                     :: p0,dp0dt
    real(rp), intent(in ), dimension(     0:,     0:,     0:) :: vof,kappa
    real(rp), intent(in ), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: s 
    real(rp), intent(out), dimension(     0:,     0:,     0:) :: div_th
    !
    real(rp) :: dsdx,dsdy,dsdz,diff, &
                dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(rp) :: kappaxp,kappaxm,kappayp,kappaym, &
                kappazp,kappazm
    integer  :: i,j,k,im,jm,km,ip,jp,kp
    !
    !@cuf attributes(managed) :: vof, kappa, div_th, s
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(im,ip,jm,jp,km,kp) &
    !$OMP PRIVATE(dsdxm,dsdxp,dsdym,dsdyp,dsdzm,dsdzp) &
    !$OMP PRIVATE(kappaxm,kappaxp,kappaym,kappayp,kappazm,kappazp) &
    !$OMP PRIVATE(diff) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,gam_g,dp0dt,p0,vof,kappa,s)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1
          !
          dsdxm = (s(i ,j,k)-s(im,j,k))*dxi
          dsdxp = (s(ip,j,k)-s(i ,j,k))*dxi
          dsdym = (s(i,j ,k)-s(i,jm,k))*dyi
          dsdyp = (s(i,jp,k)-s(i,j ,k))*dyi
          dsdzm = (s(i,j,k )-s(i,j,km))*dzi
          dsdzp = (s(i,j,kp)-s(i,j ,k))*dzi
          !
          kappaxm = 0.5_rp*(kappa(i,j,k)+kappa(im,j,k))
          kappaxp = 0.5_rp*(kappa(i,j,k)+kappa(ip,j,k))
          kappaym = 0.5_rp*(kappa(i,j,k)+kappa(i,jm,k))
          kappayp = 0.5_rp*(kappa(i,j,k)+kappa(i,jp,k))
          kappazm = 0.5_rp*(kappa(i,j,k)+kappa(i,j,km))
          kappazp = 0.5_rp*(kappa(i,j,k)+kappa(i,j,kp))
          ! 
          diff = (dsdxp*kappaxp-dsdxm*kappaxm)*dxi + &
                 (dsdyp*kappayp-dsdym*kappaym)*dyi + &
                 (dsdzp*kappazp-dsdzm*kappazm)*dzi
          ! 
          div_th(i,j,k) = (((gam_g-1._rp)/gam_g*diff - &
                                   (dp0dt)/gam_g)/p0)*vof(i,j,k)
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
  end subroutine cmpt_divth
  !
  subroutine cmpt_density(nx,ny,nz,dx,dy,dz,nh_t,p0,ri,rho2,psi,tmp,gas_mass,flag,rho0,rho)
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dx,dy,dz
    integer , intent(in )                                     :: nh_t
    real(rp), intent(in )                                     :: p0,ri,rho2
    real(rp), intent(in ), dimension(     0:,     0:,     0:) :: psi
    real(rp), intent(in ), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: tmp 
    real(rp), intent(out)                                     :: gas_mass
    logical , intent(in)                                      :: flag
    real(rp), intent(out)                                     :: rho0
    real(rp), intent(out), dimension(     0:,     0:,     0:) :: rho
    !
    integer :: i,j,k
    real(rp):: rho_gas,rho_liq
    !@cuf attributes(managed) :: psi, tmp, rho
    !
    if(flag) gas_mass = 0._rp
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(rho_gas,rho_liq) &
    !$OMP SHARED(flag,nx,ny,nz,dx,dy,dz,p0,ri,tmp,rho2,psi,rho)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          rho_gas    = p0/(ri*tmp(i,j,k))
          rho_liq    = rho2
          rho(i,j,k) = psi(i,j,k)*rho_gas+(1._rp-psi(i,j,k))*rho_liq
          if(flag) gas_mass   = gas_mass + psi(i,j,k)*rho_gas*dx*dy*dz
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
    rho0 = minval(rho(1:nx,1:ny,1:nz))
    call mpi_allreduce(MPI_IN_PLACE,rho0,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
    if(flag) call mpi_allreduce(MPI_IN_PLACE,gas_mass,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    return
  end subroutine cmpt_density
  !
  subroutine cmpt_e_div(nx,ny,nz,dx,dy,dz,lx,ly,lz,nh_u,div_th,u,v,w,e_div,e_div_mean)
    !
    implicit none
    !
    integer,  intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dx,dy,dz
    real(rp), intent(in )                                     :: lx,ly,lz
    integer , intent(in )                                     :: nh_u
    real(rp), intent(in ), dimension(      0:,    0:,     0:) :: div_th
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(out)                                     :: e_div,e_div_mean
    !
    integer :: i,j,k
    real(rp):: lxi,lyi,lzi,dxi,dyi,dzi
    real(rp):: div_s,e_div_v
    !
    !@cuf attributes(managed) :: div_th, u, v, w
    !
    lxi = 1._rp/lx
    lyi = 1._rp/ly
    lzi = 1._rp/lz
    dxi = 1._rp/dx
    dyi = 1._rp/dy
    dzi = 1._rp/dz
    !
    e_div = 0._rp
    e_div_mean = 0._rp
    !
#if defined(_OPENACC)
    ! [TODO] Add collapse(3) here or switch to kernels?
    !$acc parallel loop reduction(+:e_div_mean)
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(div_s,e_div_v) &
    !$OMP SHARED(nx,ny,nz,lx,ly,lz,dxi,dyi,xzi,div_th,e_div,e_div_mean)&
    !$OMP REDUCTION(+:e_div_mean)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          div_s = (u(i,j,k)-u(i-1,j,k))*dxi + &
                  (v(i,j,k)-v(i,j-1,k))*dyi + &
                  (w(i,j,k)-w(i,j,k-1))*dzi
          !
          e_div_v    = abs(div_th(i,j,k)-div_s)
          e_div      = max(e_div_v,e_div)
          e_div_mean = e_div_mean + e_div_v*dx*dy*dz/lx/ly/lz
          !
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end parallel loop
#else
    !$OMP END PARALLEL DO
#endif
    call mpi_allreduce(MPI_IN_PLACE,e_div     ,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,e_div_mean,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    return
  end subroutine cmpt_e_div
  !
end module mod_cmpt_divth
