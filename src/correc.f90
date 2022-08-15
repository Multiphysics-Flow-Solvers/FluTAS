!
! SPDX-License-Identifier: MIT
!
module mod_correc
  !
  use mod_types     , only: rp
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: correc
  !
  contains
  !
  subroutine correc(nx,ny,nz,nh_d,nh_u,dxi,dyi,dzi,dzci,dt,rho0,p,u,v,w,rho)
    !
    ! corrects the velocity so that the prescribed divergence is imposed
    !
    implicit none
    !
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   )                                     :: nh_d,nh_u
    real(rp), intent(in   )                                     :: dxi,dyi,dzi
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzci
    real(rp), intent(in   )                                     :: dt,rho0
    real(rp), intent(in   ), dimension(0:,0:,0:)                :: p
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in   ), dimension(0:,0:,0:)                :: rho
    !
    real(rp) :: factori,factorj
    real(rp) :: rhox,rhoy,rhoz
    real(rp) :: rho0i
    integer  :: i,j,k,ip,jp,kp
    !@cuf attributes(managed) :: u, v, w, rho, p, dzci
    !
    rho0i = 1._rp/rho0
    !
    factori = dt*dxi
    factorj = dt*dyi
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,factori,factorj,dzci,dt,p,u,v,w,rho) &
    !$OMP PRIVATE(i,j,k,ip,jp,kp,rho0i,rhox,rhoy,rhoz)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1
          jp = j+1
          kp = k+1
          !
#if defined(_CONSTANT_COEFFS_POISSON)
          u(i,j,k) = u(i,j,k) - factori*(    p( ip,j,k)-p( i,j,k) )*rho0i
          v(i,j,k) = v(i,j,k) - factorj*(    p( i,jp,k)-p( i,j,k) )*rho0i
          w(i,j,k) = w(i,j,k) - dt*dzci(k)*( p( i,j,kp)-p( i,j,k) )*rho0i
#else
          rhox = 0.5_rp*(rho(ip,j,k)+rho(i,j,k))
          rhoy = 0.5_rp*(rho(i,jp,k)+rho(i,j,k))
          rhoz = 0.5_rp*(rho(i,j,kp)+rho(i,j,k))
          !
          u(i,j,k) = u(i,j,k) - factori*(    p( ip,j,k)-p( i,j,k) )/rhox
          v(i,j,k) = v(i,j,k) - factorj*(    p( i,jp,k)-p( i,j,k) )/rhoy
          w(i,j,k) = w(i,j,k) - dt*dzci(k)*( p( i,j,kp)-p( i,j,k) )/rhoz
#endif
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
  end subroutine correc
  !
end module mod_correc
