!
! SPDX-License-Identifier: MIT
!
module mod_fillps
  !
  use mod_types, only: rp
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: fillps
  !
  contains
  !
  subroutine fillps(nx,ny,nz,nh_d,nh_u,dxi,dyi,dzi,dzfi,dti,rho0,u,v,w,p)
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzfi
    real(rp), intent(in )                                     :: dti,rho0
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(out), dimension(0:,0:,0:)                :: p
    !
    real(rp) :: dtidxi,dtidyi
    integer  :: i,j,k,im,jm,km
    !@cuf attributes(managed) :: p, u, v, w, dzfi
    !
    dtidxi = dti*dxi
    dtidyi = dti*dyi
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,rho0,p,u,v,w,dtidxi,dtidyi,dti,dzci) &
    !$OMP PRIVATE(i,j,k,im,jm,km)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          im = i-1
          jm = j-1
          km = k-1
          !
          p(i,j,k) = &
                     ((w(i,j,k)-w(i,j,km))*dti*dzfi(k) + &
                      (v(i,j,k)-v(i,jm,k))*dtidyi      + &
                      (u(i,j,k)-u(im,j,k))*dtidxi      )
          !
#if defined(_CONSTANT_COEFFS_POISSON)
          p(i,j,k) = p(i,j,k)*rho0
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
  end subroutine fillps
  !
end module mod_fillps
