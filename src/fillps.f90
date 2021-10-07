!
! SPDX-License-Identifier: MIT
!
module mod_fillps
  !
  use mod_types, only: rp
#if defined(_OPENACC)
  use cudafor
#endif
  !
  implicit none
  !
  private
  public  :: fillps
  !
  contains
  !
  subroutine fillps(nx,ny,nz,nh_d,nh_u,dxi,dyi,dzi,dzci,dzfi,dti,rho0,up,vp,wp,p)
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzci,dzfi
    real(rp), intent(in )                                     :: dti,rho0
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: up,vp,wp
    real(rp), intent(out), dimension(0:,0:,0:)                :: p
    !
    real(rp), dimension(1-nh_d:nz+nh_d) :: dtidzfi
    real(rp) :: dtidxi,dtidyi
    integer  :: i,j,k,im,jm,km
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: p, up,vp,wp,dtidzfi,dzci,dzfi
#endif
    !
    dtidxi = dti*dxi
    dtidyi = dti*dyi
    !$acc kernels
    do k=1-nh_d,nz+nh_d
      dtidzfi(k) = dti*dzfi(k)
    enddo
    !$acc end kernels
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,up,vp,wp,dtidzfi,dtidyi,dtidxi) &
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
                     ((wp(i,j,k)-wp(i,j,km))*dtidzfi(k)+ &
                      (vp(i,j,k)-vp(i,jm,k))*dtidyi    + &
                      (up(i,j,k)-up(im,j,k))*dtidxi    )
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
    !@cuf istat=cudaDeviceSynchronize()
#else
    !$OMP END PARALLEL DO
#endif
    !
    return
  end subroutine fillps
  !
end module mod_fillps
