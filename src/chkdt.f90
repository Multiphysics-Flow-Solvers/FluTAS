!
! SPDX-License-Identifier: MIT
!
module mod_chkdt
  !
  use mpi
  use mod_param, only: rho1,rho2,mu1,mu2, &
#if defined(_HEAT_TRANSFER)
                       kappa1,kappa2,cp1,cp2, &
#endif
                       sigma,gacc
  use mod_common_mpi,  only: ierr
  use mod_types
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: chkdt
  !
  contains
  !
  subroutine chkdt(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
    !
    ! computes maximum allowed timestep
    !  as in Kang, Fedkiw and Liu, 
    !  "A boundary condition Capturing Method for Multiphase Incompressible Flow",
    !  JSC 2000
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzci,dzfi
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(out)                                     :: dtmax
    !
    real(rp) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(rp) :: dtix,dtiy,dtiz,dtic,dtiv,dtik,dtig,dti,dlmin,dlmini
#if defined(_HEAT_TRANSFER)
    real(rp) :: dtth
#endif
    !@cuf attributes(managed) :: u, v, w, dzci, dzfi
    integer :: i,j,k
    !
    dti = 0._rp
    !
#if defined(_OPENACC)
    !$acc parallel loop collapse(3) reduction(max:dti)
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,u,v,w,dxi,dyi,dzi,dzci,dzfi) &
    !$OMP PRIVATE(i,j,k,ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) &
    !$OMP REDUCTION(max:dti)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ux = abs(u(i,j,k))
          vx = 0.25_rp*abs( v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k) )
          wx = 0.25_rp*abs( w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1) )
          dtix = ux*dxi+vx*dyi+wx*dzfi(k)
          !
          uy = 0.25_rp*abs( u(i,j,k)+u(i,j+1,k)+u(i-1,j+1,k)+u(i-1,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25_rp*abs( w(i,j,k)+w(i,j+1,k)+w(i,j+1,k-1)+w(i,j,k-1) )
          dtiy = uy*dxi+vy*dyi+wy*dzfi(k)
          !
          uz = 0.25_rp*abs( u(i,j,k)+u(i-1,j,k)+u(i-1,j,k+1)+u(i,j,k+1) )
          vz = 0.25_rp*abs( v(i,j,k)+v(i,j-1,k)+v(i,j-1,k+1)+v(i,j,k+1) )
          wz = abs(w(i,j,k))
          dtiz = uz*dxi+vz*dyi+wz*dzci(k)
          !
          dti = max(dti,dtix,dtiy,dtiz)
          !
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end parallel loop 
#else
    !$OMP END PARALLEL DO
#endif
    call mpi_allreduce(MPI_IN_PLACE,dti,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(dti.eq.0._rp) dti = 1._rp
    dtic   = dti
    !
    dlmin  = min(1._rp/dxi,1._rp/dyi,1._rp/dzi)
    dlmin  = min(dlmin,minval(1._rp/dzfi(:))) ! minimum of dzf is an estimate on the safe side
    dlmini = dlmin**(-1)
    dtiv   = max(mu1/rho1,mu2/rho2)*(2._rp*(3._rp*dlmini**2))
    dtik   = sqrt(sigma/(min(rho1,rho2))*dlmini**3)
    dtig   = sqrt(maxval(abs(gacc(:)))*dlmini)
    dtmax  = 2._rp*(dtic+dtiv+sqrt((dtic+dtiv)**2+4._rp*(dtig**2+dtik**2)))**(-1)
#if defined(_HEAT_TRANSFER)
    dtth   = (max(kappa1/(rho1*cp1),kappa2/(rho2*cp2))*(2._rp*(3._rp*dlmini**2)))**(-1)
    dtmax  = min(dtmax,dtth)
#endif
    !
    return
  end subroutine chkdt
end module mod_chkdt
