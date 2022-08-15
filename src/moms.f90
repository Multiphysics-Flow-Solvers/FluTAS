!
! SPDX-License-Identifier: MIT
!
module mod_moms
  !
  use mod_common_mpi, only: comm_cart,ierr
  use mod_gradls    , only: weno5,weno5_old
  use mod_types
  use mpi
  !@cuf use cudafor
  !
  implicit none
  !
  private
#if defined(_USE_VOF)
  public  :: momtad_tw
#else
  public  :: momtad_sp
#endif
  !
  contains
  !
#if defined(_USE_VOF)
  subroutine momtad_tw(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,nh_t,kappa,cp,rho, &
                       dzci,dzfi,s,u,v,w,dsdt)
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_d,nh_u,nh_t
    real(rp), intent(in ), dimension(     0:,     0:,     0:) :: kappa,cp,rho
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzci,dzfi
    real(rp), intent(in ), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: s 
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(out), dimension(      :,      :,      :) :: dsdt
    !
    integer  :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(rp) :: kappaxp,kappaxm,kappayp,kappaym,kappazp,kappazm
    real(rp) :: rhocpci
    !@cuf attributes(managed) :: u, v, w, kappa, cp, rho, s, dsdt, dzci, dzfi
    !
    call weno5(nx,ny,nz,dxi,dyi,dzi,nh_u,.true.,s,u,v,w,dsdt)
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(kappaxp,kappaxm,kappayp,kappaym,kappazp,kappazm) &
    !$OMP PRIVATE(dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm) &
    !$OMP PRIVATE(rhocpci) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,rho,kappa,cp,s,dsdt)
    !$OMP SHARED(dzci,dzfi)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i + 1
          im = i - 1
          jp = j + 1
          jm = j - 1
          kp = k + 1
          km = k - 1
          !
          kappaxp = 0.5_rp*(kappa(ip,j,k)+kappa(i,j,k))
          kappaxm = 0.5_rp*(kappa(im,j,k)+kappa(i,j,k))
          kappayp = 0.5_rp*(kappa(i,jp,k)+kappa(i,j,k))
          kappaym = 0.5_rp*(kappa(i,jm,k)+kappa(i,j,k))
          kappazp = 0.5_rp*(kappa(i,j,kp)+kappa(i,j,k))
          kappazm = 0.5_rp*(kappa(i,j,km)+kappa(i,j,k))
          !
          rhocpci = 1._rp/(rho(i,j,k)*cp(i,j,k))
          !
          dsdxp = (s(ip,j,k)-s(i ,j,k))*dxi
          dsdxm = (s(i ,j,k)-s(im,j,k))*dxi
          dsdyp = (s(i,jp,k)-s(i,j ,k))*dyi
          dsdym = (s(i,j ,k)-s(i,jm,k))*dyi
          dsdzp = (s(i,j,kp)-s(i,j,k ))*dzci(k )
          dsdzm = (s(i,j,k )-s(i,j,km))*dzci(km)
          !
          dsdt(i,j,k) = + dsdt(i,j,k) &
                        + rhocpci*( &
                                    (kappaxp*dsdxp-kappaxm*dsdxm)*dxi +   &
                                    (kappayp*dsdyp-kappaym*dsdym)*dyi +   &
                                    (kappazp*dsdzp-kappazm*dsdzm)*dzfi(k) &
                                  ) 
          !
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end kernels 
#else
    !$OMP END PARALLEL DO
#endif
    return
    !
  end subroutine momtad_tw
#else
  subroutine momtad_sp(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,nh_t, &
                       dzci,dzfi,s,u,v,w,dsdt)
    !
    use mod_param, only: kappa_sp,rho_sp,cp_sp
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_d,nh_u,nh_t
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzci,dzfi
    real(rp), intent(in ), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: s 
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(out), dimension(      :,      :,      :) :: dsdt
    !
    integer  :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(rp) :: alpha_sp
    !@cuf attributes(managed) :: u, v, w, s, dsdt, dzci, dzfi
    !
    call weno5(nx,ny,nz,dxi,dyi,dzi,nh_u,.true.,s,u,v,w,dsdt)
    !
    alpha_sp = kappa_sp/(rho_sp*cp_sp)
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm) &
    !$OMP PRIVATE(alpha_sp) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,s,dsdt)
    !$OMP SHARED(dzci,dzfi)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i + 1
          im = i - 1
          jp = j + 1
          jm = j - 1
          kp = k + 1
          km = k - 1
          !
          dsdxp = (s(ip,j,k)-s(i ,j,k))*dxi
          dsdxm = (s(i ,j,k)-s(im,j,k))*dxi
          dsdyp = (s(i,jp,k)-s(i,j ,k))*dyi
          dsdym = (s(i,j ,k)-s(i,jm,k))*dyi
          dsdzp = (s(i,j,kp)-s(i,j,k ))*dzci(k )
          dsdzm = (s(i,j,k )-s(i,j,km))*dzci(km)
          !
          dsdt(i,j,k) = + dsdt(i,j,k) &
                        + alpha_sp*( &
                                     (dsdxp-dsdxm)*dxi +   &
                                     (dsdyp-dsdym)*dyi +   &
                                     (dsdzp-dsdzm)*dzfi(k) &
                                   ) 
          !
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end kernels 
#else
    !$OMP END PARALLEL DO
#endif
    return
    !
  end subroutine momtad_sp
#endif
  !
end module mod_moms
