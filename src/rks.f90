!
! SPDX-License-Identifier: MIT
!
module mod_rks
  !
  use mod_common_mpi
#if defined(_USE_VOF)
  use mod_moms, only: momtad_tw
#else
  use mod_moms, only: momtad_sp
#endif
  use mod_types
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: rk_sca
  !
  contains
  !
  subroutine rk_sca(f_t1,f_t2,nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,nh_t,rho,cpp,kappa,dzci,dzfi, &
                     u,v,w,tmp,dtmpdtrk,dtmpdtrkold)
    !
    ! subroutine to compute the temperature (or generic scalar)
    !
    ! Note: --> suitable for both 2nd order Adams-Bashforth and 3rd low-storage Runge-Kutta;
    !       --> source terms are included in the main program;
    !
    implicit none
    !
    real(rp), intent(in   )                                     :: f_t1,f_t2
    integer,  intent(in   )                                     :: nx,ny,nz
    real(rp), intent(in   )                                     :: dxi,dyi,dzi
    integer , intent(in   )                                     :: nh_d,nh_u,nh_t
    real(rp), intent(in   ), dimension(     0:,     0:,     0:) :: rho,cpp,kappa
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzci,dzfi
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(inout), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: tmp 
    real(rp), intent(out  ), dimension(      :,      :,      :) :: dtmpdtrk
    real(rp), intent(inout), dimension(      :,      :,      :) :: dtmpdtrkold
    ! 
    integer :: i,j,k
    !
    !@cuf attributes(managed) :: rho, cpp, kappa, u, v, w, tmp, dtmpdtrk, dtmpdtrkold
    !@cuf attributes(managed) :: dzci, dzfi
    !
#if defined(_USE_VOF)
    call momtad_tw(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,nh_t,kappa,cpp,rho, &
                   dzci,dzfi,tmp,u,v,w,dtmpdtrk)
#else
    call momtad_sp(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,nh_t, &
                   dzci,dzfi,tmp,u,v,w,dtmpdtrk)
#endif
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(nx,ny,nz,f_t1,f_t2,tmp,dtmpdtrk,dtmpdtrkold)
    !$OMP SHARED(dzci,dzfi)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          tmp(i,j,k) = tmp(i,j,k) + f_t1*dtmpdtrk(i,j,k) + f_t2*dtmpdtrkold(i,j,k)
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
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,dtmpdtrk,dtmpdtrkold)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          dtmpdtrkold(i,j,k) = dtmpdtrk(i,j,k)
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
  end subroutine rk_sca
  !
end module mod_rks
