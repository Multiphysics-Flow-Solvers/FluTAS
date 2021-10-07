!
! SPDX-License-Identifier: MIT
!
module mod_rks
  !
  use mod_common_mpi
  use mod_moms    , only: momtad
  use mod_types
#if defined(_OPENACC)
  use cudafor
#endif
  !
  implicit none
  !
  private
  public  :: ab2_tmp
  !
  contains
  !
  subroutine ab2_tmp(is_first,nx,ny,nz,dxi,dyi,dzi,nh_u,nh_t,dt,dto,rho,cpp,kappa,u,v,w, &
                     tmp,dtmpdtrk,dtmpdtrkold)
    !
    ! second order Adams-Bashforth scheme
    ! for time integration of the temperature
    !
    implicit none
    !
    logical,  intent(inout)                                     :: is_first
    integer,  intent(in   )                                     :: nx,ny,nz
    real(rp), intent(in   )                                     :: dxi,dyi,dzi
    integer , intent(in   )                                     :: nh_u,nh_t
    real(rp), intent(in   )                                     :: dt,dto
    real(rp), intent(in   ), dimension(     0:,     0:,     0:) :: rho,cpp,kappa
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(inout), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: tmp 
    real(rp), intent(out  ), dimension(      :,      :,      :) :: dtmpdtrk
    real(rp), intent(inout), dimension(      :,      :,      :) :: dtmpdtrkold
    ! 
    real(rp) :: factor1,factor2
    integer  :: i,j,k,ii
    !
#if defined(_OPENACC)
    attributes(managed) :: rho,cpp,kappa
    attributes(managed) :: u,v,w
    attributes(managed) :: tmp,dtmpdtrk,dtmpdtrkold
    integer             :: istat
#endif
    !
    factor1 = dt*(1._rp+0.5_rp*(dt/dto))
    factor2 = dt*(     -0.5_rp*(dt/dto))
    if(is_first) then
     factor1 = dt*(1._rp)
     factor2 = dt*(0._rp)
     is_first = .false.
    endif
    !
    call momtad(nx,ny,nz,dxi,dyi,dzi,nh_u,nh_t, &
                kappa,cpp,rho,tmp,u,v,w,dtmpdtrk)
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(nx,ny,nz,factor1,factor2,tmp,dtmpdtrk,dtmpdtrkold)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          tmp(i,j,k) = tmp(i,j,k) + factor1*dtmpdtrk(i,j,k) + factor2*dtmpdtrkold(i,j,k)
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
    !@cuf istat=cudaDeviceSynchronize()
#else
    !$OMP END PARALLEL DO
#endif
    !
    return
  end subroutine ab2_tmp
  !
end module mod_rks
