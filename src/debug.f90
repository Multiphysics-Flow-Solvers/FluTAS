!
! SPDX-License-Identifier: MIT
!
module mod_debug
  !
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_types
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: cmpt_mean
  !
  contains
  !
  subroutine cmpt_mean(nx,ny,nz,nh_d,nh_p,dx,dy,dzf,lx,ly,lz,p,mean)
    !
    ! compute the mean value of an observable over the entire domain
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    integer , intent(in )                                     :: nh_d,nh_p
    real(rp), intent(in )                                     :: dx,dy
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzf
    real(rp), intent(in )                                     :: lx,ly,lz
    real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    real(rp), intent(out)                                     :: mean
    !
    !@cuf attributes(managed) :: p, dzf
    integer :: i,j,k
    !
    mean = 0._rp
    !
#if defined(_OPENACC)
    !$acc parallel loop collapse(3) reduction(+:mean)
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,dzf) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP REDUCTION(+:mean)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          mean = mean + p(i,j,k)*dx*dy*dzf(k)
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end parallel loop 
#else
    !$OMP END PARALLEL DO
#endif
    !
    call mpi_allreduce(MPI_IN_PLACE,mean,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    mean = mean/(1._rp*lx*ly*lz)
    !
    return
  end subroutine cmpt_mean
  !
end module mod_debug
