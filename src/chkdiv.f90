!
! SPDX-License-Identifier: MIT
!
module mod_chkdiv
  !
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_types
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: chkdiv
  !
  contains
  !
  subroutine chkdiv(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzfi,u,v,w,divtot,divmax)
    !
    ! checks the divergence of the velocity field
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzfi
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(out)                                     :: divtot,divmax
    !
    real(rp) :: div
    !@cuf attributes(managed) :: u,v,w,dzfi
    integer :: i,j,k
    !
    divtot = 0._rp
    divmax = 0._rp
#if defined(_OPENACC)
    !$acc parallel loop collapse(3) reduction(max:divmax) reduction(+:divtot)
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,u,v,w,dxi,dyi,dzfi) &
    !$OMP PRIVATE(i,j,k,div) &
    !$OMP REDUCTION(+:divtot) &
    !$OMP REDUCTION(max:divmax)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          div = (w(i,j,k)-w(i,j,k-1))*dzfi(k) + &
                (v(i,j,k)-v(i,j-1,k))*dyi     + &
                (u(i,j,k)-u(i-1,j,k))*dxi
          divmax = max(divmax,abs(div))
          divtot = divtot + div
          !
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end parallel loop 
#else
    !$OMP END PARALLEL DO
#endif
    call mpi_allreduce(MPI_IN_PLACE,divtot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,divmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(myid.eq.0) print*, 'Total divergence = ', divtot, '| Maximum divergence = ', divmax
    !
    return
  end subroutine chkdiv
  !
end module mod_chkdiv
