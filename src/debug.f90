!
! SPDX-License-Identifier: MIT
!
module mod_debug
  !
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_types
#if defined(_OPENACC)
  use cudafor
#endif
  !
  implicit none
  !
  private
  public  :: chkmean,chk_helmholtz
  !
  contains
  !
  subroutine chkmean(nx,ny,nz,dims,nh_d,nh_p,dzlzi,p,mean)
    !
    ! compute the mean value of an observable over the entire domain
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    integer , intent(in ), dimension(3)                       :: dims
    integer , intent(in )                                     :: nh_d,nh_p
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzlzi
    real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    real(rp), intent(out)                                     :: mean
    !
#if defined(_OPENACC)
    attributes(managed) :: p,dzlzi
    integer :: istat
#endif
    integer :: i,j,k
    !
    mean = 0._rp
    !
#if defined(_OPENACC)
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,dzlzi) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP REDUCTION(+:mean)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          mean = mean + p(i,j,k)*dzlzi(k)
        enddo
      enddo
    enddo
#if !defined(_OPENACC)
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    call mpi_allreduce(MPI_IN_PLACE,mean,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    mean = mean/(1._rp*nx*dims(1)*ny*dims(2))
    !
    return
  end subroutine chkmean
  !
  subroutine chk_helmholtz(nx,ny,nz,nh_d,dli,dzci,dzfi,alpha,fp,fpp,bc,c_or_f,diffmax)
    !
    ! this subroutine checks if the implementation of implicit diffusion is
    ! correct
    !
    implicit none
    !
    integer         , intent(in )                      :: nx,ny,nz
    integer         , intent(in )                      :: nh_d
    real(rp)        , intent(in ), dimension(2)        :: dli
    real(rp)        , intent(in ), dimension(1-nh_d:)  :: dzfi,dzci
    real(rp)        , intent(in )                      :: alpha
    real(rp)        , intent(in ), dimension(0:,0:,0:) :: fp,fpp
    character(len=1), intent(in ), dimension(0:1,3)    :: bc
    character(len=1), intent(in ), dimension(3)        :: c_or_f
    real(rp)        , intent(out)                      :: diffmax
    !
    integer, dimension(3) :: q
    real(rp) :: val
    integer  :: i,j,k,im,ip,jm,jp,km,kp
    integer  :: idir
    !
    q(:) = 0
    do idir = 1,3
      if(bc(1,idir).ne.'P'.and.c_or_f(idir).eq.'f') q(idir) = 1
    enddo
    select case(c_or_f(3))
    !
    ! need to compute the maximum difference!
    !
    case('c')
      diffmax = 0._rp
      do k=1,nz-q(3)
        kp = k + 1
        km = k - 1
        do j=1,ny-q(2)
          jp = j + 1
          jm = j - 1
          do i=1,nx-q(1)
            ip = i + 1
            im = i - 1
            val =  fpp(i,j,k)+(1._rp/alpha)*( &
                  (fpp(ip,j,k)-2._rp*fpp(i,j,k)+fpp(im,j,k))*(dli(1)**2) + &
                  (fpp(i,jp,k)-2._rp*fpp(i,j,k)+fpp(i,jm,k))*(dli(2)**2) + &
                 ((fpp(i,j,kp)-fpp(i,j,k ))*dzci(k ) - &
                  (fpp(i,j,k )-fpp(i,j,km))*dzci(km))*dzfi(k) )
            val = val*alpha
            diffmax = max(diffmax,abs(val-fp(i,j,k)))
            !ii = coord(1)*nx+i
            !jj = coord(2)*ny+j
            !if(abs(val-fp(i,j,k)).gt.1.e-8) print*, 'Large difference : ', val-fp(i,j,k),ii,jj,k
          enddo
        enddo
      enddo
    case('f')
      diffmax = 0._rp
      do k=1,nz-q(3)
        kp = k + 1
        km = k - 1
        do j=1,ny-q(2)
          jp = j + 1
          jm = j - 1
          do i=1,nx-q(1)
            ip = i + 1
            im = i - 1
            val =  fpp(i,j,k)+(1._rp/alpha)*( &
                  (fpp(ip,j,k)-2._rp*fpp(i,j,k)+fpp(im,j,k))*(dli(1)**2) + &
                  (fpp(i,jp,k)-2._rp*fpp(i,j,k)+fpp(i,jm,k))*(dli(2)**2) + &
                 ((fpp(i,j,kp)-fpp(i,j,k ))*dzfi(kp) - &
                  (fpp(i,j,k )-fpp(i,j,km))*dzfi(k ))*dzci(k) )
            val = val*alpha
            diffmax = max(diffmax,abs(val-fp(i,j,k)))
            !ii = coord(1)*nx+i
            !jj = coord(2)*ny+j
            !if(abs(val-fp(i,j,k)).gt.1.e-8) print*, 'Large difference : ', val,fp(i,j,k),ii,jj,k
          enddo
        enddo
      enddo
    end select
    call mpi_allreduce(MPI_IN_PLACE,diffmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    !
    return
  end subroutine chk_helmholtz
  !
end module mod_debug
