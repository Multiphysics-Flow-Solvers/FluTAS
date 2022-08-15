!
! SPDX-License-Identifier: MIT
!
module mod_solver_cpu
  !
  use iso_c_binding , only: C_PTR
  use decomp_2d
  use mod_fft       , only: fft
  use mod_common_mpi, only: n_x,n_y,n_z
  use mod_types     , only: rp
  use profiler
  !
  implicit none
  !
  private
  public  :: solver_cpu
  !
  contains
  !
  subroutine solver_cpu(n,arrplan,normfft,lambdaxy,a,b,c,bcz,c_or_f,p)
    !
    implicit none
    !
    integer         , intent(in   ), dimension(3)             :: n
    type(C_PTR)     , intent(in   ), dimension(2,2)           :: arrplan
    real(rp)        , intent(in   )                           :: normfft
    real(rp)        , intent(in   ), dimension(n_z(1),n_z(2)) :: lambdaxy
    real(rp)        , intent(in   ), dimension(n_z(3))        :: a,b,c
    character(len=1), intent(in   ), dimension(0:1)           :: bcz
    character(len=1), intent(in   ), dimension(3)             :: c_or_f
    real(rp)        , intent(inout), dimension(0:,0:,0:)      :: p
    !
    real(rp), dimension(n_x(1),n_x(2),n_x(3)) :: px
    real(rp), dimension(n_y(1),n_y(2),n_y(3)) :: py
    real(rp), dimension(n_z(1),n_z(2),n_z(3)) :: pz
    integer :: q
    !
#if defined(_DECOMP_X)
    !$OMP WORKSHARE
    px(:,:,:) = p(1:n(1),1:n(2),1:n(3))
    !$OMP END WORKSHARE
#elif _DECOMP_Y
    !$OMP WORKSHARE
    py(:,:,:) = p(1:n(1),1:n(2),1:n(3))
    !$OMP END WORKSHARE
    call profiler_start("yc_to_x", tag = .true., tag_color = COLOR_PURPLE)
    call transpose_y_to_x(py,px)
    call profiler_stop("yc_to_x")
#else /*_DECOMP_Z*/
    !$OMP WORKSHARE
    pz(:,:,:) = p(1:n(1),1:n(2),1:n(3))
    !$OMP END WORKSHARE
    call profiler_start("zc_to_x", tag = .true., tag_color = COLOR_PURPLE)
    call transpose_z_to_x(pz,px)
    call profiler_stop("zc_to_x")
    !call transpose_z_to_y(pz,py)
    !call transpose_y_to_x(py,px)
#endif
    call fft(arrplan(1,1),px) ! fwd transform in x
    !
    call profiler_start("xc_to_y", tag = .true., tag_color = COLOR_PURPLE)
    call transpose_x_to_y(px,py)
    call profiler_stop("xc_to_y")
    !
    call fft(arrplan(1,2),py) ! fwd transform in y
    !
    call profiler_start("yc_to_z", tag = .true., tag_color = COLOR_PURPLE)
    call transpose_y_to_z(py,pz)
    call profiler_stop("yc_to_z")
    !
    q = 0
    if(c_or_f(3).eq.'f'.and.bcz(1).eq.'D') q = 1
    if(bcz(0)//bcz(1).eq.'PP') then
      call gaussel_periodic(n_z(1),n_z(2),n_z(3)-q,a,b,c,lambdaxy,pz)
    else
      call gaussel(         n_z(1),n_z(2),n_z(3)-q,a,b,c,lambdaxy,pz)
    endif
    !
    call profiler_start("zc_to_y", tag = .true., tag_color = COLOR_PURPLE)
    call transpose_z_to_y(pz,py)
    call profiler_stop("zc_to_y")
    !
    call fft(arrplan(2,2),py) ! bwd transform in y
    !
    call profiler_start("yc_to_x", tag = .true., tag_color = COLOR_PURPLE)
    call transpose_y_to_x(py,px)
    call profiler_stop("yc_to_x")
    !
    call fft(arrplan(2,1),px) ! bwd transform in x
    !
#if defined(_DECOMP_X)
    !$OMP WORKSHARE
    p(1:n(1),1:n(2),1:n(3)) = px(:,:,:)*normfft
    !$OMP END WORKSHARE
#elif _DECOMP_Y
    call profiler_start("xc_to_y", tag = .true., tag_color = COLOR_PURPLE)
    call transpose_x_to_y(px,py)
    call profiler_stop("xc_to_y")
    !$OMP WORKSHARE
    p(1:n(1),1:n(2),1:n(3)) = py(:,:,:)*normfft
    !$OMP END WORKSHARE
!#elif _DECOMP_Z
#else
    call profiler_start("xc_to_z", tag = .true., tag_color = COLOR_PURPLE)
    call transpose_x_to_z(px,pz)
    call profiler_stop("xc_to_z")
    !call transpose_x_to_y(px,py)
    !call transpose_y_to_z(py,pz)
    !$OMP WORKSHARE
    p(1:n(1),1:n(2),1:n(3)) = pz(:,:,:)*normfft
    !$OMP END WORKSHARE
#endif
    !
    return
  end subroutine solver_cpu
  !
  subroutine gaussel(nx,ny,n,a,b,c,lambdaxy,p)
    !
    implicit none
    !
    integer , intent(in   )                   :: nx,ny,n
    real(rp), intent(in   ), dimension(:)     :: a,b,c
    real(rp), intent(in   ), dimension(nx,ny) :: lambdaxy
    real(rp), intent(inout), dimension(:,:,:) :: p
    !
    real(rp), dimension(n) :: bb
    integer :: i,j
    !
    ! solve tridiagonal system
    !
    !$OMP PARALLEL DEFAULT(none) &
    !$OMP PRIVATE(i,j,bb) &
    !$OMP SHARED(nx,ny,n,a,b,c,lambdaxy,p)
    !$OMP DO COLLAPSE(2)
    do j=1,ny
      do i=1,nx
        bb(:) = b(1:n) + lambdaxy(i,j)
        call dgtsv_homebrewed(n,a,bb,c,p(i,j,1:n))
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    !
    return
  end subroutine gaussel
  !
  subroutine gaussel_periodic(nx,ny,n,a,b,c,lambdaxy,p)
    !
    implicit none
    !
    integer , intent(in   )                   :: nx,ny,n
    real(rp), intent(in   ), dimension(:)     :: a,b,c
    real(rp), intent(in   ), dimension(nx,ny) :: lambdaxy
    real(rp), intent(inout), dimension(:,:,:) :: p
    !
    real(rp), dimension(n) :: bb,p1,p2
    integer :: i,j,info
    !
    ! solve tridiagonal system when periodicity is imposed for pressure
    !
    !$OMP PARALLEL DEFAULT(none) &
    !$OMP PRIVATE(i,j,bb,p1,p2) &
    !$OMP SHARED(nx,ny,n,a,b,c,lambdaxy,p)
    !$OMP DO COLLAPSE(2)
    do j=1,ny
      do i=1,nx
        !
        bb(:) = b(:) + lambdaxy(i,j)
        p1(1:n-1) = p(i,j,1:n-1)
        call dgtsv_homebrewed(n-1,a(1:n-1),bb(1:n-1),c(1:n-1),p1(1:n-1))
        p2(:) = 0._rp
        p2(1  ) = -a(1  )
        p2(n-1) = -c(n-1)
        call dgtsv_homebrewed(n-1,a(1:n-1),bb(1:n-1),c(1:n-1),p2(1:n-1))
        p(i,j,n) = (p(i,j,n) - c(n)*p1(1) - a(n)*p1(n-1)) / &
                   (bb(   n) + c(n)*p2(1) + a(n)*p2(n-1))
        p(i,j,1:n-1) = p1(1:n-1) + p2(1:n-1)*p(i,j,n)
        !
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    !
    return
  end subroutine gaussel_periodic
  !
  subroutine dgtsv_homebrewed(n,a,b,c,p)
    !
    implicit none
    !
    integer , intent(in   )               :: n
    real(rp), intent(in   ), dimension(:) :: a,b,c
    real(rp), intent(inout), dimension(:) :: p
    !
    real(rp), dimension(n) :: d
    real(rp) :: z
    integer  :: l
    !
    ! Gauss elimination
    !
    z = 1._rp/b(1)
    d(1) = c(1)*z
    p(1) = p(1)*z
    do l=2,n-1
      z    = 1._rp/(b(l)-a(l)*d(l-1))
      d(l) = c(l)*z
      p(l) = (p(l)-a(l)*p(l-1))*z
    enddo
    z = b(n)-a(n)*d(n-1)
    if(z.ne.0._rp) then
      p(n) = (p(n)-a(n)*p(n-1))/z
    else
      p(n) = 0._rp
    endif
    !
    ! backward substitution
    !
    do l=n-1,1,-1
      p(l) = p(l) - d(l)*p(l+1)
    enddo
    !
    return
  end subroutine dgtsv_homebrewed
  !
end module mod_solver_cpu
