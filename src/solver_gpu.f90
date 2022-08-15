module mod_solver_gpu
  !
  use iso_c_binding , only: C_PTR
  use decomp_2d
  use mod_types
#if defined(_OPENACC)
  use cudafor
  use mod_fftw_param
  use mod_fft,        only: fftf_gpu, fftb_gpu, signal_processing, arr_tmp
  use mod_common_mpi, only: mydev
#else
  use mod_fft       , only: fft
#endif
  use profiler
  !
  implicit none
  !
  real(rp), allocatable, dimension(:,:,:) :: px,py,pz
#if defined(_OPENACC)
  real(rp), allocatable, dimension(:,:,:) :: pw
  !@cuf attributes(managed) :: px, py, pz, pw
  real(rp), allocatable, dimension(:,:,:), device :: py_t
  real(rp), allocatable, dimension(:,:,:), device :: pxc, pyc_t
#endif
  !
  private
  public  :: solver_gpu
  !
  contains
  !
  subroutine solver_gpu(n,dims,arrplan,normfft,lambdaxy,a,b,c,bc,c_or_f,p)
    !
    implicit none
    !
    integer         , intent(in   ), dimension(3)   :: n
    integer         , intent(in   ), dimension(3)   :: dims
    type(C_PTR)     , intent(in   ), dimension(2,2) :: arrplan
    real(rp)        , intent(in   )                 :: normfft
#if defined(_OPENACC)
    real(rp)        , intent(in   ), dimension(n(1)*dims(1)/dims(2),n(2)*dims(2)/dims(1)) :: lambdaxy
#else
    real(rp)        , intent(in   ), dimension(n(1),n(2)) :: lambdaxy
#endif
    real(rp)        , intent(in   ), dimension(n(3))      :: a,b,c
    character(len=1), intent(in   ), dimension(0:1,3)     :: bc
    character(len=1), intent(in   ), dimension(3)         :: c_or_f
    real(rp)        , intent(inout), dimension(0:,0:,0:)  :: p
    !
    integer :: i,j,k
    !@cuf attributes(managed) :: p, lambdaxy, a, b, c
    integer :: ii
    integer :: ng1, ng2, ng3
    integer :: n1, n2, n3
    integer :: dims1, dims2, dims3
    integer :: q, buff_dim1, buff_dim2, buff_dim3
    !@cuf integer :: istat
    !
    dims1 = dims(1) 
    dims2 = dims(2)
    dims3 = dims(3) 
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)
    ng1 = n1*dims1 
    ng2 = n2*dims2
    ng3 = n3 
    !
    if( .not. allocated(px) ) allocate(px(ng1,ng2/dims1,ng3/dims2))
    if( .not. allocated(py) ) allocate(py(ng1/dims1,ng2,ng3/dims2))
    if( .not. allocated(pz) ) allocate(pz(ng1/dims1,ng2/dims2,ng3))
#if defined(_OPENACC)
    if( .not. allocated(pw) ) allocate(pw(ng1/dims2,ng2/dims1,ng3))
    if( .not. allocated(pxc  )) allocate( pxc( 2*(ng1/2 + 1), ng2/dims1, ng3/dims2 ) )
    if( .not. allocated(pyc_t)) allocate( pyc_t( 2*(ng2/2 + 1), ng1/dims1, ng3/dims2 ) )
    if( .not. allocated(py_t )) then
      allocate( py_t( ng2, ng1/dims1, ng3/dims2 ) )
      !@cuf istat = cudaMemAdvise( px, size(px), cudaMemAdviseSetPreferredLocation, mydev )
      !@cuf istat = cudaMemAdvise( py, size(py), cudaMemAdviseSetPreferredLocation, mydev )
      !@cuf istat = cudaMemAdvise( pz, size(pz), cudaMemAdviseSetPreferredLocation, mydev )
      !@cuf istat = cudaMemAdvise( pw, size(pw), cudaMemAdviseSetPreferredLocation, mydev )
      !@cuf istat = cudaMemAdvise( pxc, size(pxc), cudaMemAdviseSetPreferredLocation, mydev )
      !@cuf istat = cudaMemAdvise( pyc_t, size(pyc_t), cudaMemAdviseSetPreferredLocation, mydev )
      !@cuf istat = cudaMemAdvise( py_t, size(py_t), cudaMemAdviseSetPreferredLocation, mydev )
      !!! !@cuf istat = cudaMemPrefetchAsync( px, size(px), mydev, 0)
      !!! !@cuf istat = cudaMemPrefetchAsync( py, size(py), mydev, 0)
      !!! !@cuf istat = cudaMemPrefetchAsync( pz, size(pz), mydev, 0)
      !!! !@cuf istat = cudaMemPrefetchAsync( pw, size(pw), mydev, 0)
      !!! !@cuf istat = cudaMemPrefetchAsync( p, size(p), mydev, 0)
      !!! !@cuf istat = cudaMemPrefetchAsync(lambdaxy,size(lambdaxy),mydev,0)
      !!! !@cuf istat = cudaMemPrefetchAsync( a, size(a), mydev, 0)
      !!! !@cuf istat = cudaMemPrefetchAsync( b, size(b), mydev, 0)
      !!! !@cuf istat = cudaMemPrefetchAsync( c, size(c), mydev, 0)
    endif
#endif
    !
#if defined(_OPENACC)
    !
    ! One big buffer allocated just once
    if(.not.allocated(arr_tmp)) then
      buff_dim1 = abs(max(ng1, ng2))+2
      buff_dim2 = abs(max(ng1/dims1, ng2/dims1))
      buff_dim3 = ng3/dims2
      allocate (arr_tmp(0:buff_dim1-1,buff_dim2,buff_dim3))
    endif
    !
#if defined(_DECOMP_X) || defined(_DECOMP_Y)
    if(product(dims(1:1)).eq.1) then
#else
    if(product(dims(1:2)).eq.1) then
#endif
      !
      ! case 1: only 1 GPU
      !
      ! transpose to have y values in the leading dimension
      !
      !$acc parallel loop collapse(3) present(py_t,p)
      do k=1,ng3/dims2
        do j=1,ng1/dims1
          do i=1,ng2
            py_t(i,j,k) = p(j,i,k)
          enddo
        enddo
      enddo
      !$acc end parallel loop
      !
      call signal_processing(0,'F',bc(0,2)//bc(1,2),c_or_f(2),(/ ng2  ,ng1/dims1,ng3/dims2 /),1,py_t )
      call fftf_gpu(cufft_plan_fwd_y, py_t, pyc_t)
      call signal_processing(1,'F',bc(0,2)//bc(1,2),c_or_f(2),(/ ng2+2,ng1/dims1,ng3/dims2 /),1,pyc_t)
      !@cuf istat = cudaDeviceSynchronize()
      !
      ! transpose to have x values in the leading dimension
      !
      !$acc parallel loop collapse(3) present(px,pyc_t)
      do k=1,ng3/dims2
        do i=1,ng2
          do j=1,ng1/dims1
            px(j,i,k) = pyc_t(i,j,k)
          enddo
        enddo
      enddo
      !$acc end parallel loop
      !
      call signal_processing(0,'F',bc(0,1)//bc(1,1),c_or_f(1),(/ ng1  ,ng2/dims1,ng3/dims2 /),1,px )
      call fftf_gpu(cufft_plan_fwd_x, px, pxc)
      call signal_processing(1,'F',bc(0,1)//bc(1,1),c_or_f(1),(/ ng1+2,ng2/dims1,ng3/dims2 /),1,pxc)
      !@cuf istat = cudaDeviceSynchronize()
      !
      ! dummy operation
      !
#if defined(_DECOMP_X) || defined(_DECOMP_Y)
      call profiler_start("xc_to_z", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_xc_to_z(px,pw,pxc)
      call profiler_stop("xc_to_z")
#else
      !
      !$acc parallel loop collapse(3) present(pw,pxc)
      do k=1,ng3/dims2
        do j=1,ng2/dims1
          do i=1,ng1
            pw(i,j,k) = pxc(i,j,k)
          enddo
        enddo
      enddo
      !$acc end parallel loop
#endif
      !
      q = 0
      if(c_or_f(3).eq.'f'.and.bc(1,3).eq.'D') q = 1
      !
      if(bc(0,3)//bc(1,3).eq.'PP') then
        call gaussel_periodic_gpu(ng1/dims2,ng2/dims1,n3-q,a,b,c,lambdaxy,pw,px,py,pxc,pyc_t)
      else
        call gaussel_gpu(         ng1/dims2,ng2/dims1,n3-q,a,b,c,lambdaxy,pw,px,py)
      endif
      !@cuf istat = cudaDeviceSynchronize()
      !
      ! dummy operation
      !
#if defined(_DECOMP_X) || defined(_DECOMP_Y)
      call profiler_start("z_to_xc", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_z_to_xc(pw,px,pxc)
      call profiler_stop("z_to_xc")
#else
      !
      !$acc parallel loop collapse(3) present(pxc,pw)
      do k=1,ng3/dims2
        do j=1,ng2/dims1
          do i=1,ng1
            pxc(i,j,k) = pw(i,j,k)
          enddo
        enddo
      enddo
      !$acc end parallel loop
#endif
      !
      call signal_processing(0,'B',bc(0,1)//bc(1,1),c_or_f(1),(/ ng1+2,ng2/dims1,ng3/dims2 /),1,pxc)
      call fftb_gpu(cufft_plan_bwd_x, pxc, px)
      call signal_processing(1,'B',bc(0,1)//bc(1,1),c_or_f(1),(/ ng1  ,ng2/dims1,ng3/dims2 /),1,px )
      !@cuf istat = cudaDeviceSynchronize()
      !
      ! transpose to have y values in the leading dimension
      !
      !$acc parallel loop collapse(3) present(pyc_t,px)
      do k=1,ng3/dims2
        do j=1,ng1/dims1
          do i=1,ng2
            pyc_t(i,j,k) = px(j,i,k)
          enddo
        enddo
      enddo
      !$acc end parallel loop
      !
      call signal_processing(0,'B',bc(0,2)//bc(1,2),c_or_f(2),(/ ng2+2,ng1/dims1,ng3/dims2 /),1,pyc_t)
      call fftb_gpu(cufft_plan_bwd_y, pyc_t, py_t)
      call signal_processing(1,'B',bc(0,2)//bc(1,2),c_or_f(2),(/ ng2  ,ng1/dims1,ng3/dims2 /),1,py_t )
      !@cuf istat = cudaDeviceSynchronize()
      !
      ! transpose to have x values in the leading dimension
      !
      !$acc parallel loop collapse(3) present(py_t,p)
      do k=1,ng3/dims2
        do j=1,ng2
          do i=1,ng1/dims1
            p(i,j,k) = py_t(j,i,k)*normfft
          enddo
        enddo
      enddo
      !$acc end parallel loop
      !
    else 
#endif
      !
      ! case 2: many GPUs
      !
      ! [NOTE] this is an interesting loop, why is it swapped this way?
#if !defined(_EPHC)
      !$acc parallel loop collapse(3) present(pz,p)
      do k=1,ng3
        do j=1,ng2/dims2
          do i=1,ng1/dims1
            pz(i,j,k) = p(i,j,k)
          enddo
        enddo
      enddo
      !$acc end parallel loop
      !
      call profiler_start("z_to_y", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_z_to_y(pz,py)
      call profiler_stop("z_to_y")
#endif
      !
#if defined(_OPENACC)
#if !defined(_EPHC)
      !$acc parallel loop collapse(3) present(py_t,py)
      do k=1,ng3/dims2
        do j=1,ng1/dims1
          do i=1,ng2
            py_t(i,j,k) = py(j,i,k)
          enddo
        enddo
      enddo
      !$acc end parallel loop
#else
      call profiler_start("zp_to_yt", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_zp_to_yt(pz,py,p,py_t)
      call profiler_stop("zp_to_yt")
#endif
      !
      call signal_processing(0,'F',bc(0,2)//bc(1,2),c_or_f(2),(/ ng2  ,ng1/dims1,ng3/dims2 /),1,py_t )
      call fftf_gpu(cufft_plan_fwd_y, py_t, pyc_t)
      call signal_processing(1,'F',bc(0,2)//bc(1,2),c_or_f(2),(/ ng2+2,ng1/dims1,ng3/dims2 /),1,pyc_t)
      !@cuf istat = cudaDeviceSynchronize()
      !
#if !defined(_EPHC)
      !$acc parallel loop collapse(3) present(pyc_t,py)
      do k=1,ng3/dims2
        do i=1,ng2
          do j=1,ng1/dims1
            py(j,i,k) = pyc_t(i,j,k)
          enddo
        enddo
      enddo
      !$acc end parallel loop
#endif
      !
#else
      call fft(arrplan(1,2),py) ! fwd transform in y
#endif
      !
#if !defined(_EPHC)
      call profiler_start("y_to_x", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_y_to_x(py,px)
      call profiler_stop("y_to_x")
#endif
      !
#if defined(_OPENACC)
#if defined(_EPHC)
      call profiler_start("yct_to_x", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_yct_to_x(py,px,pyc_t)
      call profiler_stop("yct_to_x")
#endif
      !
      call signal_processing(0,'F',bc(0,1)//bc(1,1),c_or_f(1),(/ ng1  ,ng2/dims1,ng3/dims2 /),1,px )
      call fftf_gpu(cufft_plan_fwd_x, px, pxc)
      call signal_processing(1,'F',bc(0,1)//bc(1,1),c_or_f(1),(/ ng1+2,ng2/dims1,ng3/dims2 /),1,pxc)
      !@cuf istat = cudaDeviceSynchronize()
      !
#if !defined(_EPHC)
      !$acc parallel loop collapse(3) present(px,pxc)
      do k=1,ng3/dims2
        do j=1,ng2/dims1
          do i=1,ng1
            px(i,j,k) = pxc(i,j,k)
          enddo
        enddo
      enddo
      !$acc end parallel loop
#endif
      !
#else
      call fft(arrplan(1,1),px) ! fwd transform in x
#endif
      !
#if defined(_OPENACC)
#if !defined(_EPHC)
      call profiler_start("x_to_z", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_x_to_z(px,pw)
      call profiler_stop("x_to_z")
#else
      call profiler_start("xc_to_z", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_xc_to_z(px,pw,pxc)
      call profiler_stop("xc_to_z")
#endif
#else
      call profiler_start("x_to_y", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_x_to_y(px,py)
      call profiler_stop("x_to_y")
      call profiler_start("y_to_z", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_y_to_z(py,pz)
      call profiler_stop("y_to_z")
#endif
      !
      q = 0
      if(c_or_f(3).eq.'f'.and.bc(1,3).eq.'D') q = 1
      if(bc(0,3)//bc(1,3).eq.'PP') then
#if defined(_OPENACC)
        call gaussel_periodic_gpu(ng1/dims2,ng2/dims1,n3-q,a,b,c,lambdaxy,pw,px,py,pxc,pyc_t)
#else
        call gaussel_periodic(n1,n2,n3-q,a,b,c,lambdaxy,pz)
#endif
      else
#if defined(_OPENACC)
        call gaussel_gpu(         ng1/dims2,ng2/dims1,n3-q,a,b,c,lambdaxy,pw,px,py)
#else
        call gaussel(         n1,n2,n3-q,a,b,c,lambdaxy,pz)
#endif
      endif
      !@cuf istat = cudaDeviceSynchronize()
      !
#if defined(_OPENACC)
#if !defined(_EPHC)
      call profiler_start("z_to_x", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_z_to_x(pw,px)
      call profiler_stop("z_to_x")
#else
      call profiler_start("z_to_xc", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_z_to_xc(pw,px,pxc)
      call profiler_stop("z_to_xc")
#endif
#else
      call profiler_start("z_to_y", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_z_to_y(pz,py)
      call profiler_stop("z_to_y")
      call profiler_start("y_to_x", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_y_to_x(py,px)
      call profiler_stop("y_to_x")
#endif
      !
#if defined(_OPENACC)
#if !defined(_EPHC)
      !$acc parallel loop collapse(3) present(pxc,px)
      do k=1,ng3/dims2
        do j=1,ng2/dims1
          do i=1,ng1
            pxc(i,j,k) = px(i,j,k)
          enddo
        enddo
      enddo
      !$acc end parallel loop
#endif
      call signal_processing(0,'B',bc(0,1)//bc(1,1),c_or_f(1),(/ ng1+2,ng2/dims1,ng3/dims2 /),1,pxc)
      call fftb_gpu(cufft_plan_bwd_x, pxc, px)
      call signal_processing(1,'B',bc(0,1)//bc(1,1),c_or_f(1),(/ ng1  ,ng2/dims1,ng3/dims2 /),1,px )
      !@cuf istat = cudaDeviceSynchronize()
      !
#else
      call fft(arrplan(2,1),px) ! bwd transform in x
#endif
      !
#if !defined(_EPHC)
      call profiler_start("x_to_y", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_x_to_y(px,py)
      call profiler_stop("x_to_y")
#else
      call profiler_start("x_to_yct", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_x_to_yct(px,py,pyc_t)
      call profiler_stop("x_to_yct")
#endif
      !
#if defined(_OPENACC)
      !
#if !defined(_EPHC)
      !$acc parallel loop collapse(3) present(pyc_t,py)
      do k=1,ng3/dims2
        do j=1,ng1/dims1
          do i=1,ng2
            pyc_t(i,j,k) = py(j,i,k)
          enddo
        enddo
      enddo
      !$acc end parallel loop
#endif
      !
      call signal_processing(0,'B',bc(0,2)//bc(1,2),c_or_f(2),(/ ng2+2,ng1/dims1,ng3/dims2 /),1,pyc_t)
      call fftb_gpu(cufft_plan_bwd_y, pyc_t, py_t)
      call signal_processing(1,'B',bc(0,2)//bc(1,2),c_or_f(2),(/ ng2  ,ng1/dims1,ng3/dims2 /),1,py_t)
      !@cuf istat = cudaDeviceSynchronize()
      !
#if !defined(_EPHC)
      !$acc parallel loop collapse(3) present(py_t,py)
      do k=1,ng3/dims2
        do i=1,ng2
          do j=1,ng1/dims1
            py(j,i,k) = py_t(i,j,k)
          enddo
        enddo
      enddo
      !$acc end parallel loop
#endif
#else
      call fft(arrplan(2,2),py) ! bwd transform in y
#endif
      !
#if !defined(_EPHC)
      call profiler_start("y_to_z", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_y_to_z(py,pz)
      call profiler_stop("y_to_z")
#else
      call profiler_start("yt_to_zp", tag = .true., tag_color = COLOR_PURPLE)
      call transpose_yt_to_zp( py, pz, py_t, p, normfft)
      call profiler_stop("yt_to_zp")
#endif
      !
#if !defined(_EPHC)
      !$acc parallel loop collapse(3) present(p, pz)
      do k=1,ng3
        do j=1,ng2/dims2
          do i=1,ng1/dims1
            p(i,j,k) = pz(i,j,k)*normfft
          enddo
        enddo
      enddo
      !$acc end parallel loop
      !
#endif
#if defined(_OPENACC)
    endif
#endif
    !
    return
    !
  end subroutine solver_gpu
  !
#if defined(_OPENACC)
  subroutine gaussel_gpu(nx,ny,n,a,b,c,lambdaxy,p,bb,d)
    !
    implicit none
    !
    integer , intent(in   )                     :: nx,ny,n
    real(rp), intent(in   ), dimension(:)       :: a,b,c
    real(rp), intent(in   ), dimension(nx,ny)   :: lambdaxy
    real(rp), intent(inout), dimension(:,:,:)   :: p
    real(rp), intent(inout), dimension(nx,ny,n) :: bb
    real(rp), intent(inout), dimension(nx,ny,n) :: d
    !
    !@cuf attributes(managed) :: a, b, c, lambdaxy, p
    !@cuf attributes(device)  :: bb, d
    real(rp) :: z
    integer  :: i,j,k
    !
    ! solve tridiagonal system
    !
    !$acc parallel loop collapse(2) present(d, bb)
    do j=1,ny
      do i=1,nx
        !
        !$acc loop seq
        do k=1,n
          bb(i,j,k) = b(k) + lambdaxy(i,j)
        enddo
        !$acc end loop
        z = 1._rp/bb(i,j,1)
        d(i,j,1) = c(1)*z
        p(i,j,1) = p(i,j,1)*z
        !$acc loop seq
        do k=2,n-1
          z = 1._rp/(bb(i,j,k)-a(k)*d(i,j,k-1))
          d(i,j,k) = c(k)*z
          p(i,j,k) = (p(i,j,k)-a(k)*p(i,j,k-1))*z
        enddo
        !$acc end loop
        z = bb(i,j,n)-a(n)*d(i,j,n-1)
        if(z.ne.0._rp) then
          p(i,j,n) = (p(i,j,n)-a(n)*p(i,j,n-1))/z
        else
          p(i,j,n) = 0._rp
        endif
        !
        ! backward substitution
        !
        !$acc loop seq
        do k=n-1,1,-1
          p(i,j,k) = p(i,j,k) - d(i,j,k)*p(i,j,k+1)
        enddo
        !$acc end loop
        !
      enddo
    enddo
    !$acc end parallel loop
    !
    return
  end subroutine gaussel_gpu
#endif
  !
#if defined(_OPENACC)
  subroutine gaussel_periodic_gpu(nx,ny,n,a,b,c,lambdaxy,p,bb,d,p1,p2)
    !
    implicit none
    !
    integer , intent(in   )                     :: nx,ny,n
    real(rp), intent(in   ), dimension(:)       :: a, b, c
    real(rp), intent(in   ), dimension(nx,ny)   :: lambdaxy
    real(rp), intent(inout), dimension(:,:,:)   :: p
    real(rp), intent(inout), dimension(nx,ny,n) :: bb, d, p1, p2
    !
    !@cuf attributes(managed) :: a, b, c, lambdaxy, p
    !@cuf attributes(device)  :: bb, d, p1, p2
    real(rp) :: z
    integer  :: i,j,k
    !
    ! solve tridiagonal system
    !
    !$acc parallel loop collapse(2) present(d, bb, p1, p2)
    do j=1,ny
      do i=1,nx
        !
        !$acc loop seq
        do k=1,n
          bb(i,j,k) = b(k) + lambdaxy(i,j)
        enddo
        !$acc end loop
        !$acc loop seq
        do k=1,n-1
          p1(i,j,k) = p(i,j,k)
        enddo
        ! call dgtsv_homebrewed(n-1,a(1:n-1),bb(1:n-1),c(1:n-2),p1(1:n-1))
        z = 1._rp/bb(i,j,1)
        d(i,j,1) = c(1)*z
        p1(i,j,1) = p1(i,j,1)*z
        !$acc loop seq
        do k=2,n-2
          z = 1._rp/(bb(i,j,k)-a(k)*d(i,j,k-1))
          d(i,j,k) = c(k)*z
          p1(i,j,k) = (p1(i,j,k)-a(k)*p1(i,j,k-1))*z
        enddo
        !$acc end loop
        z = bb(i,j,n-1)-a(n-1)*d(i,j,n-2)
        if(z.ne.0.) then
          p1(i,j,n-1) = (p1(i,j,n-1)-a(n-1)*p1(i,j,n-2))/z
        else
          p1(i,j,n-1) = 0._rp
        endif
        !
        ! backward substitution
        !
        !$acc loop seq
        do k=n-2,1,-1
          p1(i,j,k) = p1(i,j,k) - d(i,j,k)*p1(i,j,k+1)
        enddo
        !$acc end loop
        !
        !$acc loop seq
        do k=1,n
          p2(i,j,k) = 0._rp
        enddo
        !$acc end loop
        !
        p2(i,j,1  ) = -a(1  )
        p2(i,j,n-1) = -c(n-1)
        !call dgtsv_homebrewed(n-1,a(2:n-1),bb(1:n-1),c(1:n-2),p2(1:n-1))
        z = 1._rp/bb(i,j,1)
        d(i,j,1) = c(1)*z
        p2(i,j,1) = p2(i,j,1)*z
        !$acc loop seq
        do k=2,n-2
          z = 1._rp/(bb(i,j,k)-a(k+1)*d(i,j,k-1))
          d(i,j,k) = c(k)*z
          p2(i,j,k) = (p2(i,j,k)-a(k+1)*p2(i,j,k-1))*z
        enddo
        !$acc end loop
        z = bb(i,j,n-1)-a(n)*d(i,j,n-2)
        if(z.ne.0.) then
          p2(i,j,n-1) = (p2(i,j,n-1)-a(n)*p2(i,j,n-2))/z
        else
          p2(i,j,n-1) = 0._rp
        endif
        !
        ! backward substitution
        !
        !$acc loop seq
        do k=n-2,1,-1
          p2(i,j,k) = p2(i,j,k) - d(i,j,k)*p2(i,j,k+1)
        enddo
        !$acc end loop
        p(i,j,n) = (p(i,j,n) - c(n)*p1(i,j,1) - a(n)*p1(i,j,n-1)) / &
                   (bb(i,j,n) + c(n)*p2(i,j,1) + a(n)*p2(i,j,n-1))
        !$acc loop seq
        do k=1,n-1
          p(i,j,k) = p1(i,j,k) + p2(i,j,k)*p(i,j,n)
        enddo
        !$acc end loop
        !
      enddo
    enddo
    !$acc end parallel loop
    !
    return
  end subroutine gaussel_periodic_gpu
#endif
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
        !
        bb(:) = b(1:n) + lambdaxy(i,j)
        call dgtsv_homebrewed(n,a,bb,c,p(i,j,1:n))
        !
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
        bb(:)  = b(:) + lambdaxy(i,j)
        p1(1:n-1) = p(i,j,1:n-1)
        call dgtsv_homebrewed(n-1,a(1:n-1),bb(1:n-1),c(1:n-1),p1(1:n-1))
        p2(:)   = 0._rp
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
end module mod_solver_gpu
