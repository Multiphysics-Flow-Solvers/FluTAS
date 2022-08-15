!
! SPDX-License-Identifier: MIT
!
module mod_fft
  !
  use iso_c_binding , only: C_INT
  use mod_common_mpi, only: ierr,myid
  use mod_fftw_param
  use mod_param     , only: pi
  use mod_types
  !@cuf use cudafor
  !
  !$ use omp_lib
  real(rp), allocatable, dimension(:,:,:) :: arr_tmp
  !@cuf attributes(device) :: arr_tmp
  !
  public  :: fftini,fftend,fft
#if defined(_OPENACC)
  public  :: fftf_gpu,fftb_gpu,signal_processing
#endif
  !
  contains
  ! 
  subroutine fftini(n_x,n_y,bcxy,c_or_f,arrplan,normfft)
    !
    implicit none
    !
    integer         , intent(in ), dimension(3)     :: n_x,n_y
    character(len=1), intent(in ), dimension(0:1,2) :: bcxy
    character(len=1), intent(in ), dimension(2)     :: c_or_f
    type(C_PTR)     , intent(out), dimension(2,2)   :: arrplan
    real(rp)        , intent(out)                   :: normfft
    !
    real(rp), dimension(n_x(1),n_x(2),n_x(3))  :: arrx
    real(rp), dimension(n_y(1),n_y(2),n_y(3))  :: arry
    type(C_PTR) :: plan_fwd_x,plan_bwd_x, &
                   plan_fwd_y,plan_bwd_y
    type(fftw_iodim), dimension(1) :: iodim
    type(fftw_iodim), dimension(2) :: iodim_howmany
    integer :: kind_fwd,kind_bwd
    real(rp), dimension(2) :: norm
    integer(C_INT) :: nx_x,ny_x,nz_x, &
                      nx_y,ny_y,nz_y
    integer :: ix,iy
#if defined(_OPENACC)
    integer :: istat
    integer(int_ptr_kind()) :: worksize, max_worksize
    integer, pointer :: null_fptr
    call c_f_pointer( c_null_ptr, null_fptr )
    max_worksize = 0
#endif
#if defined(_SINGLE_PRECISION)
    !$ call sfftw_init_threads(ierr)
    !$ call sfftw_plan_with_nthreads(omp_get_max_threads())
#else
    !$ call dfftw_init_threads(ierr)
    !$ call dfftw_plan_with_nthreads(omp_get_max_threads())
#endif
    !
    ! fft in x
    !
    ! prepare plans with guru interface
    !
    nx_x = n_x(1)
    ny_x = n_x(2)
    nz_x = n_x(3)
    nx_y = n_y(1)
    ny_y = n_y(2)
    nz_y = n_y(3)
    !
    normfft = 1._rp
    ix = 0
    ! size of transform reduced by 1 point with Dirichlet BC in face
    if (bcxy(0, 1)//bcxy(1, 1).eq.'DD' .and. c_or_f(1).eq.'f') ix = 1
    iodim(1)%n = nx_x - ix
    iodim(1)%is = 1
    iodim(1)%os = 1
    iodim_howmany(1)%n = ny_x
    iodim_howmany(1)%is = nx_x
    iodim_howmany(1)%os = nx_x
    iodim_howmany(2)%n = nz_x
    iodim_howmany(2)%is = nx_x*ny_x
    iodim_howmany(2)%os = nx_x*ny_x
    call find_fft(bcxy(:, 1), c_or_f(1), kind_fwd, kind_bwd, norm)
    plan_fwd_x = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, arrx, arrx, kind_fwd, FFTW_ESTIMATE)
    plan_bwd_x = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, arrx, arrx, kind_bwd, FFTW_ESTIMATE)
    normfft = normfft*norm(1)*(nx_x + norm(2) - ix)
#if defined(_OPENACC)
    if(.not. allocated(cufft_workspace)) then
      batch = ny_x*nz_x
      !istat = cufftPlan1D(cufft_plan_fwd_x,nx_x,CUFFT_FWD_TYPE,batch)
      istat = cufftCreate(cufft_plan_fwd_x)
      istat = cufftSetAutoAllocation(cufft_plan_fwd_x, 0)
      ! Then very first time this routine is caslled it consumes a lot of time... why?
      istat = cufftMakePlanMany(cufft_plan_fwd_x, 1, nx_x, null_fptr, 1, nx_x, null_fptr, 1, nx_x, CUFFT_FWD_TYPE, batch, worksize)
      max_worksize = max(worksize, max_worksize)
      !
      !istat = cufftPlan1D(cufft_plan_bwd_x,nx_x,CUFFT_BWD_TYPE,batch);
      istat = cufftCreate(cufft_plan_bwd_x)
      istat = cufftSetAutoAllocation(cufft_plan_bwd_x, 0)
      istat = cufftMakePlanMany(cufft_plan_bwd_x, 1, nx_x, null_fptr, 1, nx_x, null_fptr, 1, nx_x, CUFFT_BWD_TYPE, batch, worksize)
      max_worksize = max(worksize, max_worksize)
    end if
#endif
    !
    ! fft in y
    !
    ! prepare plans with guru interface
    !
    iy = 0
    ! size of transform reduced by 1 point with Dirichlet BC in face
    if (bcxy(0, 2)//bcxy(1, 2).eq.'DD' .and. c_or_f(2).eq.'f') iy = 1
    iodim(1)%n = ny_y - iy
    iodim(1)%is = nx_y
    iodim(1)%os = nx_y
    iodim_howmany(1)%n = nx_y
    iodim_howmany(1)%is = 1
    iodim_howmany(1)%os = 1
    iodim_howmany(2)%n = nz_y
    iodim_howmany(2)%is = nx_y*ny_y
    iodim_howmany(2)%os = nx_y*ny_y
    call find_fft(bcxy(:, 2), c_or_f(2), kind_fwd, kind_bwd, norm)
    plan_fwd_y = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, arry, arry, kind_fwd, FFTW_ESTIMATE)
    plan_bwd_y = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, arry, arry, kind_bwd, FFTW_ESTIMATE)
    normfft = normfft*norm(1)*(ny_y + norm(2) - iy)
#if defined(_OPENACC)
    if(.not. allocated(cufft_workspace)) then
      batch = nx_y*nz_y
      !istat = cufftPlan1D(cufft_plan_fwd_y,ny_y,CUFFT_FWD_TYPE,batch)
      istat = cufftCreate(cufft_plan_fwd_y)
      istat = cufftSetAutoAllocation(cufft_plan_fwd_y, 0)
      istat = cufftMakePlanMany(cufft_plan_fwd_y, 1, ny_y, null_fptr, 1, ny_y, null_fptr, 1, ny_y, CUFFT_FWD_TYPE, batch, worksize)
      max_worksize = max(worksize, max_worksize)
      !
      !istat = cufftPlan1D(cufft_plan_bwd_y,ny_y,CUFFT_BWD_TYPE,batch)
      istat = cufftCreate(cufft_plan_bwd_y)
      istat = cufftSetAutoAllocation(cufft_plan_bwd_y, 0)
      istat = cufftMakePlanMany(cufft_plan_bwd_y, 1, ny_y, null_fptr, 1, ny_y, null_fptr, 1, ny_y, CUFFT_BWD_TYPE, batch, worksize)
      max_worksize = max(worksize, max_worksize)
      !
      allocate (cufft_workspace(max_worksize/(2*sizeof(1._rp))))
      !
      istat = cufftSetWorkArea(cufft_plan_fwd_x, cufft_workspace)
      istat = cufftSetWorkArea(cufft_plan_bwd_x, cufft_workspace)
      istat = cufftSetWorkArea(cufft_plan_fwd_y, cufft_workspace)
      istat = cufftSetWorkArea(cufft_plan_bwd_y, cufft_workspace)
    endif
#endif
    !
    normfft = normfft**(-1)
    arrplan(1, 1) = plan_fwd_x
    arrplan(2, 1) = plan_bwd_x
    arrplan(1, 2) = plan_fwd_y
    arrplan(2, 2) = plan_bwd_y
    !
    return
  end subroutine fftini
  !
  subroutine fftend(arrplan)
    !
    implicit none
    type(C_PTR), intent(in), dimension(2,2) :: arrplan
    !
#if defined(_SINGLE_PRECISION)
    call sfftw_destroy_plan(arrplan(1, 1))
    call sfftw_destroy_plan(arrplan(1, 2))
    call sfftw_destroy_plan(arrplan(2, 1))
    call sfftw_destroy_plan(arrplan(2, 2))
    !$call sfftw_cleanup_threads(ierr)
#else
    call dfftw_destroy_plan(arrplan(1, 1))
    call dfftw_destroy_plan(arrplan(1, 2))
    call dfftw_destroy_plan(arrplan(2, 1))
    call dfftw_destroy_plan(arrplan(2, 2))
    !$call dfftw_cleanup_threads(ierr)
#endif
    !
    return
  end subroutine fftend
  !
  subroutine fft(plan,arr)
    !
    implicit none
    !
    type(C_PTR), intent(in   )                   :: plan
    real(rp)   , intent(inout), dimension(:,:,:) :: arr
#if defined(_SINGLE_PRECISION)
    call sfftw_execute_r2r(plan, arr, arr)
#else
    call dfftw_execute_r2r(plan, arr, arr)
#endif
    return
  end subroutine fft
  !
#if defined(_OPENACC)
  subroutine fftf_gpu(plan, arrin, arrout)
    !
    implicit none
    !
    integer , intent(in )                           :: plan
    real(rp), intent(in ), dimension(:,:,:), device :: arrin
    real(rp), intent(out), dimension(:,:,:), device :: arrout
    !
    integer :: istat
#if defined(_SINGLE_PRECISION)
    istat = cufftExecR2C(plan, arrin, arrout)
#else
    istat = cufftExecD2Z(plan, arrin, arrout)
#endif
    !
    return
  end subroutine fftf_gpu
  !
  subroutine fftb_gpu(plan, arrin, arrout)
    !
    implicit none
    !
    integer , intent(in )                           :: plan
    real(rp), intent(in ), dimension(:,:,:), device :: arrin
    real(rp), intent(out), dimension(:,:,:), device :: arrout
    !
    integer :: istat
#if defined(_SINGLE_PRECISION)
    istat = cufftExecC2R(plan, arrin, arrout)
#else
    istat = cufftExecZ2D(plan, arrin, arrout)
#endif
    !
    return
  end subroutine fftb_gpu
#endif
  !
  subroutine find_fft(bc,c_or_f,kind_fwd,kind_bwd,norm)
    !
    implicit none
    !
    character(len=1), intent(in ), dimension(0:1) :: bc
    character(len=1), intent(in )                 :: c_or_f
    integer         , intent(out)                 :: kind_fwd, kind_bwd
    real(rp)        , intent(out), dimension(2)   :: norm
    !
    if(c_or_f.eq.'c') then
      select case (bc(0)//bc(1))
      case ('PP')
         kind_fwd = FFTW_R2HC
         kind_bwd = FFTW_HC2R
         norm = (/1._rp, 0._rp/)
      case ('NN')
         kind_fwd = FFTW_REDFT10
         kind_bwd = FFTW_REDFT01
         norm = (/2._rp, 0._rp/)
      case ('DD')
         kind_fwd = FFTW_RODFT10
         kind_bwd = FFTW_RODFT01
         norm = (/2._rp, 0._rp/)
      case ('ND')
         kind_fwd = FFTW_REDFT11
         kind_bwd = FFTW_REDFT11
         norm = (/2._rp, 0._rp/)
      case ('DN')
         kind_fwd = FFTW_RODFT11
         kind_bwd = FFTW_RODFT11
         norm = (/2._rp, 0._rp/)
      end select
    elseif (c_or_f.eq.'f') then
      select case (bc(0)//bc(1))
      case ('PP')
         kind_fwd = FFTW_R2HC
         kind_bwd = FFTW_HC2R
         norm = (/1._rp, 0._rp/)
      case ('NN')
         kind_fwd = FFTW_REDFT00
         kind_bwd = FFTW_REDFT00
         norm = (/2._rp, -1._rp/)
      case ('DD')
         kind_fwd = FFTW_RODFT00
         kind_bwd = FFTW_RODFT00
         norm = (/2._rp, 1._rp/)
      case ('ND')
         kind_fwd = FFTW_REDFT10
         kind_bwd = FFTW_REDFT01
         norm = (/2._rp, 0._rp/)
      case ('DN')
         kind_fwd = FFTW_RODFT01
         kind_bwd = FFTW_RODFT10
         norm = (/2._rp, 0._rp/)
      end select
    end if
    !
    return
  end subroutine find_fft
  !
#if defined(_OPENACC)
  subroutine posp_fftf(n1, n2, n3, idir, arr)
    !
    ! post-processing of a signal following a forward FFT
    ! to order the data as follows:
    ! (r[0],r[n],r[1],i[1],...,r[n-1],i[n-1])
    !
    implicit none
    !
    integer, intent(in) :: n1, n2, n3       ! dimensions of input/output array
    integer , intent(in   )                   :: idir ! direction where the transform is taken
    real(rp), intent(inout), dimension(:,:,:) :: arr  ! input/output array
    integer :: j, k, nn
!@cuf attributes(device) :: arr
    !
    !nn = n(idir) - 2
    select case(idir)
    case(1)
      nn = n1 - 2
      !$acc parallel loop collapse(2) present(arr)
      do k=1,n3
        do j=1,n2
          arr(2,j,k) = arr(nn+1,j,k)
         enddo
      enddo
      !$acc end parallel loop
      !
    end select
    !
    return
  end subroutine posp_fftf
  !
  subroutine prep_fftb(n1, n2, n3, idir, arr)
    !
    ! pre-processing of a signal preciding a backward FFT
    ! to order the data as follows:
    ! (r[0],i[0],r[1],i[1],...,r[n-1],i[n-1],r[n],i[n])
    ! note that i[0] = i[n] = 0
    !
    implicit none
    !
    integer, intent(in) :: n1, n2, n3       ! dimensions of input/output array
    integer , intent(in   )                   :: idir ! direction where the transform is taken
    real(rp), intent(inout), dimension(:,:,:) :: arr  ! input/output array
    !
    integer :: j,k,nn
    !@cuf attributes(device) :: arr
    !
    !nn = n(idir) - 2
    select case (idir)
    case(1)
      nn = n1 - 2
      !$acc parallel loop collapse(2) present(arr)
      do k=1,n3
        do j=1,n2
          arr(nn+1,j,k) = arr(2,j,k)
          arr(2,j,k)    = 0._rp
        end do
      end do
      !$acc end parallel loop
      !
    end select
    !
    return
  end subroutine prep_fftb
  !
  subroutine prep_dctiif(n1, n2, n3,idir,arr,is_swap_order,is_negate_even)
    !
    ! pre-processing of a signal to perform a fast forward
    ! discrete cosine transform (DCT) with FFTs (see Makhoul 1980)
    !
    ! the input signal x(n) is pre-processed into a signal v(n)
    ! as follows:
    !
    ! v(n) = x(2n       ),              0 <= n <= floor((N-1)/2)
    !      = x(2N -2n -1), floor((N+1)/2) <= n <= N-1
    ! with n = 0,...,N-1 and N is the total number of elements
    ! of the signal.
    !
    ! pre-processing required for computing the corresponding
    ! discrete sine transform (DST) may also be performed
    ! if one of the last boolean input variables are .true.
    !
    implicit none
    !
    integer, intent(in) :: n1, n2, n3       ! dimensions of input/output array
    integer , intent(in   )                   :: idir                  ! array direction where the transform is taken
    real(rp), intent(inout), dimension(:,:,:) :: arr ! input/output array
    logical , intent(in   )                   :: is_swap_order  ! swap order of the elements of the input array? (for DST)
    logical , intent(in   )                   :: is_negate_even ! negate every other element of the input array?
    !
    integer :: i,j,k,nn,ii
    !@cuf attributes(device) :: arr
    !
    !nn = n(idir)
    select case (idir)
    case(1)
      nn = n1
      if(.not.allocated(arr_tmp)) allocate (arr_tmp(0:n1-1,n2,n3))
      !
      if(is_swap_order) then
        !$acc parallel loop collapse(3) present(arr_tmp, arr)
        do k=1,n3
          do j=1,n2
            do i =1,nn
              ii =i-1
              arr_tmp(ii,j,k) = arr(i,j,k)
            end do
          end do
        end do
        !$acc end parallel loop
        !
        !$acc parallel loop collapse(3) present(arr_tmp, arr)
        do k=1,n3
          do j =1,n2
            do i =1,nn
              ii = i-1
              arr(i,j,k) = arr_tmp(nn-1-ii,j,k)
            end do
          end do
        end do
        !$acc end parallel loop
      end if
      !
      if(is_negate_even) then
        !$acc parallel loop collapse(3) present(arr)
        do k=1,n3
          do j=1,n2
            do i=1,nn
              if(mod(i,2).eq.0) then
                arr(i,j,k) = -arr(i,j,k)
              end if
            end do
          end do
        end do
        !$acc end parallel loop
      end if
      !
      !$acc parallel loop collapse(3) present(arr_tmp, arr)
      do k = 1, n3
         do j = 1, n2
            do i = 1, n1
               ii = i - 1
               if (ii.le.(nn - 1)/2) then
                  arr_tmp(ii, j, k) = arr(2*ii + 1, j, k)
               else
                  arr_tmp(ii, j, k) = arr(2*(nn - ii) - 1 + 1, j, k)
               end if
            end do
         end do
      end do
      !$acc end parallel loop
      !
      !$acc parallel loop collapse(3) present(arr_tmp, arr)
      do k = 1, n3
         do j = 1, n2
            do i = 1, n1
               ii = i - 1
               arr(i, j, k) = arr_tmp(ii, j, k)
            end do
         end do
      end do
      !$acc end parallel loop
      !deallocate(arr_tmp)
    end select
    !
    return
  end subroutine prep_dctiif
  !
  subroutine posp_dctiif(n1, n2, n3,idir,arr,is_swap_order,is_negate_even)
    !
    ! post-processing of a signal to perform a fast forward discrete
    ! cosine transform with FFTs (see Makhoul 1980)
    !
    ! post-processing required for computing the corresponding
    ! discrete sine transform (DST) may also be performed
    ! if one of the last boolean input variables are .true.
    !
    implicit none
    !
    integer, intent(in) :: n1, n2, n3       ! dimensions of input/output array
    integer , intent(in   )                   :: idir
    real(rp), intent(inout), dimension(:,:,:) :: arr
    logical , intent(in   )                   :: is_swap_order  ! swap order of the elements of the input array? (for DST)
    logical , intent(in   )                   :: is_negate_even ! negate every other element of the input array?
    !
!!!    real(rp), allocatable, dimension(:,:,:) :: arr_tmp
!!!    real(rp), allocatable, dimension(:,:)   :: arr_sincos
    integer :: i, j, k, ii, nn
    !@cuf attributes(device) :: arr
    real(rp) :: arg, carg, sarg
    !
    !nn = n(idir) - 2
    select case (idir)
    case(1)
      nn = n1 - 2
      if(.not. allocated(arr_tmp)) allocate (arr_tmp(0:n1 - 1, n2, n3))
      !
!!!      if(.not. allocated(arr_sincos)) allocate (arr_sincos(2, 0:nn/2))
!!!      !$acc parallel loop collapse(1) present(arr_sincos)
!!!      do i = 1, nn + 2, 2
!!!         ii = (i - 1)/2
!!!         arg = -pi*ii/(2.*nn)
!!!         arr_sincos(1, ii) = sin(arg)
!!!         arr_sincos(2, ii) = cos(arg)
!!!       end do
      !
      !$acc parallel loop collapse(3) present(arr_tmp, arr)
      do k=1,n3
        do j=1,n2
           do i=1,nn+2,2 
             !
             ii = (i - 1)/2
             !!! arr_tmp(ii   ,j,k) =    real( &
             !!!                          2.*exp(-ri_unit*pi*ii/(2.*nn))*cmplx(arr(i,j,k),arr(i+1,j,k),rp) &
             !!!                         )
             !!! arr_tmp(nn-ii,j,k) = - aimag( &
             !!!                          2.*exp(-ri_unit*pi*ii/(2.*nn))*cmplx(arr(i,j,k),arr(i+1,j,k),rp) &
             !!!                         ) ! = 0 for ii=0
             arg = -pi*ii/(2.*nn)
             !!! call sincos(arg,sarg,carg)
             !!! carg = arr_sincos(2, ii)!cos(arg)
             !!! sarg = arr_sincos(1, ii)!sin(arg)
             carg = cos(arg)
             sarg = sin(arg)
             arr_tmp(ii, j, k) = 2.*(carg*arr(i, j, k) - sarg*arr(i + 1, j, k))
             arr_tmp(nn - ii, j, k) = -2.*(sarg*arr(i, j, k) + carg*arr(i + 1, j, k))
             !
            end do
         end do
      end do 
      !$acc end parallel loop
      !
      !$acc parallel loop collapse(3) present(arr_tmp, arr)
      do k=1,n3
        do j=1,n2
          do i=1,nn
            ii = i - 1
            arr(i, j, k) = arr_tmp(ii, j, k)
          end do
        end do
      end do
      !$acc end parallel loop
      !
      if(is_swap_order) then
        !
        !$acc parallel loop collapse(3) present(arr_tmp, arr)
        do k=1,n3
          do j=1,n2
            do i=1,nn
              ii = i - 1
              arr_tmp(ii, j, k) = arr(i, j, k) ! redundant
            end do
          end do
        end do
        !$acc end parallel loop
        !
        !$acc parallel loop collapse(3) present(arr_tmp, arr)
        do k = 1, n3
          do j = 1, n2
            do i = 1, nn
              ii = i - 1
              arr(i, j, k) = arr_tmp(nn - 1 - ii, j, k)
            end do
          end do
        end do
        !$acc end parallel loop
        !
      end if
      !
      if(is_negate_even) then
        !
        !$acc parallel loop collapse(3) present(arr)
        do k = 1, n3
          do j = 1, n2
            do i = 1, nn
              if(mod(i,2).eq.0) then
                arr(i, j, k) = -arr(i, j, k)
              end if
            end do
          end do
        end do
        !$acc end parallel loop
        !
      end if
    end select
    !
    return
  end subroutine posp_dctiif
  !
  subroutine prep_dctiib(n1, n2, n3, idir, arr, is_swap_order, is_negate_even)
    !
    ! pre-processing of a signal to perform a fast backward
    ! discrete cosine transform (DST) with FFTs (see Makhoul 1980)
    !
    ! pre-processing required for computing the corresponding
    ! discrete sine transform (DST) may also be performed
    ! if one of the last boolean input variables are .true.
    !
    implicit none
    integer, intent(in) :: n1, n2, n3       ! dimensions of input/output array
    integer, intent(in) :: idir
    real(rp), intent(inout), dimension(:, :, :) :: arr
    logical, intent(in) :: is_swap_order  ! swap order of the elements of the input array? (for DST)
    logical, intent(in) :: is_negate_even ! negate every other element of the input array?
!!!     real(rp), allocatable, dimension(:, :, :) :: arr_tmp
!!!     real(rp), allocatable, dimension(:, :) :: arr_sincos
    integer :: i, j, k, nn, ii
!@cuf attributes(device) :: arr
!!! !@cuf attributes(device) :: arr_tmp
!!! !@cuf attributes(device) :: arr_sincos
    real(rp) :: arg, carg, sarg
    !
    !nn = n(idir) - 2
    select case (idir)
    case (1)
       nn = n1 - 2
       if (.not. allocated(arr_tmp)) allocate (arr_tmp(0:n1 - 1, n2, n3))
!!       if (.not. allocated(arr_sincos)) allocate (arr_sincos(2, 0:nn/2))
       !
       if (is_swap_order) then
          !$acc parallel loop collapse(3) present(arr_tmp, arr)
          do k = 1, n3
             do j = 1, n2
                do i = 1, nn
                   ii = i - 1
                   arr_tmp(ii, j, k) = arr(i, j, k)
                end do
             end do
          end do
          !$acc parallel loop collapse(3) present(arr_tmp, arr)
          do k = 1, n3
             do j = 1, n2
                do i = 1, nn
                   ii = i - 1
                   arr(i, j, k) = arr_tmp(nn - 1 - ii, j, k)
                end do
             end do
          end do
       end if
       !
       if (is_negate_even) then
          !$acc parallel loop collapse(3) present(arr_tmp, arr)
          do k = 1, n3
             do j = 1, n2
                do i = 1, nn
                   if (mod(i, 2).eq.0) then
                      arr(i, j, k) = -arr(i, j, k)
                   end if
                end do
             end do
          end do
       end if
       !
       !$acc parallel loop collapse(2) present(arr)
       do k = 1, n3
          do j = 1, n2
             arr(nn + 1, j, k) = 0.
             arr(nn + 2, j, k) = 0.
          end do
       end do
       !
!!       !$acc parallel loop collapse(1) present(arr_sincos)
!!       do i = 1, nn + 2, 2
!!          ii = (i - 1)/2
!!          arg = pi*ii/(2.*nn)
!!          arr_sincos(1, ii) = sin(arg)
!!          arr_sincos(2, ii) = cos(arg)
!!       end do
       !
       !$acc parallel loop collapse(3) present(arr_tmp, arr)
       do k = 1, n3
          do j = 1, n2
             do i = 1, nn + 2, 2
                ii = (i - 1)/2
                !arr_tmp(2*ii  ,j,k)  = real( 1.*exp(ri_unit*pi*ii/(2.*nn))*(arr(ii+1,j,k)-ri_unit*arr(nn-ii+1,j,k)))
                !arr_tmp(2*ii+1,j,k)  = aimag(1.*exp(ri_unit*pi*ii/(2.*nn))*(arr(ii+1,j,k)-ri_unit*arr(nn-ii+1,j,k)))
                arg = pi*ii/(2.*nn)
                !call sincos(arg,sarg,carg)
                carg = cos(arg)
                sarg = sin(arg)
                arr_tmp(2*ii, j, k) = 1.*(carg*arr(ii + 1, j, k) + sarg*arr(nn - ii + 1, j, k))
                arr_tmp(2*ii + 1, j, k) = 1.*(sarg*arr(ii + 1, j, k) - carg*arr(nn - ii + 1, j, k))
             end do
          end do
       end do
       !
       !$acc parallel loop collapse(3) present(arr_tmp, arr)
       do k = 1, n3
          do j = 1, n2
             do i = 1, nn + 2
                ii = i - 1
                arr(i, j, k) = arr_tmp(ii, j, k)
             end do
          end do
       end do
       !deallocate(arr_tmp)
    end select
    !
    return
  end subroutine prep_dctiib
  !
  subroutine posp_dctiib(n1, n2, n3, idir, arr, is_swap_order, is_negate_even)
    !
    ! post-processing of a signal to perform a fast forward
    ! discrete cosine transform (DCT) with FFTs (see Makhoul 1980)
    !
    ! the input signal v(n) is post-processed into a signal x(n)
    ! as follows:
    !
    ! v(n) = x(2n       ),              0 <= n <= floor((N-1)/2)
    !      = x(2N -2n -1), floor((N+1)/2) <= n <= N-1
    ! with n = 0,...,N-1 and N is the total number of elements
    ! of the signal.
    !
    ! post-processing required for computing the corresponding
    ! discrete sine transform (DST) may also be performed
    ! if one of the last boolean input variables are .true.
    !
    implicit none
    integer, intent(in) :: n1, n2, n3       ! dimensions of input/output array
    integer, intent(in) :: idir                  ! array direction where the transform is taken
    real(rp), intent(inout), dimension(:, :, :) :: arr ! input/output array
    logical, intent(in) :: is_swap_order  ! swap order of the elements of the input array? (for DST)
    logical, intent(in) :: is_negate_even ! negate every other element of the input array?
!!!     real(rp), allocatable, dimension(:, :, :) :: arr_tmp
    integer :: i, j, k, nn, ii
!@cuf attributes(device) :: arr
!!! !@cuf attributes(device) :: arr_tmp
    !
    !nn = n(idir)
    select case (idir)
    case (1)
       nn = n1
       if (.not. allocated(arr_tmp)) allocate (arr_tmp(0:n1 - 1, n2, n3))
       !
       !$acc parallel loop collapse(3) present(arr_tmp, arr)
       do k = 1, n3
          do j = 1, n2
             do i = 1, n1
                ii = i - 1
                arr_tmp(ii, j, k) = arr(i, j, k)
             end do
          end do
       end do
       !$acc parallel loop collapse(3) present(arr_tmp, arr)
       do k = 1, n3
          do j = 1, n2
             do i = 1, n1
                ii = i - 1
                if (ii.le.(nn - 1)/2) then
                   arr(2*ii + 1, j, k) = arr_tmp(ii, j, k)
                else
                   arr(2*(nn - ii) - 1 + 1, j, k) = arr_tmp(ii, j, k)
                end if
             end do
          end do
       end do
       !
       if (is_swap_order) then
          !$acc parallel loop collapse(3) present(arr)
          do k = 1, n3
             do j = 1, n2
                do i = 1, nn
                   ii = i - 1
                   arr_tmp(ii, j, k) = arr(i, j, k)
                end do
             end do
          end do
          !
          !$acc parallel loop collapse(3) present(arr_tmp, arr)
          do k = 1, n3
             do j = 1, n2
                do i = 1, nn
                   ii = i - 1
                   arr(i, j, k) = arr_tmp(nn - 1 - ii, j, k)
                end do
             end do
          end do
       end if
       !
       if (is_negate_even) then
          !$acc parallel loop collapse(3) present(arr)
          do k = 1, n3
             do j = 1, n2
                do i = 1, nn
                   if (mod(i, 2).eq.0) then
                      arr(i, j, k) = -arr(i, j, k)
                   end if
                end do
             end do
          end do
       end if
       !deallocate(arr_tmp)
    end select
    !
    return
  end subroutine posp_dctiib
  !
  subroutine signal_processing(pre_or_pos,f_or_b,cbc,c_or_f,n,idir,arr)
    !
    ! wrapper subroutine for signal processing to compute FFT-based transforms
    ! (can also be done with pointers to a subroutine like in initgrid.f90)
    !
    implicit none
    !
    integer         , intent(in   )                   :: pre_or_pos ! prior (0) or after (1) fft
    character(len=1), intent(in   )                   :: f_or_b     ! forward or backward transform
    character(len=2), intent(in   )                   :: cbc        ! type of boundary condition
    character(len=1), intent(in   )                   :: c_or_f     ! cell- or face-centred BC?
    integer         , intent(in   ), dimension(3)     :: n
    integer         , intent(in   )                   :: idir
    real(rp)        , intent(inout), dimension(:,:,:) :: arr
    !
!@cuf attributes(device) :: arr
    integer :: n1, n2, n3
!@cuf integer :: istat
    !
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)
    !
    select case(cbc)
    case('PP')
      !
      select case(f_or_b)
      case('F')
        if(    pre_or_pos.eq.0) then
          return
        elseif(pre_or_pos.eq.1) then
          call posp_fftf(n1, n2, n3, idir, arr)
        else
        end if
      case('B')
        if(    pre_or_pos.eq.0) then
          call prep_fftb(n1, n2, n3, idir, arr)
        elseif(pre_or_pos.eq.1) then
          return
        else
        end if
      end select
      return
      !
    case('NN')
      !
      if(c_or_f.eq.'c') then
        select case (f_or_b)
        case('F')
          ! I AM HERE
          if(pre_or_pos.eq.0) then
            call prep_dctiif(n1, n2, n3, idir, arr,.false.,.false.)
          elseif (pre_or_pos.eq.1) then
            call posp_dctiif(n1, n2, n3, idir, arr,.false.,.false.)
          else
          end if
        case('B')
          if(pre_or_pos.eq.0) then
            call prep_dctiib(n1, n2, n3, idir, arr,.false.,.false.)
          elseif(pre_or_pos.eq.1) then
            call posp_dctiib(n1, n2, n3, idir, arr,.false.,.false.)
          else
          end if
        case default
        end select
      end if
      !
    case('DD')
      if(c_or_f.eq.'c') then
        select case (f_or_b)
        case('F')
          if(pre_or_pos.eq.0) then
            call prep_dctiif(n1, n2, n3, idir, arr,.false.,.true. )
          elseif (pre_or_pos.eq.1) then
            call posp_dctiif(n1, n2, n3, idir, arr,.true. ,.false.)
          else
          end if
        case('B')
          if(pre_or_pos.eq.0) then
            call prep_dctiib(n1, n2, n3, idir, arr,.true. ,.false.)
          elseif (pre_or_pos.eq.1) then
            call posp_dctiib(n1, n2, n3, idir, arr,.false.,.true. )
          else
          end if
        case default
        end select
      end if
    case default
      if(myid.eq.0) print*, "In GPU, only PP, DD, and NN boundary conditions for pressure are supported"
      if(myid.eq.0) print*, "Check dns.in"
      call MPI_FINALIZE(ierr)
      call exit
    end select
    !
    return
  end subroutine signal_processing
#endif
  !
end module mod_fft
