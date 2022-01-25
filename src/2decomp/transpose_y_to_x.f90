!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contains the routines that transpose data from Y to X pencil

#if defined(_OPENACC)
  subroutine transpose_yct_to_x(src, dst, yct, opt_decomp)
    implicit none
    real(mytype), dimension(:,:,:), intent(INOUT) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    real(mytype), dimension(:,:,:), intent(INOUT) :: yct
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    attributes( managed ) :: src, dst
    attributes( device ) :: yct
    TYPE(DECOMP_INFO) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror, istat, m, i1, i2, pos
    integer :: iter, dest, sorc, pow2
    integer :: i,j,k,jstart,jend,idx,decomp_y1distm
    integer ::       istart,iend,    decomp_x1distm
#if defined(_USE_NVTX_FFT)
    call profiler_start("tranYctX", tag = .true., tag_color = COLOR_WHITE)
#endif
    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

#if defined(_EPA2A)

#if !defined(_EPHC)
    !$cuf kernel do(3) <<<*,(8,8,8)>>>
    do k=1,s3
    do j=1,s2
    do i=1,s1
      src(i,j,k) = yct(j,i,k)
    end do
    end do
    end do
#endif

    if(IAND(dims(1),dims(1)-1)==0) then
      pow2 = 1
    else
      pow2 = 0
    endif

    ! rearrange source array as send buffer
    do iter=1,dims(1)-1
       if( pow2 ) then
         dest = IEOR(col_rank,iter)
       else
         dest = mod(col_rank + iter, dims(1))
       endif
       m = dest
       pos = decomp%y1disp(m) + 1

#if !defined(_EPHC)
       istat = cudaMemcpy2DAsync( work1_r_d(pos), s1*(decomp%y1dist(m)), src(1,decomp%y1idx(m),1), s1*s2, s1*(decomp%y1dist(m)),s3, stream=a2a_d2h )
#else
       jstart = decomp%y1idx(m)
       jend   = jstart + decomp%y1dist(m) - 1
       decomp_y1distm = decomp%y1dist(m)
       !$cuf kernel do(3) <<<*,(8,8,8),stream=a2a_d2h>>>
       do k=1,s3
       do j=jstart, jend
       do i=1,s1
         idx = i + s1*(j-jstart) + s1*decomp_y1distm*(k-1)
         work1_r_d(pos + idx - 1) = yct(j,i,k)
       enddo
       enddo
       enddo
#endif
       istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
    end do

    ! self
    m = col_rank
    pos = decomp%x1disp(m) + 1
    !if( dims(1) .eq. 1 ) then
    !  istat = cudaMemcpy2DAsync( dst, s1*s2, src, s1*s2, s1*s2, s3, stream=a2a_comp )
    !else
#if !defined(_EPHC)
      !TODO: replace these two copy with a 3D copy or custom kernel for direct src => dst
      istat = cudaMemcpy2DAsync( work2_r_d(pos), s1*(decomp%y1dist(m)), src(1,decomp%y1idx(m),1), s1*s2, s1*(decomp%y1dist(m)),s3, stream=a2a_comp )
      istat = cudaMemcpy2DAsync( dst(decomp%x1idx(m),1,1), d1, work2_r_d(pos), decomp%x1dist(m), decomp%x1dist(m), d2*d3, stream=a2a_comp )
#else
      istart = decomp%x1idx(m)
      iend   = istart + decomp%x1dist(m) - 1
      decomp_x1distm = decomp%x1dist(m)

      jstart = decomp%y1idx(m)
      jend   = jstart + decomp%y1dist(m) - 1
      decomp_y1distm = decomp%y1dist(m)
      !$cuf kernel do(3) <<<*,(8,8,8),stream=a2a_comp>>>
      do k=1,s3
      do j=jstart, jend
      do i=istart, iend
        dst(i,(j-jstart+1),k) = yct(j,(i-istart+1),k)
      enddo
      enddo
      enddo
#endif
    !endif

    do iter=1,dims(1)-1
      if( pow2 ) then
        sorc = IEOR(col_rank,iter)
      else
        sorc = mod(col_rank - iter + dims(1), dims(1))
      endif
      m = sorc
      call MPI_IRECV( work2_r_d(decomp%x1disp(m)+1), decomp%x1cnts(m), real_type, m, 0, DECOMP_2D_COMM_COL, a2a_requests(iter),ierror)
    end do

    do iter=1,dims(1)-1
       if( pow2 ) then
          dest = IEOR(col_rank,iter)
          sorc = dest
       else
          dest = mod(col_rank + iter, dims(1))
          sorc = mod(col_rank - iter + dims(1), dims(1))
       endif
       m = dest
       istat = cudaEventSynchronize( a2a_event(iter) )
       call MPI_SEND( work1_r_d(decomp%y1disp(m)+1), decomp%y1cnts(m), real_type, m, 0, DECOMP_2D_COMM_COL, ierror)
       call MPI_WAIT(a2a_requests(iter), MPI_STATUS_IGNORE, ierror)
       m = sorc
       pos = decomp%x1disp(m) + 1
       istat = cudaMemcpy2DAsync( dst(decomp%x1idx(m),1,1), d1, work2_r_d(pos), decomp%x1dist(m), decomp%x1dist(m), d2*d3,stream=a2a_comp )
    end do
    istat = cudaEventRecord( a2a_event(0), 0 )
    istat = cudaEventSynchronize( a2a_event(0) )
#else
    ! NOT SUPPORTED YET
    STOP

    ! rearrange source array as send buffer
    do m=0,dims(1)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%y1dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%y1dist(m)-1
       end if
       pos = decomp%y1disp(m) + 1
       istat = cudaMemcpy2D( work1_r(pos), s1*(i2-i1+1), src(1,i1,1), s1*s2, s1*(i2-i1+1), s3, cudaMemcpyDeviceToHost )
    end do

    call MPI_ALLTOALLV(work1_r, decomp%y1cnts, decomp%y1disp, &
         real_type, work2_r, decomp%x1cnts, decomp%x1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)

    ! rearrange receive buffer
    !call mem_merge_yx_real(work2_r, d1, d2, d3, dst, dims(1), &
    !     decomp%x1dist, decomp)
    do m=0,dims(1)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%x1dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%x1dist(m)-1
       end if
       pos = decomp%x1disp(m) + 1
       istat = cudaMemcpy2D( dst(i1,1,1), d1, work2_r(pos), i2-i1+1, i2-i1+1, d2*d3, cudaMemcpyHostToDevice )
    end do
#endif

#if defined(_USE_NVTX_FFT)
    call profiler_stop("tranYctX")
#endif

    return
  end subroutine transpose_yct_to_x

  subroutine transpose_y_to_x_real_d(src, dst, opt_decomp)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    attributes( managed ) :: src, dst
    TYPE(DECOMP_INFO) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror, istat, m, i1, i2, pos
    integer :: iter, dest, sorc, pow2
#if defined(_USE_NVTX_FFT)
    call profiler_start("tranYX", tag = .true., tag_color = COLOR_WHITE)
#endif
    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

#if defined(_EPA2A)
    if(IAND(dims(1),dims(1)-1)==0) then
      pow2 = 1
    else
      pow2 = 0
    endif

    ! rearrange source array as send buffer
    do iter=1,dims(1)-1
       if( pow2 ) then
         dest = IEOR(col_rank,iter)
       else
         dest = mod(col_rank + iter, dims(1))
       endif
       m = dest
       pos = decomp%y1disp(m) + 1
       istat = cudaMemcpy2DAsync( work1_r_d(pos), s1*(decomp%y1dist(m)), src(1,decomp%y1idx(m),1), s1*s2, s1*(decomp%y1dist(m)),s3, stream=a2a_d2h )
       istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
    end do

    ! self
    m = col_rank
    pos = decomp%x1disp(m) + 1
    if( dims(1) .eq. 1 ) then
      istat = cudaMemcpy2DAsync( dst, s1*s2, src, s1*s2, s1*s2, s3, stream=a2a_comp )
    else
      !TODO: replace these two copy with a 3D copy or custom kernel for direct src => dst
      istat = cudaMemcpy2DAsync( work2_r_d(pos), s1*(decomp%y1dist(m)), src(1,decomp%y1idx(m),1), s1*s2, s1*(decomp%y1dist(m)),s3, stream=a2a_comp )
      istat = cudaMemcpy2DAsync( dst(decomp%x1idx(m),1,1), d1, work2_r_d(pos), decomp%x1dist(m), decomp%x1dist(m), d2*d3, stream=a2a_comp )
    endif

    do iter=1,dims(1)-1
      if( pow2 ) then
        sorc = IEOR(col_rank,iter)
      else
        sorc = mod(col_rank - iter + dims(1), dims(1))
      endif
      m = sorc
      call MPI_IRECV( work2_r_d(decomp%x1disp(m)+1), decomp%x1cnts(m), real_type, m, 0, DECOMP_2D_COMM_COL, a2a_requests(iter),ierror)
    end do

    do iter=1,dims(1)-1
       if( pow2 ) then
          dest = IEOR(col_rank,iter)
          sorc = dest
       else
          dest = mod(col_rank + iter, dims(1))
          sorc = mod(col_rank - iter + dims(1), dims(1))
       endif
       m = dest
       istat = cudaEventSynchronize( a2a_event(iter) )
       call MPI_SEND( work1_r_d(decomp%y1disp(m)+1), decomp%y1cnts(m), real_type, m, 0, DECOMP_2D_COMM_COL, ierror)
       call MPI_WAIT(a2a_requests(iter), MPI_STATUS_IGNORE, ierror)
       m = sorc
       pos = decomp%x1disp(m) + 1
       istat = cudaMemcpy2DAsync( dst(decomp%x1idx(m),1,1), d1, work2_r_d(pos), decomp%x1dist(m), decomp%x1dist(m), d2*d3,stream=a2a_comp )
    end do
    istat = cudaEventRecord( a2a_event(0), 0 )
    istat = cudaEventSynchronize( a2a_event(0) )
#else
    ! rearrange source array as send buffer
    do m=0,dims(1)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%y1dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%y1dist(m)-1
       end if
       pos = decomp%y1disp(m) + 1
       istat = cudaMemcpy2D( work1_r(pos), s1*(i2-i1+1), src(1,i1,1), s1*s2, s1*(i2-i1+1), s3, cudaMemcpyDeviceToHost )
    end do

    call MPI_ALLTOALLV(work1_r, decomp%y1cnts, decomp%y1disp, &
         real_type, work2_r, decomp%x1cnts, decomp%x1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)

    ! rearrange receive buffer
    !call mem_merge_yx_real(work2_r, d1, d2, d3, dst, dims(1), &
    !     decomp%x1dist, decomp)
    do m=0,dims(1)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%x1dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%x1dist(m)-1
       end if
       pos = decomp%x1disp(m) + 1
       istat = cudaMemcpy2D( dst(i1,1,1), d1, work2_r(pos), i2-i1+1, i2-i1+1, d2*d3, cudaMemcpyHostToDevice )
    end do
#endif

#if defined(_USE_NVTX_FFT)
    call profiler_stop("tranYX")
#endif

    return
  end subroutine transpose_y_to_x_real_d
#endif

  subroutine transpose_y_to_x_real(src, dst, opt_decomp)

    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

#if defined(_SHM)
    real(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    ! rearrange source array as send buffer
#if defined(_SHM)
    work1_p = decomp%COL_INFO%SND_P
    call mem_split_yx_real(src, s1, s2, s3, work1, dims(1), &
         decomp%y1dist, decomp)
#else
    call mem_split_yx_real(src, s1, s2, s3, work1_r, dims(1), &
         decomp%y1dist, decomp)
#endif

    ! define receive buffer
#if defined(_SHM)
    work2_p = decomp%COL_INFO%RCV_P
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
    ! transpose using MPI_ALLTOALL(V)
#if defined(_SHM)
    if (decomp%COL_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%y1cnts_s, decomp%y1disp_s, &
            real_type, work2, decomp%x1cnts_s, decomp%x1disp_s, &
            real_type, decomp%COL_INFO%SMP_COMM, ierror)
    end if
#else
#if defined(_EVEN)
    call MPI_ALLTOALL(work1_r, decomp%y1count, &
         real_type, work2_r, decomp%x1count, &
         real_type, DECOMP_2D_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1_r, decomp%y1cnts, decomp%y1disp, &
         real_type, work2_r, decomp%x1cnts, decomp%x1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)
#endif
#endif

    ! rearrange receive buffer
#if defined(_SHM)
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
    call mem_merge_yx_real(work2, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#else
    call mem_merge_yx_real(work2_r, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#endif
    
    return
  end subroutine transpose_y_to_x_real


#if defined(_OCC)
  subroutine transpose_y_to_x_real_start(handle, src, dst, sbuf, rbuf, &
       opt_decomp)

    implicit none
    
    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: s1,s2,s3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    
    ! rearrange source array as send buffer
    call mem_split_yx_real(src, s1, s2, s3, sbuf, dims(1), &
         decomp%y1dist, decomp)

#if defined(_EVEN)
    call NBC_IALLTOALL(sbuf, decomp%y1count, real_type, &
         rbuf, decomp%x1count, real_type, &
         DECOMP_2D_COMM_COL, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%y1cnts, decomp%y1disp, real_type, &
         rbuf, decomp%x1cnts, decomp%x1disp, real_type, &
         DECOMP_2D_COMM_COL, handle, ierror)
#endif

    return
  end subroutine transpose_y_to_x_real_start


  subroutine transpose_y_to_x_real_wait(handle, src, dst, sbuf, rbuf, &
       opt_decomp)

    implicit none
    
    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    call NBC_WAIT(handle, ierror)

    ! rearrange receive buffer
    call mem_merge_yx_real(rbuf, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)

    return
  end subroutine transpose_y_to_x_real_wait
#endif


  subroutine transpose_y_to_x_complex(src, dst, opt_decomp)

    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

#if defined(_SHM)
    complex(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! rearrange source array as send buffer
#if defined(_SHM)
    work1_p = decomp%COL_INFO%SND_P_c
    call mem_split_yx_complex(src, s1, s2, s3, work1, dims(1), &
         decomp%y1dist, decomp)
#else
    call mem_split_yx_complex(src, s1, s2, s3, work1_c, dims(1), &
         decomp%y1dist, decomp)
#endif
    
    ! define receive buffer
#if defined(_SHM)
    work2_p = decomp%COL_INFO%RCV_P_c
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
    ! transpose using MPI_ALLTOALL(V)
#if defined(_SHM)
    if (decomp%COL_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%y1cnts_s, decomp%y1disp_s, &
            complex_type, work2, decomp%x1cnts_s, decomp%x1disp_s, &
            complex_type, decomp%COL_INFO%SMP_COMM, ierror)
    end if
#else
#if defined(_EVEN)
    call MPI_ALLTOALL(work1_c, decomp%y1count, &
         complex_type, work2_c, decomp%x1count, &
         complex_type, DECOMP_2D_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1_c, decomp%y1cnts, decomp%y1disp, &
         complex_type, work2_c, decomp%x1cnts, decomp%x1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)
#endif
#endif

    ! rearrange receive buffer
#if defined(_SHM)
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
    call mem_merge_yx_complex(work2, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#else
    call mem_merge_yx_complex(work2_c, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#endif

    return
  end subroutine transpose_y_to_x_complex


#if defined(_OCC)
  subroutine transpose_y_to_x_complex_start(handle, src, dst, sbuf, &
       rbuf, opt_decomp)

    implicit none
    
    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: s1,s2,s3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    
    ! rearrange source array as send buffer
    call mem_split_yx_complex(src, s1, s2, s3, sbuf, dims(1), &
         decomp%y1dist, decomp)

#if defined(_EVEN)
    call NBC_IALLTOALL(sbuf, decomp%y1count, &
         complex_type, rbuf, decomp%x1count, &
         complex_type, DECOMP_2D_COMM_COL, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%y1cnts, decomp%y1disp, &
         complex_type, rbuf, decomp%x1cnts, decomp%x1disp, &
         complex_type, DECOMP_2D_COMM_COL, handle, ierror)
#endif

    return
  end subroutine transpose_y_to_x_complex_start


  subroutine transpose_y_to_x_complex_wait(handle, src, dst, sbuf, &
       rbuf, opt_decomp)

    implicit none
    
    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    call NBC_WAIT(handle, ierror)

    ! rearrange receive buffer
    call mem_merge_yx_complex(rbuf, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)

    return
  end subroutine transpose_y_to_x_complex_wait
#endif


  ! pack/unpack ALLTOALL(V) buffers

  subroutine mem_split_yx_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(n1,n2,n3), intent(IN) :: in
    real(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#if defined(_SHM)
       pos = decomp%y1disp_o(m) + 1
#else
#if defined(_EVEN)
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_yx_real


  subroutine mem_split_yx_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(n1,n2,n3), intent(IN) :: in
    complex(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#if defined(_SHM)
       pos = decomp%y1disp_o(m) + 1
#else
#if defined(_EVEN)
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_yx_complex


  subroutine mem_merge_yx_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(*), intent(IN) :: in
    real(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#if defined(_SHM)
       pos = decomp%x1disp_o(m) + 1
#else
#if defined(_EVEN)
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_yx_real


  subroutine mem_merge_yx_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(*), intent(IN) :: in
    complex(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#if defined(_SHM)
       pos = decomp%x1disp_o(m) + 1
#else
#if defined(_EVEN)
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_yx_complex
