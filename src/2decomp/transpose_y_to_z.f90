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

! This file contains the routines that transpose data from Y to Z pencil

#if defined(_OPENACC)
  subroutine transpose_yt_to_zp(src, dst, yt, zp, normfft, opt_decomp)

    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT) :: src
    real(mytype), dimension(:,:,:), intent(INOUT) :: dst
    real(mytype), dimension(:,:,:), intent(INOUT) :: yt
    real(mytype), dimension(0:,0:,0:), intent(INOUT) :: zp
    real(mytype), intent(IN) :: normfft

    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    attributes( managed ) :: src, dst
    attributes( device  ) :: yt
    attributes( managed ) :: zp
    TYPE(DECOMP_INFO) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror, istat, m, i1, i2, pos
    integer :: iter, dest, sorc, pow2
    integer :: i,j,k,kstart,kend,jstart,jend,decomp_y2distm, idx
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
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,s3
    do i=1,s2
    do j=1,s1
      src(j,i,k) = yt(i,j,k)
    enddo
    enddo
    enddo
#endif

    if(IAND(dims(2),dims(2)-1)==0) then
      pow2 = 1
    else
      pow2 = 0
    endif

    ! rearrange source array as send buffer
    do iter=1,dims(2)-1
       if( pow2 ) then
         dest = IEOR(row_rank,iter)
       else
         dest = mod(row_rank + iter, dims(2))
       endif
       m = dest
       pos = decomp%y2disp(m) + 1
#if !defined(_EPHC)
       istat = cudaMemcpy2DAsync( work1_r_d(pos), s1*(decomp%y2dist(m)), src(1,decomp%y2idx(m),1), s1*s2, s1*(decomp%y2dist(m)), s3, stream=a2a_d2h )
#else
       jstart = decomp%y2idx(m)
       jend   = jstart + decomp%y2dist(m) - 1
       decomp_y2distm = decomp%y2dist(m)
       !$cuf kernel do(3) <<<*,*,stream=a2a_d2h>>>
       do k=1,s3
       do j=jstart, jend
       do i=1,s1
         idx = i + s1*(j-jstart) + s1*decomp_y2distm*(k-1)
         work1_r_d(pos + idx - 1) = yt( j, i, k ) 
       enddo
       enddo
       enddo
#endif
       istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
    end do

    ! self
    m = row_rank
#if !defined(_EPHC)
    istat = cudaMemcpy2DAsync( dst(1,1,decomp%z2idx(m)), s1*(decomp%y2dist(m)), src(1,decomp%y2idx(m),1), s1*s2, s1*(decomp%y2dist(m)), s3, stream=a2a_comp )
#else
    kstart = decomp%z2idx(m)
    jstart = decomp%y2idx(m)
    !$cuf kernel do(3) <<<*,(8,8,8),stream=a2a_comp>>>
    do k=1,s3
    do j=1,d2
    do i=1,s1
      zp(i,j,k+kstart-1) = yt(j+jstart-1,i,k)*normfft
    enddo
    enddo
    enddo
#endif

    do iter=1,dims(2)-1
      if( pow2 ) then
        sorc = IEOR(row_rank,iter)
      else
        sorc = mod(row_rank - iter + dims(2), dims(2))
      endif
      m = sorc
      call MPI_IRECV( work2_r_d(decomp%z2disp(m)+1), decomp%z2cnts(m), real_type, m, 0, DECOMP_2D_COMM_ROW, a2a_requests(iter), ierror)
    end do

    do iter=1,dims(2)-1
       if( pow2 ) then
          dest = IEOR(row_rank,iter)
          sorc = dest
       else
          dest = mod(row_rank + iter, dims(2))
          sorc = mod(row_rank - iter + dims(2), dims(2))
       endif
       m = dest
       istat = cudaEventSynchronize( a2a_event(iter) )
       call MPI_SEND( work1_r_d(decomp%y2disp(m)+1), decomp%y2cnts(m), real_type, m, 0, DECOMP_2D_COMM_ROW, ierror)
       call MPI_WAIT(a2a_requests(iter), MPI_STATUS_IGNORE, ierror)
       m = sorc
       pos = decomp%z2disp(m) + 1
#if !defined(_EPHC)
       istat = cudaMemcpyAsync( dst(1,1,decomp%z2idx(m)), work2_r_d(pos), decomp%z2cnts(m), a2a_h2d )
#else
       kstart = decomp%z2idx(m)
       kend = kstart + decomp%z2dist(m) - 1

       !$cuf kernel do(3) <<<*,*,stream=a2a_h2d>>>
       do k=kstart, kend
       do j=1,d2
       do i=1,d1
         idx = i + d1*(j-1) + d1*d2*(k-1)
         zp( i, j, k ) = work2_r_d( idx )*normfft
       enddo
       enddo
       enddo
#endif
    end do

#if !defined(_EPHC)
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,d3
    do j=1,d2
    do i=1,d1
      zp(i,j,k) = dst(i,j,k)*normfft
    enddo
    enddo
    enddo
#endif

    istat = cudaEventRecord( a2a_event(0), 0 )
    istat = cudaEventSynchronize( a2a_event(0) )
#else
    ! rearrange source array as send buffer
    do m=0,dims(2)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%y2dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%y2dist(m)-1
       end if
       pos = decomp%y2disp(m) + 1
       istat = cudaMemcpy2D( work1_r(pos), s1*(i2-i1+1), src(1,i1,1), s1*s2, s1*(i2-i1+1), s3, cudaMemcpyDeviceToHost )
    end do

    call MPI_ALLTOALLV(work1_r, decomp%y2cnts, decomp%y2disp, &
         real_type, work2_r, decomp%z2cnts, decomp%z2disp, &
         real_type, DECOMP_2D_COMM_ROW, ierror)

    istat = cudaMemcpy( dst, work2_r, d1*d2*d3, cudaMemcpyHostToDevice )
#endif


    return
  end subroutine transpose_yt_to_zp

  subroutine transpose_y_to_z_real_d(src, dst, opt_decomp)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    attributes( managed ) :: src, dst
    TYPE(DECOMP_INFO) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror, istat, m, i1, i2, pos
    integer :: iter, dest, sorc, pow2
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
    if(IAND(dims(2),dims(2)-1)==0) then
      pow2 = 1
    else
      pow2 = 0
    endif

    ! rearrange source array as send buffer
    do iter=1,dims(2)-1
       if( pow2 ) then
         dest = IEOR(row_rank,iter)
       else
         dest = mod(row_rank + iter, dims(2))
       endif
       m = dest
       pos = decomp%y2disp(m) + 1
       istat = cudaMemcpy2DAsync( work1_r_d(pos), s1*(decomp%y2dist(m)), src(1,decomp%y2idx(m),1), s1*s2, s1*(decomp%y2dist(m)), s3, stream=a2a_d2h )
       istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
    end do

    ! self
    m = row_rank
    istat = cudaMemcpy2DAsync( dst(1,1,decomp%z2idx(m)), s1*(decomp%y2dist(m)), src(1,decomp%y2idx(m),1), s1*s2, s1*(decomp%y2dist(m)), s3, stream=a2a_comp )

    do iter=1,dims(2)-1
      if( pow2 ) then
        sorc = IEOR(row_rank,iter)
      else
        sorc = mod(row_rank - iter + dims(2), dims(2))
      endif
      m = sorc
      call MPI_IRECV( work2_r_d(decomp%z2disp(m)+1), decomp%z2cnts(m), real_type, m, 0, DECOMP_2D_COMM_ROW, a2a_requests(iter), ierror)
    end do

    do iter=1,dims(2)-1
       if( pow2 ) then
          dest = IEOR(row_rank,iter)
          sorc = dest
       else
          dest = mod(row_rank + iter, dims(2))
          sorc = mod(row_rank - iter + dims(2), dims(2))
       endif
       m = dest
       istat = cudaEventSynchronize( a2a_event(iter) )
       call MPI_SEND( work1_r_d(decomp%y2disp(m)+1), decomp%y2cnts(m), real_type, m, 0, DECOMP_2D_COMM_ROW, ierror)
       call MPI_WAIT(a2a_requests(iter), MPI_STATUS_IGNORE, ierror)
       m = sorc
       pos = decomp%z2disp(m) + 1
       istat = cudaMemcpyAsync( dst(1,1,decomp%z2idx(m)), work2_r_d(pos), decomp%z2cnts(m), a2a_h2d )
    end do
    istat = cudaEventRecord( a2a_event(0), 0 )
    istat = cudaEventSynchronize( a2a_event(0) )
#else
    ! rearrange source array as send buffer
    do m=0,dims(2)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%y2dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%y2dist(m)-1
       end if
       pos = decomp%y2disp(m) + 1
       istat = cudaMemcpy2D( work1_r(pos), s1*(i2-i1+1), src(1,i1,1), s1*s2, s1*(i2-i1+1), s3, cudaMemcpyDeviceToHost )
    end do

    call MPI_ALLTOALLV(work1_r, decomp%y2cnts, decomp%y2disp, &
         real_type, work2_r, decomp%z2cnts, decomp%z2disp, &
         real_type, DECOMP_2D_COMM_ROW, ierror)

    istat = cudaMemcpy( dst, work2_r, d1*d2*d3, cudaMemcpyHostToDevice )
#endif


    return
  end subroutine transpose_y_to_z_real_d
#endif

  subroutine transpose_y_to_z_real(src, dst, opt_decomp)

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
    work1_p = decomp%ROW_INFO%SND_P
    call mem_split_yz_real(src, s1, s2, s3, work1, dims(2), &
         decomp%y2dist, decomp)
#else
    call mem_split_yz_real(src, s1, s2, s3, work1_r, dims(2), &
         decomp%y2dist, decomp)
#endif

    ! define receive buffer
#if defined(_SHM)
    work2_p = decomp%ROW_INFO%RCV_P
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
#if defined(_SHM)
    if (decomp%ROW_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%y2cnts_s, decomp%y2disp_s, &
            real_type, work2, decomp%z2cnts_s, decomp%z2disp_s, &
            real_type, decomp%ROW_INFO%SMP_COMM, ierror)
    end if
#else
#if defined(_EVEN)
    if (decomp%even) then
       call MPI_ALLTOALL(work1_r, decomp%y2count, &
            real_type, dst, decomp%z2count, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
    else
       call MPI_ALLTOALL(work1_r, decomp%y2count, &
            real_type, work2_r, decomp%z2count, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
    end if
#else
    call MPI_ALLTOALLV(work1_r, decomp%y2cnts, decomp%y2disp, &
         real_type, dst, decomp%z2cnts, decomp%z2disp, &
         real_type, DECOMP_2D_COMM_ROW, ierror)
#endif
#endif

    ! rearrange receive buffer
#if defined(_SHM)
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
    call mem_merge_yz_real(work2, d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
#else
#if defined(_EVEN)
    if (.not. decomp%even) then
       call mem_merge_yz_real(work2_r, d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
    end if
#else
    ! note the receive buffer is already in natural (i,j,k) order
    ! so no merge operation needed
#endif
#endif
    
    return
  end subroutine transpose_y_to_z_real


#if defined(_OCC)
  subroutine transpose_y_to_z_real_start(handle, src, dst, sbuf, rbuf, &
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
    call mem_split_yz_real(src, s1, s2, s3, sbuf, dims(2), &
         decomp%y2dist, decomp)

#if defined(_EVEN)
    call NBC_IALLTOALL(sbuf, decomp%y2count, real_type, &
         rbuf, decomp%z2count, real_type, &
         DECOMP_2D_COMM_ROW, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%y2cnts, decomp%y2disp, real_type, &
         rbuf, decomp%z2cnts, decomp%z2disp, real_type, &
         DECOMP_2D_COMM_ROW, handle, ierror)
#endif

    return
  end subroutine transpose_y_to_z_real_start


  subroutine transpose_y_to_z_real_wait(handle, src, dst, sbuf, rbuf, &
       opt_decomp)

    implicit none
    
    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    call NBC_WAIT(handle, ierror)

    dst = rbuf

    return
  end subroutine transpose_y_to_z_real_wait
#endif


  subroutine transpose_y_to_z_complex(src, dst, opt_decomp)

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
    work1_p = decomp%ROW_INFO%SND_P_c
    call mem_split_yz_complex(src, s1, s2, s3, work1, dims(2), &
         decomp%y2dist, decomp)
#else
    call mem_split_yz_complex(src, s1, s2, s3, work1_c, dims(2), &
         decomp%y2dist, decomp)
#endif
    
    ! define receive buffer
#if defined(_SHM)
    work2_p = decomp%ROW_INFO%RCV_P_c
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif
    
#if defined(_SHM)
    if (decomp%ROW_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%y2cnts_s, decomp%y2disp_s, &
            complex_type, work2, decomp%z2cnts_s, decomp%z2disp_s, &
            complex_type, decomp%ROW_INFO%SMP_COMM, ierror)
    end if
#else
#if defined(_EVEN)
    if (decomp%even) then
       call MPI_ALLTOALL(work1_c, decomp%y2count, &
            complex_type, dst, decomp%z2count, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
    else
       call MPI_ALLTOALL(work1_c, decomp%y2count, &
            complex_type, work2_c, decomp%z2count, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
    end if
#else
    call MPI_ALLTOALLV(work1_c, decomp%y2cnts, decomp%y2disp, &
         complex_type, dst, decomp%z2cnts, decomp%z2disp, &
         complex_type, DECOMP_2D_COMM_ROW, ierror)
#endif
#endif

    ! rearrange receive buffer
#if defined(_SHM)
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
    call mem_merge_yz_complex(work2, d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
#else
#if defined(_EVEN)
    if (.not. decomp%even) then
       call mem_merge_yz_complex(work2_c, d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
    end if
#else
    ! note the receive buffer is already in natural (i,j,k) order
    ! so no merge operation needed
#endif
#endif

    return
  end subroutine transpose_y_to_z_complex


#if defined(_OCC)
  subroutine transpose_y_to_z_complex_start(handle, src, dst, sbuf, &
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
    call mem_split_yz_complex(src, s1, s2, s3, sbuf, dims(2), &
         decomp%y2dist, decomp)

#if defined(_EVEN)
    call NBC_IALLTOALL(sbuf, decomp%y2count, &
         complex_type, rbuf, decomp%z2count, &
         complex_type, DECOMP_2D_COMM_ROW, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%y2cnts, decomp%y2disp, &
         complex_type, rbuf,decomp%z2cnts, decomp%z2disp, &
         complex_type, DECOMP_2D_COMM_ROW, handle, ierror)
#endif

    return
  end subroutine transpose_y_to_z_complex_start


  subroutine transpose_y_to_z_complex_wait(handle, src, dst, sbuf, &
       rbuf, opt_decomp)

    implicit none
    
    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    call NBC_WAIT(handle, ierror)

    dst = rbuf

    return
  end subroutine transpose_y_to_z_complex_wait
#endif


  ! pack/unpack ALLTOALL(V) buffers

  subroutine mem_split_yz_real(in,n1,n2,n3,out,iproc,dist,decomp)

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
       pos = decomp%y2disp_o(m) + 1
#else
#if defined(_EVEN)
       pos = m * decomp%y2count + 1
#else
       pos = decomp%y2disp(m) + 1
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
  end subroutine mem_split_yz_real


  subroutine mem_split_yz_complex(in,n1,n2,n3,out,iproc,dist,decomp)

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
       pos = decomp%y2disp_o(m) + 1
#else
#if defined(_EVEN)
       pos = m * decomp%y2count + 1
#else
       pos = decomp%y2disp(m) + 1
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
  end subroutine mem_split_yz_complex


  subroutine mem_merge_yz_real(in,n1,n2,n3,out,iproc,dist,decomp)

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
       pos = decomp%z2disp_o(m) + 1
#else
#if defined(_EVEN)
       pos = m * decomp%z2count + 1
#else
       pos = decomp%z2disp(m) + 1
#endif
#endif

       do k=i1,i2
          do j=1,n2
             do i=1,n1
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_yz_real


  subroutine mem_merge_yz_complex(in,n1,n2,n3,out,iproc,dist,decomp)
    
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
       pos = decomp%z2disp_o(m) + 1
#else
#if defined(_EVEN)
       pos = m * decomp%z2count + 1
#else
       pos = decomp%z2disp(m) + 1
#endif
#endif

       do k=i1,i2
          do j=1,n2
             do i=1,n1
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_yz_complex
