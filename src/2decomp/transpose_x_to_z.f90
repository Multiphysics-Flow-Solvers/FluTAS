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

! This file contains the routines that transpose data from X to Z pencil

#if defined(_OPENACC)
  subroutine transpose_xc_to_z(src, dst, pxc, opt_decomp)

    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    real(mytype), dimension(:,:,:), intent(INOUT) :: pxc
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    attributes( managed ) :: src, dst
    attributes( device ) :: pxc
    TYPE(DECOMP_INFO) :: decomp

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror, istat, m, i1, i2, pos
    integer :: iter, dest, sorc, pow2
    integer :: i,j,k,idx,istart,iend,decomp_x2distm
    integer ::           kstart,kend
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
    do j=1,s2
    do i=1,s1
      src(i,j,k) = pxc(i,j,k)
    end do
    end do
    end do
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
       pos = decomp%x2disp(m) + 1

#if !defined(_EPHC)
       istat = cudaMemcpy2DAsync( work1_r_d(pos), decomp%x2dist(m), src(decomp%x2idx(m),1,1), s1, decomp%x2dist(m), s2*s3, stream=a2a_d2h )
#else
       istart = decomp%x2idx(m)
       iend   = istart + decomp%x2dist(m) - 1
       decomp_x2distm = decomp%x2dist(m)
       !$cuf kernel do(3) <<<*,*,stream=a2a_d2h>>>
       do k=1,s3
       do j=1,s2
       do i=istart, iend
         idx = (i-istart+1) + decomp_x2distm*(j-1) + decomp_x2distm*s2*(k-1)
         work1_r_d(pos + idx - 1) = pxc(i,j,k)
       enddo
       enddo
       enddo
#endif
       istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
    end do

    ! self
    m = row_rank
#if !defined(_EPHC)
    istat = cudaMemcpy2DAsync( dst(1,1,decomp%z2idx(m)), d1, src(decomp%x2idx(m),1,1), s1, decomp%x2dist(m), s2*s3, stream=a2a_comp )
#else
      istart = decomp%x2idx(m)
      iend   = istart + decomp%x2dist(m) - 1
      decomp_x2distm = decomp%x2dist(m)

      kstart = decomp%z2idx(m)
      kend   = kstart + decomp%z2dist(m) - 1

      !$cuf kernel do(3) <<<*,*,stream=a2a_comp>>>
      do k=kstart, kend
      do j=1,s2
      do i=istart, iend
        dst((i-istart+1),j,k) = pxc(i,j,(k-kstart+1))
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
      call MPI_IRECV( work2_r_d(decomp%w2disp(m)+1), decomp%w2cnts(m), real_type, m, 0, DECOMP_2D_COMM_ROW, a2a_requests(iter), ierror)
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
       call MPI_SEND( work1_r_d(decomp%x2disp(m)+1), decomp%x2cnts(m), real_type, m, 0, DECOMP_2D_COMM_ROW, ierror)
       call MPI_WAIT(a2a_requests(iter), MPI_STATUS_IGNORE, ierror)
       m = sorc
       pos = decomp%w2disp(m) + 1
       istat = cudaMemcpyAsync( dst(1,1,decomp%z2idx(m)), work2_r_d(pos), decomp%w2cnts(m), a2a_h2d )
    end do
    istat = cudaEventRecord( a2a_event(0), 0 )
    istat = cudaEventSynchronize( a2a_event(0) )
#else
    ! rearrange source array as send buffer
    do m=0,dims(2)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%x2dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%x2dist(m)-1
       end if
       pos = decomp%x2disp(m) + 1
       istat = cudaMemcpy2D( work1_r(pos), i2-i1+1, src(i1,1,1), s1, i2-i1+1,s2*s3, cudaMemcpyDevicetoHost )
    end do

    call MPI_ALLTOALLV(work1_r, decomp%x2cnts, decomp%x2disp, &
         real_type, work2_r, decomp%w2cnts, decomp%w2disp, &
         real_type, DECOMP_2D_COMM_ROW, ierror)

    istat = cudaMemcpy( dst, work2_r, d1*d2*d3, cudaMemcpyHostToDevice )
#endif


    return
  end subroutine transpose_xc_to_z

  subroutine transpose_x_to_z_real_d(src, dst, opt_decomp)

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
       pos = decomp%x2disp(m) + 1
       istat = cudaMemcpy2DAsync( work1_r_d(pos), decomp%x2dist(m), src(decomp%x2idx(m),1,1), s1, decomp%x2dist(m), s2*s3, stream=a2a_d2h )
       istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
    end do

    ! self
    m = row_rank
    istat = cudaMemcpy2DAsync( dst(1,1,decomp%z2idx(m)), d1, src(decomp%x2idx(m),1,1), s1, decomp%x2dist(m), s2*s3, stream=a2a_comp )

    do iter=1,dims(2)-1
      if( pow2 ) then
        sorc = IEOR(row_rank,iter)
      else
        sorc = mod(row_rank - iter + dims(2), dims(2))
      endif
      m = sorc
      call MPI_IRECV( work2_r_d(decomp%w2disp(m)+1), decomp%w2cnts(m), real_type, m, 0, DECOMP_2D_COMM_ROW, a2a_requests(iter), ierror)
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
       call MPI_SEND( work1_r_d(decomp%x2disp(m)+1), decomp%x2cnts(m), real_type, m, 0, DECOMP_2D_COMM_ROW, ierror)
       call MPI_WAIT(a2a_requests(iter), MPI_STATUS_IGNORE, ierror)
       m = sorc
       pos = decomp%w2disp(m) + 1
       istat = cudaMemcpyAsync( dst(1,1,decomp%z2idx(m)), work2_r_d(pos), decomp%w2cnts(m), a2a_h2d )
    end do
    istat = cudaEventRecord( a2a_event(0), 0 )
    istat = cudaEventSynchronize( a2a_event(0) )
#else
    ! rearrange source array as send buffer
    do m=0,dims(2)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%x2dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%x2dist(m)-1
       end if
       pos = decomp%x2disp(m) + 1
       istat = cudaMemcpy2D( work1_r(pos), i2-i1+1, src(i1,1,1), s1, i2-i1+1,s2*s3, cudaMemcpyDevicetoHost )
    end do

    call MPI_ALLTOALLV(work1_r, decomp%x2cnts, decomp%x2disp, &
         real_type, work2_r, decomp%w2cnts, decomp%w2disp, &
         real_type, DECOMP_2D_COMM_ROW, ierror)

    istat = cudaMemcpy( dst, work2_r, d1*d2*d3, cudaMemcpyHostToDevice )
#endif

    return
  end subroutine transpose_x_to_z_real_d
#endif

  subroutine transpose_x_to_z_real(src, dst, opt_decomp)

    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

#if !defined(_MPI3)
    call MPI_Alltoallw(src,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xz, &
      dst,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xz,MPI_COMM_WORLD,ierror)
#endif

#if defined(_MPI3)
    call MPI_Neighbor_alltoallw( &
      src,decomp%xcnts_xz(decomp%xranks),decomp%xdispls_xz(decomp%xranks),decomp%xtypes_xz(decomp%xranks), &
      dst,decomp%zcnts_xz(decomp%zranks),decomp%zdispls_xz(decomp%zranks),decomp%ztypes_xz(decomp%zranks), &
      decomp%xtozNeighborComm,ierror)
#endif

    return
  end subroutine transpose_x_to_z_real!complex

