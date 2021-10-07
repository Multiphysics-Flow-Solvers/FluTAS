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

! This file contains the routines that transpose data from Z to X pencil

#if defined(_OPENACC)
  subroutine transpose_z_to_xc(src, dst, pxc, opt_decomp)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN) :: src
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
    integer ::           kstart,kend,decomp_z2distm
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
       pos = decomp%w2disp(m) + 1
       istat = cudaMemcpyAsync( work1_r_d(pos), src(1,1,decomp%z2idx(m)), decomp%w2cnts(m), a2a_d2h )
       istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
    end do

    ! self
    m = row_rank
#if !defined(_EPHC)
    istat = cudaMemcpy2DAsync( dst(decomp%x2idx(m),1,1), d1, src(1,1,decomp%z2idx(m)), s1, decomp%x2dist(m), d2*d3, stream=a2a_comp )
#else
       istart = decomp%x2idx(m)
       iend   = istart + decomp%x2dist(m) - 1
       decomp_x2distm = decomp%x2dist(m)
       kstart = decomp%z2idx(m)
       kend   = kstart + decomp%z2dist(m) - 1
       !$cuf kernel do(3) <<<*,*,stream=a2a_comp>>>
       do k=kstart,kend
       do j=1,d2
       do i=istart, iend
         pxc(i,j,(k-kstart+1)) = src((i-istart+1),j,k) 
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
      call MPI_IRECV( work2_r_d(decomp%x2disp(m)+1), decomp%x2cnts(m), real_type, m, 0, DECOMP_2D_COMM_ROW, a2a_requests(iter), ierror)
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
       call MPI_SEND( work1_r_d(decomp%w2disp(m)+1), decomp%w2cnts(m), real_type, m, 0, DECOMP_2D_COMM_ROW, ierror)
       call MPI_WAIT(a2a_requests(iter), MPI_STATUS_IGNORE, ierror)
       m = sorc
       pos = decomp%x2disp(m) + 1
#if !defined(_EPHC)
       istat = cudaMemcpy2DAsync( dst(decomp%x2idx(m),1,1), d1, work2_r_d(pos), decomp%x2dist(m), decomp%x2dist(m), d2*d3, stream=a2a_comp )
#else
       istart = decomp%x2idx(m)
       iend   = istart + decomp%x2dist(m) - 1
       decomp_x2distm = decomp%x2dist(m)
       !$cuf kernel do(3) <<<*,*,stream=a2a_comp>>>
       do k=1,d3
       do j=1,d2
       do i=istart, iend
         idx = (i-istart+1) + decomp_x2distm*(j-1) + decomp_x2distm*d2*(k-1)
         pxc(i,j,k) = work2_r_d(pos + idx - 1) 
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
      pxc(i,j,k) = dst(i,j,k)
    end do
    end do
    end do
#endif

    istat = cudaEventRecord( a2a_event(0), 0 )
    istat = cudaEventSynchronize( a2a_event(0) )
#else
    istat = cudaMemcpy( work1_r, src, s1*s2*s3, cudaMemcpyDeviceToHost )

    call MPI_ALLTOALLV(work1_r, decomp%w2cnts, decomp%w2disp, &
         real_type, work2_r, decomp%x2cnts, decomp%x2disp, &
         real_type, DECOMP_2D_COMM_ROW, ierror)

    ! rearrange receive buffer
    do m=0,dims(2)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%x2dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%x2dist(m)-1
       end if
       pos = decomp%x2disp(m) + 1
       istat = cudaMemcpy2D( dst(i1,1,1), d1, work2_r(pos), i2-i1+1, i2-i1+1, d2*d3, cudaMemcpyHostToDevice )
    end do
#endif


    return
  end subroutine transpose_z_to_xc

  subroutine transpose_z_to_x_real_d(src, dst, opt_decomp)

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
       pos = decomp%w2disp(m) + 1
       istat = cudaMemcpyAsync( work1_r_d(pos), src(1,1,decomp%z2idx(m)), decomp%w2cnts(m), a2a_d2h )
       istat = cudaEventRecord( a2a_event(iter), a2a_d2h )
    end do

    ! self
    m = row_rank
    istat = cudaMemcpy2DAsync( dst(decomp%x2idx(m),1,1), d1, src(1,1,decomp%z2idx(m)), s1, decomp%x2dist(m), d2*d3, stream=a2a_comp )

    do iter=1,dims(2)-1
      if( pow2 ) then
        sorc = IEOR(row_rank,iter)
      else
        sorc = mod(row_rank - iter + dims(2), dims(2))
      endif
      m = sorc
      call MPI_IRECV( work2_r_d(decomp%x2disp(m)+1), decomp%x2cnts(m), real_type, m, 0, DECOMP_2D_COMM_ROW, a2a_requests(iter), ierror)
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
       call MPI_SEND( work1_r_d(decomp%w2disp(m)+1), decomp%w2cnts(m), real_type, m, 0, DECOMP_2D_COMM_ROW, ierror)
       call MPI_WAIT(a2a_requests(iter), MPI_STATUS_IGNORE, ierror)
       m = sorc
       pos = decomp%x2disp(m) + 1
       istat = cudaMemcpy2DAsync( dst(decomp%x2idx(m),1,1), d1, work2_r_d(pos), decomp%x2dist(m), decomp%x2dist(m), d2*d3, stream=a2a_comp )
    end do
    istat = cudaEventRecord( a2a_event(0), 0 )
    istat = cudaEventSynchronize( a2a_event(0) )
#else
    istat = cudaMemcpy( work1_r, src, s1*s2*s3, cudaMemcpyDeviceToHost )

    call MPI_ALLTOALLV(work1_r, decomp%w2cnts, decomp%w2disp, &
         real_type, work2_r, decomp%x2cnts, decomp%x2disp, &
         real_type, DECOMP_2D_COMM_ROW, ierror)

    ! rearrange receive buffer
    do m=0,dims(2)-1
       if (m==0) then
          i1 = 1
          i2 = decomp%x2dist(0)
       else
          i1 = i2+1
          i2 = i1+decomp%x2dist(m)-1
       end if
       pos = decomp%x2disp(m) + 1
       istat = cudaMemcpy2D( dst(i1,1,1), d1, work2_r(pos), i2-i1+1, i2-i1+1, d2*d3, cudaMemcpyHostToDevice )
    end do
#endif


    return
  end subroutine transpose_z_to_x_real_d
#endif

  subroutine transpose_z_to_x_real(src, dst, opt_decomp)

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
    call MPI_Alltoallw(src,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xz, &
      dst,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xz,MPI_COMM_WORLD,ierror)
#endif

#if defined(_MPI3)
    call MPI_Neighbor_alltoallw( &
      src,decomp%zcnts_xz(decomp%zranks),decomp%zdispls_xz(decomp%zranks),decomp%ztypes_xz(decomp%zranks), &
      dst,decomp%xcnts_xz(decomp%xranks),decomp%xdispls_xz(decomp%xranks),decomp%xtypes_xz(decomp%xranks), &
      decomp%ztoxNeighborComm,ierror)
#endif

    return
  end subroutine transpose_z_to_x_real!complex

