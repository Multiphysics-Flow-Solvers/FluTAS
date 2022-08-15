!
! SPDX-License-Identifier: MIT
!
module mod_initmpi
  !
  use mpi
  use decomp_2d
#if defined(_OPENACC)
  use mod_common_mpi, only: mydev
  use mod_common_mpi, only: xsl_buf, xrl_buf, xsr_buf, xrr_buf, &
                            ysr_buf, yrr_buf, ysl_buf, yrl_buf, &
                            zsr_buf, zrr_buf, zsl_buf, zrl_buf
#endif
  use mod_common_mpi, only: myid,comm_cart,ierr, &
                            left,right,front,back,top,bottom
  use mod_common_mpi, only: ijk_start,ijk_start_x,ijk_start_y,ijk_start_z, &
                            n_x,n_y,n_z,ipencil
  use mod_types
  !@cuf use cudafor
  !
  implicit none
  !
  private
#if defined(_OPENACC)
  public :: initmpi,halo_gen,alloc_buf
#else
  public :: initmpi,halo_gen
#endif
  !
  contains
  !
  subroutine initmpi(auto_tuning,ng,bc,dims_in,dims_xyz,dims,n)
    !
    implicit none
    !
    logical         , intent(in )                   :: auto_tuning
    integer         , intent(in ), dimension(3)     :: ng
    character(len=1), intent(in ), dimension(0:1,3) :: bc
    integer         , intent(in ), dimension(2)     :: dims_in
    integer         , intent(out), dimension(3,3)   :: dims_xyz
    integer         , intent(out), dimension(3)     :: dims
    integer         , intent(out), dimension(3)     :: n
    !
    logical, dimension(3) :: periods
    integer, dimension(2) :: dims_in_aux
    integer :: ntx,nty,ntz
#if defined(_OPENACC)
    integer :: dev,local_rank,local_comm,istat
#endif
    !
    ! [NOTE] I think it is much better if this is entirely controlled externally 
    !        by masking GPU devices using CUDA_VISIBLE_DEVICES   
#if defined(_OPENACC)
    ! 
    !  assign a different GPU to each MPI rank
    !  note: all the memory allocation should be dynamic, otherwise all the arrays
    !  will be allocated on device 0
    !
    call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
         MPI_INFO_NULL, local_comm, ierr)
    call MPI_Comm_rank(local_comm, mydev, ierr)
    ierr = cudaSetDevice(mydev)
    print *, " MPI rank", mydev, " using GPU", mydev
    !
#endif
    !
    periods(:) = .false.
    if(bc(0,1)//bc(1,1).eq.'PP') periods(1) = .true.
    if(bc(0,2)//bc(1,2).eq.'PP') periods(2) = .true.
    if(bc(0,3)//bc(1,3).eq.'PP') periods(3) = .true.
    !
    ! define dims_in_aux according if there is or not autotuning!
    !
    if(auto_tuning) then
      dims_in_aux(1:2) = 0
    else
      dims_in_aux(1:2) = dims_in(1:2)
    endif
    !
    ! initialize decomp_2d library (e.g., xstart, ystart, zstart)
    !
    call decomp_2d_init(ng(1),ng(2),ng(3),dims_in_aux(1),dims_in_aux(2),periods)
    myid = nrank
    ! 
    ! define dims based on the three possible decompositions
    !
    dims_xyz(1,1) = 1
    dims_xyz(2,1) = dims_in_aux(1)
    dims_xyz(3,1) = dims_in_aux(2)
    !
    dims_xyz(1,2) = dims_in_aux(1)
    dims_xyz(2,2) = 1
    dims_xyz(3,2) = dims_in_aux(2)
    !
    dims_xyz(1,3) = dims_in_aux(1)
    dims_xyz(2,3) = dims_in_aux(2)
    dims_xyz(3,3) = 1
    !
    ! define sizes of computational subdomains
    ! for the three pencil configurations
    !
    n_x(:) = ng(:)/dims_xyz(:,1)
    n_y(:) = ng(:)/dims_xyz(:,2)
    n_z(:) = ng(:)/dims_xyz(:,3)
    !
    ! define initial global index 
    ! for the three pencil configurations
    !
    ijk_start_x(:) = xstart(:) - 1
    ijk_start_y(:) = ystart(:) - 1
    ijk_start_z(:) = zstart(:) - 1
    !
    ! Three possible pencil decomposition
    !  a. _DECOMP_X: decomposition along y and z (NOT x!)
    !  b. _DECOMP_Y: decomposition along x and z (NOT y!)
    !  c. _DECOMP_Z: decomposition along x and y (NOT z!)
    !
    !  note: we typically use _DECOMP_X as it minimizes the number of transpose operations
    !
#if defined(_DECOMP_X)
    dims(:)      = dims_xyz(:,1)
    ijk_start(:) = xstart(:) - 1
    comm_cart    = DECOMP_2D_COMM_CART_X
    left         = MPI_PROC_NULL
    right        = MPI_PROC_NULL
    call MPI_CART_SHIFT(comm_cart,0,1,front,back,ierr)
    call MPI_CART_SHIFT(comm_cart,1,1,bottom,top,ierr)
    ipencil      = 1
#elif _DECOMP_Y
    dims(:)      = dims_xyz(:,2)
    ijk_start(:) = ystart(:) - 1
    comm_cart    = DECOMP_2D_COMM_CART_Y
    call MPI_CART_SHIFT(comm_cart,0,1,left,right,ierr)
    front        = MPI_PROC_NULL
    back         = MPI_PROC_NULL
    call MPI_CART_SHIFT(comm_cart,1,1,bottom,top,ierr)
    ipencil      = 2
#else /*_DECOMP_Z*/
    dims(:)      = dims_xyz(:,3)
    ijk_start(:) = zstart(:) - 1
    comm_cart    = DECOMP_2D_COMM_CART_Z
    call MPI_CART_SHIFT(comm_cart,0,1,left,right,ierr)
    call MPI_CART_SHIFT(comm_cart,1,1,front,back,ierr)
    bottom       = MPI_PROC_NULL
    top          = MPI_PROC_NULL
    ipencil      = 3
#endif
    !
    ! calculation of the local indexes
    !
    n(1) = ng(1)/dims(1)
    n(2) = ng(2)/dims(2)
    n(3) = ng(3)/dims(3)
    !
    return
  end subroutine initmpi
  !
  subroutine halo_gen(n,nh,halos)
    !
    implicit none
    !
    integer, intent(in ), dimension(3) :: n
    integer, intent(in )               :: nh
    integer, intent(out), dimension(3) :: halos
    !
    integer :: idir,ntx,nty,ntz
    !
    ntx = n(1)+2*nh
    nty = n(2)+2*nh
    ntz = n(3)+2*nh
    !
    call MPI_TYPE_VECTOR(nty*ntz,nh        ,ntx        ,MPI_REAL_RP,halos(1),ierr)
    call MPI_TYPE_VECTOR(ntz    ,nh*ntx    ,ntx*nty    ,MPI_REAL_RP,halos(2),ierr)
    call MPI_TYPE_VECTOR(1      ,nh*ntx*nty,ntx*nty*ntz,MPI_REAL_RP,halos(3),ierr)
    !
    call MPI_TYPE_COMMIT(halos(1),ierr)
    call MPI_TYPE_COMMIT(halos(2),ierr)
    call MPI_TYPE_COMMIT(halos(3),ierr)
    !
    return 
  end subroutine halo_gen
  !
#if defined(_OPENACC)
  subroutine alloc_buf(n,nh)
    !
    implicit none
    !
    integer, intent(in), dimension(3) :: n
    integer, intent(in)               :: nh
    !
#if !defined(_DECOMP_X)
    allocate( xsl_buf( 1-nh:n(2)+nh, 1-nh:n(3)+nh ) )
    allocate( xrl_buf( 1-nh:n(2)+nh, 1-nh:n(3)+nh ) )
    allocate( xsr_buf( 1-nh:n(2)+nh, 1-nh:n(3)+nh ) )
    allocate( xrr_buf( 1-nh:n(2)+nh, 1-nh:n(3)+nh ) )
#endif
    !
#if !defined(_DECOMP_Y)
    allocate( ysl_buf( 1-nh:n(1)+nh, 1-nh:n(3)+nh ) )
    allocate( yrl_buf( 1-nh:n(1)+nh, 1-nh:n(3)+nh ) )
    allocate( ysr_buf( 1-nh:n(1)+nh, 1-nh:n(3)+nh ) )
    allocate( yrr_buf( 1-nh:n(1)+nh, 1-nh:n(3)+nh ) )
#endif    
    !
#if !defined(_DECOMP_Z)
    allocate( zsl_buf( 1-nh:n(1)+nh, 1-nh:n(2)+nh ) )
    allocate( zrl_buf( 1-nh:n(1)+nh, 1-nh:n(2)+nh ) )
    allocate( zsr_buf( 1-nh:n(1)+nh, 1-nh:n(2)+nh ) )
    allocate( zrr_buf( 1-nh:n(1)+nh, 1-nh:n(2)+nh ) )
#endif
    !
    ! [TODO]: wrap those in a high-level module
    ! [TODO]: to review how effective these API are
    !
#if defined(_OPENACC) && !defined(_GPU_MPI)
    !
#if !defined(_DECOMP_X)
    istat = cudaMemAdvise(        xsl_buf, size(xsl_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
    istat = cudaMemAdvise(        xsl_buf, size(xsl_buf), cudaMemAdviseSetAccessedBy, mydev )
    istat = cudaMemPrefetchAsync( xsl_buf, size(xsl_buf), cudaCpuDeviceId, 0 )
    istat = cudaMemAdvise(        xrl_buf, size(xrl_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
    istat = cudaMemAdvise(        xrl_buf, size(xrl_buf), cudaMemAdviseSetAccessedBy, mydev )
    istat = cudaMemPrefetchAsync( xrl_buf, size(xrl_buf), cudaCpuDeviceId, 0 )
    istat = cudaMemAdvise(        xsr_buf, size(xsr_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
    istat = cudaMemAdvise(        xsr_buf, size(xsr_buf), cudaMemAdviseSetAccessedBy, mydev )
    istat = cudaMemPrefetchAsync( xsr_buf, size(xsr_buf), cudaCpuDeviceId, 0 )
    istat = cudaMemAdvise(        xrr_buf, size(xrr_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
    istat = cudaMemAdvise(        xrr_buf, size(xrr_buf), cudaMemAdviseSetAccessedBy, mydev )
    istat = cudaMemPrefetchAsync( xrr_buf, size(xrr_buf), cudaCpuDeviceId, 0 )
#endif
    !
#if !defined(_DECOMP_Y)
    istat = cudaMemAdvise(        ysl_buf, size(ysl_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
    istat = cudaMemAdvise(        ysl_buf, size(ysl_buf), cudaMemAdviseSetAccessedBy, mydev )
    istat = cudaMemPrefetchAsync( ysl_buf, size(ysl_buf), cudaCpuDeviceId, 0 )
    istat = cudaMemAdvise(        yrl_buf, size(yrl_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
    istat = cudaMemAdvise(        yrl_buf, size(yrl_buf), cudaMemAdviseSetAccessedBy, mydev )
    istat = cudaMemPrefetchAsync( yrl_buf, size(yrl_buf), cudaCpuDeviceId, 0 )
    istat = cudaMemAdvise(        ysr_buf, size(ysr_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
    istat = cudaMemAdvise(        ysr_buf, size(ysr_buf), cudaMemAdviseSetAccessedBy, mydev )
    istat = cudaMemPrefetchAsync( ysr_buf, size(ysr_buf), cudaCpuDeviceId, 0 )
    istat = cudaMemAdvise(        yrr_buf, size(yrr_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
    istat = cudaMemAdvise(        yrr_buf, size(yrr_buf), cudaMemAdviseSetAccessedBy, mydev )
    istat = cudaMemPrefetchAsync( yrr_buf, size(yrr_buf), cudaCpuDeviceId, 0 )
#endif
    !
#if !defined(_DECOMP_Z)
    istat = cudaMemAdvise(        zsl_buf, size(zsl_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
    istat = cudaMemAdvise(        zsl_buf, size(zsl_buf), cudaMemAdviseSetAccessedBy, mydev )
    istat = cudaMemPrefetchAsync( zsl_buf, size(zsl_buf), cudaCpuDeviceId, 0 )
    istat = cudaMemAdvise(        zrl_buf, size(zrl_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
    istat = cudaMemAdvise(        zrl_buf, size(zrl_buf), cudaMemAdviseSetAccessedBy, mydev )
    istat = cudaMemPrefetchAsync( zrl_buf, size(zrl_buf), cudaCpuDeviceId, 0 )
    istat = cudaMemAdvise(        zsr_buf, size(zsr_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
    istat = cudaMemAdvise(        zsr_buf, size(zsr_buf), cudaMemAdviseSetAccessedBy, mydev )
    istat = cudaMemPrefetchAsync( zsr_buf, size(zsr_buf), cudaCpuDeviceId, 0 )
    istat = cudaMemAdvise(        zrr_buf, size(zrr_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
    istat = cudaMemAdvise(        zrr_buf, size(zrr_buf), cudaMemAdviseSetAccessedBy, mydev )
    istat = cudaMemPrefetchAsync( zrr_buf, size(zrr_buf), cudaCpuDeviceId, 0 )
#endif
    !
#endif
    !
    return
  end subroutine alloc_buf
#endif
  !
end module mod_initmpi
