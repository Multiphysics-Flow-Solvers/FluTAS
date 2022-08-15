!
! SPDX-License-Identifier: MIT
!
module mod_param
  !
  ! for easy to read we recommend to:
  !    1. First, declare the variables to be determined from '*.in' files;
  !    2. Next, declare the auxiliaries variable not determined from '*.in' files;
  !    3. Finally, declare the parameter variables.
  !
  use mod_types
  !
  implicit none
  !
  public
  !
  ! 1. variables to be determined from the input files '*.in'
  !
  !  --> from dns.in (always defined in this app)
  !
  integer                              :: itot,jtot,ktot
  real(rp)                             :: lx,ly,lz
  real(rp)                             :: gr
  real(rp)                             :: cfl,dt_input
  logical                              :: constant_dt
  character(len=3)                     :: time_scheme,space_scheme_mom
  real(rp)                             :: rho_sp,mu_sp
  character(len=100)                   :: inivel
  logical                              :: is_wallturb
  character(len=3)                     :: wallturb_type
  character(len=3)                     :: bulk_ftype
  integer                              :: nstep
  real(rp)                             :: time_max,tw_max
  logical, dimension(3)                :: stop_type
  logical                              :: restart
  integer                              :: icheck,iout0d,iout1d,iout2d,iout3d,isave
  character(len=1), dimension(0:1,3,3) :: cbcvel
  real(rp)        , dimension(0:1,3,3) :: bcvel
  character(len=1), dimension(0:1,3)   :: cbcpre
  real(rp)        , dimension(0:1,3)   :: bcpre
  logical         , dimension(3)       :: is_forced
  real(rp)                             :: gacc_x,gacc_y,gacc_z
  real(rp)                             :: bvel_x,bvel_y,bvel_z
  real(rp)                             :: dpdl_x,dpdl_y,dpdl_z
  logical, dimension(0:1,3)            :: is_outflow
  integer, dimension(2)                :: dims_in  ! input (user-choice for a 2D decomp.)
  integer                              :: nthreadsmax
  !
  !  --> from vof.in (always defined in this app)
  !
  real(rp) :: rho1,rho2,mu1,mu2
  character(len=3) :: inivof
  real(rp),  allocatable, dimension(:) :: xc, yc, zc, r
  character(len=1), dimension(0:1,3) :: cbcvof
  real(rp)        , dimension(0:1,3) ::  bcvof
  real(rp) :: sigma
  logical  :: late_init
  integer  :: nbub=0, i_b,i_late_init
  !
#if defined(_TURB_FORCING)
  !
  !  --> from forcing.in (if defined)
  !
  character(len=3) :: turb_type
  real(rp)         :: u0_t, f0_t, k0_t
  real(rp)         :: abc_x,abc_y,abc_z
  logical          :: add_noise_abc
  !
#endif
  !
#if defined(_USE_VOF) && (_DO_POSTPROC) 
  !
  !  --> from post_vof.in (if defined)
  !
  logical :: do_tagging
  integer :: iout0d_ta
  !
#endif
  !
  ! 2. auxiliaries variables not determined from the input files
  !
  logical :: exists
  integer , dimension(3) :: ng
  integer , dimension(3) :: n
  real(rp), dimension(3) :: l
  real(rp), dimension(3) :: dl
  real(rp), dimension(3) :: dli
  real(rp) :: dxi,dyi,dzi,dx,dy,dz
  real(rp) :: rho0
  integer  :: n_stage
  real(rp) :: cfl_c,cfl_d
  !
  ! 3. parameters, e.g., pi, RK coefficients, small, universal constants etc.
  !
  real(rp), parameter :: pi = acos(-1._rp)
  real(rp), parameter :: small = epsilon(pi)*10**(precision(pi)/2._rp)
  logical , parameter, dimension(2,3) :: no_outflow = & 
      reshape((/.false.,.false.,   & ! no outflow in x lower,upper bound
                .false.,.false.,   & ! no outflow in y lower,upper bound
                .false.,.false./), & ! no outflow in z lower,upper bound
                shape(no_outflow))
  real(rp), parameter, dimension(2,3) :: rkcoeff = reshape((/ 32._rp/60._rp,  0._rp        , &
                                                              25._rp/60._rp, -17._rp/60._rp, &
                                                              45._rp/60._rp, -25._rp/60._rp/), shape(rkcoeff))
  real(rp), parameter, dimension(3)   :: rkcoeff12 = rkcoeff(1,:)+rkcoeff(2,:)
  !
  contains
  ! 
  subroutine read_input(myid)
    !
    use mpi
    !
    implicit none
    !
    integer, intent(in) :: myid
    integer :: iunit,ierr
    !
    ! load input files from dns.in
    !
    open(newunit=iunit,file='dns.in',status='old',action='read',iostat=ierr)
    if( ierr.eq.0 ) then
      read(iunit,*) itot,jtot,ktot
      read(iunit,*) lx,ly,lz
      read(iunit,*) gr
      read(iunit,*) cfl,dt_input
      read(iunit,*) constant_dt
      read(iunit,*) time_scheme, space_scheme_mom
      read(iunit,*) rho_sp, mu_sp
      read(iunit,*) inivel
      read(iunit,*) is_wallturb,wallturb_type
      read(iunit,*) bulk_ftype
      read(iunit,*) nstep,time_max,tw_max
      read(iunit,*) stop_type(1),stop_type(2),stop_type(3)
      read(iunit,*) restart
      read(iunit,*) icheck,iout0d,iout1d,iout2d,iout3d,isave
      read(iunit,*) cbcvel(0,1,1),cbcvel(1,1,1),cbcvel(0,2,1),cbcvel(1,2,1),cbcvel(0,3,1),cbcvel(1,3,1)
      read(iunit,*) cbcvel(0,1,2),cbcvel(1,1,2),cbcvel(0,2,2),cbcvel(1,2,2),cbcvel(0,3,2),cbcvel(1,3,2)
      read(iunit,*) cbcvel(0,1,3),cbcvel(1,1,3),cbcvel(0,2,3),cbcvel(1,2,3),cbcvel(0,3,3),cbcvel(1,3,3)
      read(iunit,*) cbcpre(0,1  ),cbcpre(1,1  ),cbcpre(0,2  ),cbcpre(1,2  ),cbcpre(0,3  ),cbcpre(1,3  )
      read(iunit,*)  bcvel(0,1,1), bcvel(1,1,1), bcvel(0,2,1), bcvel(1,2,1), bcvel(0,3,1), bcvel(1,3,1)
      read(iunit,*)  bcvel(0,1,2), bcvel(1,1,2), bcvel(0,2,2), bcvel(1,2,2), bcvel(0,3,2), bcvel(1,3,2)
      read(iunit,*)  bcvel(0,1,3), bcvel(1,1,3), bcvel(0,2,3), bcvel(1,2,3), bcvel(0,3,3), bcvel(1,3,3)
      read(iunit,*)  bcpre(0,1  ), bcpre(1,1  ), bcpre(0,2  ), bcpre(1,2  ), bcpre(0,3  ), bcpre(1,3  )
      read(iunit,*)  is_forced(1),is_forced(2),is_forced(3)
      read(iunit,*)  gacc_x,gacc_y,gacc_z
      read(iunit,*)  bvel_x,bvel_y,bvel_z
      read(iunit,*)  dpdl_x,dpdl_y,dpdl_z
      read(iunit,*)  is_outflow(0,1),is_outflow(1,1),is_outflow(0,2),is_outflow(1,2),is_outflow(0,3),is_outflow(1,3)
      read(iunit,*) dims_in(1),dims_in(2)
      read(iunit,*) nthreadsmax
    else
      if(myid.eq.0) print*, 'Error reading the dns.in input file' 
      if(myid.eq.0) print*, 'Input file missing or incomplete' 
      if(myid.eq.0) print*, 'Aborting...'
      call MPI_FINALIZE(ierr)
      call exit
    endif
    close(iunit)
    !  
    ! load two-phase input files
    !
    open(newunit=iunit,file='vof.in',status='old',action='read',iostat=ierr)
    if( ierr.eq.0 ) then
      read(iunit,*) rho1, rho2, mu1, mu2
      read(iunit,*) inivof
      read(iunit,*) nbub
      allocate(xc(nbub), yc(nbub), zc(nbub), r(nbub))
      if(nbub.eq.1) then
        read(iunit,*) xc, yc, zc, r
      else
        read(iunit,*) xc(1), yc(1), zc(1), r(1)
      endif
      read(iunit,*) cbcvof(0,1 ),cbcvof(1,1 ),cbcvof(0,2 ),cbcvof(1,2 ),cbcvof(0,3 ),cbcvof(1,3 )
      read(iunit,*) bcvof(0,1  ),bcvof(1,1  ),bcvof(0,2  ),bcvof(1,2  ),bcvof(0,3  ),bcvof(1,3  )
      read(iunit,*) sigma
      read(iunit,*) late_init, i_late_init
    else
      if(myid.eq.0) print*, 'Error reading the vof.in input file' 
      if(myid.eq.0) print*, 'Input file missing or incomplete' 
      if(myid.eq.0) print*, 'Aborting...'
      call MPI_FINALIZE(ierr)
      call exit
    endif
    rho0 = min(rho1,rho2)
    close(iunit)
    !
    if(nbub.gt.1) then
      open(newunit=iunit,file='bub.in',status='old',action='read',iostat=ierr)
      if( ierr.eq.0 ) then    
        do i_b=1,nbub
          read(iunit,*) xc(i_b), yc(i_b), zc(i_b), r(i_b)
        enddo
      else
        if(myid.eq.0) print*, 'Error reading the bub.in input file' 
        if(myid.eq.0) print*, 'Input file missing or incomplete' 
        if(myid.eq.0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        call exit
      endif
      close(iunit)
    endif
    !
    ! load input files from single-phase postprocessing (if defined)
    !
#if defined(_USE_VOF) && defined(_DO_POSTPROC)
    inquire(file='post_vof.in',exist=exists)
    if(exists) then
      open(newunit=iunit,file='post_vof.in',status='old',action='read',iostat=ierr)
      if( ierr.eq.0 ) then
        read(iunit,*) do_tagging,iout0d_ta
      else
        if(myid.eq.0) print*, 'Error reading the post_vof.in input file' 
        if(myid.eq.0) print*, 'Input file missing or incomplete' 
        if(myid.eq.0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        call exit
      endif
    else
      do_tagging = .false.
      iout0d_ta  = 1
    endif
    close(iunit)
#endif
    !  
    ! load turbulent forcing input files (if defined)
    !
#if defined(_TURB_FORCING)
    open(newunit=iunit,file='forcing.in',status='old',action='read',iostat=ierr)
    if( ierr.eq.0 ) then
      read(iunit,*) turb_type
      read(iunit,*) u0_t   
      read(iunit,*) f0_t   
      read(iunit,*) k0_t  
      read(iunit,*) abc_x, abc_y, abc_z
      read(iunit,*) add_noise_abc
    else
      if(myid.eq.0) print*, 'Error reading the forcing.in input file' 
      if(myid.eq.0) print*, 'Input file missing or incomplete' 
      if(myid.eq.0) print*, 'Aborting...'
      call MPI_FINALIZE(ierr)
      call exit
    endif
    close(iunit)
#endif
    !  
    ! compute the cartesian spacings
    !
    dx  = lx/(1._rp*itot)
    dy  = ly/(1._rp*jtot)
    dz  = lz/(1._rp*ktot)
    dxi = dx**(-1)
    dyi = dy**(-1)
    dzi = dz**(-1)
    !
    ng  = (/itot,jtot,ktot/)
    l   = (/lx  ,ly  ,lz  /)
    dl  = (/dx  ,dy  ,dz  /)
    dli = (/dxi ,dyi ,dzi /)
    !
    ! for consistency, we replace rho_sp and mu_sp with rho2 and mu2.
    ! Anyway, rho_sp and mu_sp are not employed if (_USE_VOF) flag is active
    !
    rho_sp = rho2 
    mu_sp  = mu2
    !
    ! compute the cfl and other stuff related to
    ! the time discretization
    !
    if(    time_scheme.eq.'ab2') then
      n_stage = 1
      cfl_c   = 1._rp
      cfl_d   = 1._rp/6._rp
    elseif(time_scheme.eq.'rk3') then
      n_stage = 3
      cfl_c   = sqrt(3._rp)
      cfl_d   = 1.65_rp/12._rp
    endif
    !
    return
  end subroutine read_input
  !
end module mod_param
