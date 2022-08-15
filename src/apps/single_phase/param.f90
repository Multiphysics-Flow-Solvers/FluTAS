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
#if defined(_HEAT_TRANSFER)
  !
  !  --> from heat_transfer_sp.in (if defined)
  !
  character(len=3) :: initmp
  real(rp) :: tmp0
  real(rp) :: cp_sp,cv_sp,kappa_sp
#if defined(_BOUSSINESQ)
  real(rp) :: beta_sp_th
#endif
  character(len=1), dimension(0:1,3) :: cbctmp
  real(rp)        , dimension(0:1,3) ::  bctmp
  !
#endif
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
#if defined(_DO_POSTPROC)
  !
  !  --> from post_sp.in (if defined)
  !
  real(rp) :: deltaT
  logical  :: do_avg,do_favre,do_wall
  integer  :: avg_dir,time_deltai,wall_deltai
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
  character(len=3) :: time_scheme,space_scheme_mom
  integer  :: n_stage
  real(rp) :: cfl_c,cfl_d
#if defined(_HEAT_TRANSFER)
  real(rp) :: gam_g
#endif
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
    ! load input files from heat_transfer_sp.in (if defined)
    !
#if defined(_HEAT_TRANSFER)
    open(newunit=iunit,file='heat_transfer_sp.in',status='old',action='read',iostat=ierr)
    if( ierr.eq.0 ) then
      read(iunit,*) initmp
      read(iunit,*) tmp0
      read(iunit,*) cp_sp
      read(iunit,*) cv_sp
      read(iunit,*) kappa_sp
      read(iunit,*) cbctmp(0,1 ),cbctmp(1,1 ),cbctmp(0,2 ),cbctmp(1,2 ),cbctmp(0,3 ),cbctmp(1,3 )
      read(iunit,*) bctmp(0,1  ),bctmp(1,1  ),bctmp(0,2  ),bctmp(1,2  ),bctmp(0,3  ),bctmp(1,3  )
#if defined(_BOUSSINESQ)
      read(iunit,*) beta_sp_th
#endif
      else
        if(myid.eq.0) print*, 'Error reading the heat_transfer_sp.in input file' 
        if(myid.eq.0) print*, 'Input file missing or incomplete' 
        if(myid.eq.0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        call exit
      endif
    close(iunit)
    gam_g = cp_sp/cv_sp
#endif
    !
    ! load input files from single-phase postprocessing (if defined)
    ! 
#if defined(_DO_POSTPROC)
    open(newunit=iunit,file='post_sp.in',status='old',action='read',iostat=ierr)
    if( ierr.eq.0 ) then
      read(iunit,*) do_avg
      read(iunit,*) avg_dir 
      read(iunit,*) time_deltai
      read(iunit,*) do_favre
      read(iunit,*) do_wall
      read(iunit,*) wall_deltai
    else
      if(myid.eq.0) print*, 'Error reading the post.in input file' 
      if(myid.eq.0) print*, 'Input file missing or incomplete' 
      if(myid.eq.0) print*, 'Aborting...'
      call MPI_FINALIZE(ierr)
      call exit
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
#if defined(_DO_POSTPROC)
    if(    avg_dir.eq.1) then
      deltaT = abs(bctmp(0,1)-bctmp(1,1))
    elseif(avg_dir.eq.2) then
      deltaT = abs(bctmp(0,2)-bctmp(1,2))
    elseif(avg_dir.eq.3) then
      deltaT = abs(bctmp(0,3)-bctmp(1,3))
    else
      if(myid.eq.0) print*, 'Invalid averaging direction. Please provide one in post_sp.in' 
      call MPI_FINALIZE(ierr)
      call exit
    endif
#endif
    !
    return
  end subroutine read_input
  !
end module mod_param
