!
! SPDX-License-Identifier: MIT
!
module mod_param
  !
  use mod_types
  !
  implicit none
  !
  public
  !
  ! parameters
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
  ! variables to be determined from the input files '*.in'
  !
  integer                              :: itot,jtot,ktot
  real(rp)                             :: lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,gr
  real(rp)                             :: cfl,dt_input
  logical                              :: constant_dt
  real(rp)                             :: visc
  character(len=100)                   :: inivel
  integer                              :: is_wallturb
  integer                              :: nstep
  real(rp)                             :: time_max,tw_max
  logical, dimension(3)                :: stop_type
  logical                              :: restart
  integer                              :: icheck,iout0d,iout1d,iout2d,iout3d,isave
  integer, dimension(2)                :: dims_in  ! input (user-choice for a 2D decomp.)
  integer                              :: nthreadsmax
  character(len=1), dimension(0:1,3,3) :: cbcvel
  real(rp)        , dimension(0:1,3,3) :: bcvel
  character(len=1), dimension(0:1,3)   :: cbcpre
  real(rp)        , dimension(0:1,3)   :: bcpre
  !
  real(rp), dimension(3) :: bforce
  real(rp) :: bulk_velx,bulk_vely,bulk_velz
  logical , dimension(3) :: is_forced
  logical , dimension(0:1,3) :: is_outflow
  !
  integer , dimension(3) :: ng
  integer , dimension(3) :: n
  real(rp), dimension(3) :: l
  real(rp), dimension(3) :: dl
  real(rp), dimension(3) :: dli
  ! 
#if defined(_DO_POSTPROC)
  logical  :: do_avg,do_favre,do_wall
  real(rp) :: deltaT
  integer  :: avg_dir,time_deltai,wall_deltai
#endif
  ! 
#if defined(_USE_VOF)
  real(rp),  dimension(3) :: gacc 
  character(len=8) :: inivof
  character(len=1), dimension(0:1,3) ::  cbcvof
  real(rp),  allocatable, dimension(:) :: xc, yc, zc, r
  real(rp) :: rho1,rho2,mu1,mu2,sigma,rho0
  real(rp), dimension(0:1,3) ::   bcvof
  real(rp),  dimension(3) :: dpdl
  logical :: late_init
  integer :: nbub=0, i_b,i_late_init
#endif
  !
#if defined(_TURB_FORCING)
  character(len=3)       :: turbType
  real(rp)               :: u0_t, f0_t, k0_t
  real(rp), dimension(3) :: abc
  logical :: add_noise_abc
#endif
#if defined(_HEAT_TRANSFER)
  character(len=3)   :: initmp
  real(rp)           :: cp1,cv1,kappa1
  real(rp)           :: cp2,cv2,kappa2
  real(rp)           :: tl0,tg0,gam_g
#if defined(_BOUSSINESQ)
  real(rp)           :: tmp0,beta_th
#endif
  character(len=1), dimension(0:1,3) :: cbctmp
  real(rp)        , dimension(0:1,3) ::  bctmp
#endif
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
    ! load single-phase input files
    !
    open(newunit=iunit,file='dns.in',status='old',action='read',iostat=ierr)
    if( ierr.eq.0 ) then
      read(iunit,*) itot,jtot,ktot
      read(iunit,*) lx,ly,lz
      read(iunit,*) gr
      read(iunit,*) cfl,dt_input
      read(iunit,*) constant_dt
      read(iunit,*) visc
      read(iunit,*) inivel
      read(iunit,*) is_wallturb
      read(iunit,*) nstep, time_max,tw_max
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
      read(iunit,*)  bforce(1),bforce(2),bforce(3)
      read(iunit,*)  is_forced(1),is_forced(2),is_forced(3)
      read(iunit,*)  bulk_velx,bulk_vely,bulk_velz
      read(iunit,*)  is_outflow(0,1),is_outflow(1,1),is_outflow(0,2),is_outflow(1,2),is_outflow(0,3),is_outflow(1,3)
      read(iunit,*) dims_in(1),dims_in(2)
      read(iunit,*) nthreadsmax
    else
      if(myid.eq.0) print*, 'Error reading the dns.in input file' 
      if(myid.eq.0) print*, 'Aborting...'
      call MPI_FINALIZE(ierr)
      call exit
    endif
    close(iunit)
    !  
    ! load two-phase input files
    !
#if defined(_USE_VOF)
    open(newunit=iunit,file='vof.in',status='old',action='read',iostat=ierr)
    if( ierr.eq.0 ) then
      read(iunit,*) rho1, rho2, mu1,mu2
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
      read(iunit,*) gacc(1), gacc(2), gacc(3)
      read(iunit,*) sigma
      read(iunit,*) dpdl(1), dpdl(2), dpdl(3)
      read(iunit,*) late_init, i_late_init
    else
      if(myid.eq.0) print*, 'Error reading the vof.in input file' 
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
        if(myid.eq.0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        call exit
      endif
      close(iunit)
    endif
#endif
    !  
    ! load post-processing input files
    !
#if defined(_DO_POSTPROC)
    open(newunit=iunit,file='post.in',status='old',action='read',iostat=ierr)
    if( ierr.eq.0 ) then
      read(iunit,*) do_avg
      read(iunit,*) avg_dir 
      read(iunit,*) time_deltai
      read(iunit,*) do_favre
      read(iunit,*) do_wall
      read(iunit,*) wall_deltai
    else
      if(myid.eq.0) print*, 'Error reading the post.in input file' 
      if(myid.eq.0) print*, 'Aborting...'
      call MPI_FINALIZE(ierr)
      call exit
    endif
    close(iunit)
#endif
    !  
    ! load turbulent forcing input files
    !
#if defined(_TURB_FORCING)
    open(newunit=iunit,file='forcing.in',status='old',action='read',iostat=ierr)
    if( ierr.eq.0 ) then
      read(iunit,*) turbType
      read(iunit,*) u0_t   
      read(iunit,*) f0_t   
      read(iunit,*) k0_t  
      read(iunit,*) abc(1), abc(2), abc(3)
      read(iunit,*) add_noise_abc
    else
      if(myid.eq.0) print*, 'Error reading the forcing.in input file' 
      if(myid.eq.0) print*, 'Aborting...'
      call MPI_FINALIZE(ierr)
      call exit
    endif
    close(iunit)
#endif
    !  
    ! load heat transfer input files
    !
#if defined(_HEAT_TRANSFER)
    open(newunit=iunit,file='heat_transfer.in',status='old',action='read',iostat=ierr)
    if( ierr.eq.0 ) then
      read(iunit,*) initmp
      read(iunit,*) tl0,tg0
      read(iunit,*) cp1,cp2
      read(iunit,*) cv1,cv2
      read(iunit,*) kappa1,kappa2
      read(iunit,*) cbctmp(0,1 ),cbctmp(1,1 ),cbctmp(0,2 ),cbctmp(1,2 ),cbctmp(0,3 ),cbctmp(1,3 )
      read(iunit,*) bctmp(0,1  ),bctmp(1,1  ),bctmp(0,2  ),bctmp(1,2  ),bctmp(0,3  ),bctmp(1,3  )
#if defined(_BOUSSINESQ)
      read(iunit,*) tmp0,beta_th
#endif
      gam_g = cp1/cv1
#if defined(_DO_POSTPROC)
      if(avg_dir.eq.1) then
        deltaT=abs(bctmp(0,1)-bctmp(1,1))
      elseif(avg_dir.eq.2) then
        deltaT=abs(bctmp(0,2)-bctmp(1,2))
      else
        deltaT=abs(bctmp(0,3)-bctmp(1,3))
      endif
#endif
      else
        if(myid.eq.0) print*, 'Error reading the heat_transfer.in input file' 
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
    ! if _USE_VOF is active, we replace visc with the minimum viscosity between the two phases
    ! note: in this case, please fill dns.in and vof.in consistently
    !
#if defined(_USE_VOF)
    visc = min(mu1/rho1,rho2/mu2)
#endif
    !
    return
  end subroutine read_input
  !
end module mod_param
