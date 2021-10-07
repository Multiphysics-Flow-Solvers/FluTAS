!
! SPDX-License-Identifier: MIT
!
module mod_sanity
  !
  use iso_c_binding , only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound     , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv    , only: chkdiv
  use mod_common_mpi, only: myid,ierr,n_z
  use mod_correc    , only: correc
  use mod_fft       , only: fftend
  use mod_fillps    , only: fillps
  use mod_initflow  , only: add_noise
  use mod_initmpi   , only: initmpi
  use mod_initsolver, only: initsolver
  use mod_param     , only: small,rho1,rho2,mu1,mu2
#if defined(_OPENACC)
  use mod_solver    , only: solver
#else
  use mod_solver_cpu, only: solver_cpu
#endif
  use mod_types
  !
  implicit none
  !
  private
  public  :: test_sanity
  !
  contains
  !
  subroutine test_sanity(ng,n,dims,ipencil,nh_d,nh_u,nh_p,halo_d,halo_u,halo_p,stop_type, &
                         cbcvel,cbcpre,bcvel,bcpre,is_outflow,is_forced,dli,dzci,dzfi)
    !
    ! performs some a priori checks of the input files before the calculation starts
    !
    implicit none
    !
    integer         , intent(in), dimension(3)                :: ng,n
    integer         , intent(in), dimension(3)                :: dims
    integer         , intent(in)                              :: ipencil
    integer         , intent(in)                              :: nh_d,nh_u,nh_p
    integer         , intent(in), dimension(2)                :: halo_d,halo_u,halo_p
    logical         , intent(in), dimension(3)                :: stop_type
    character(len=1), intent(in), dimension(0:1,3,3)          :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)            :: cbcpre
    real(rp)        , intent(in), dimension(0:1,3,3)          :: bcvel
    real(rp)        , intent(in), dimension(0:1,3)            :: bcpre
    logical         , intent(in), dimension(0:1,3)            :: is_outflow
    logical         , intent(in), dimension(3)                :: is_forced
    real(rp)        , intent(in), dimension(3)                :: dli
    real(rp)        , intent(in), dimension(1-nh_d:n(3)+nh_d) :: dzci,dzfi
    !
    logical :: passed
    !
    call chk_dims(ng,dims,ipencil,passed);         if(.not.passed) call abortit
    call chk_stop_type(stop_type,passed);          if(.not.passed) call abortit
    call chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed); if(.not.passed) call abortit
    call chk_outflow(cbcpre,is_outflow,passed);    if(.not.passed) call abortit
    call chk_forcing(cbcpre,is_forced,passed);     if(.not.passed) call abortit
    ! 
    return
  end subroutine test_sanity
  !
  subroutine chk_stop_type(stop_type,passed)
    !
    implicit none
    !
    logical, intent(in ), dimension(3) :: stop_type
    logical, intent(out)               :: passed
    !
    passed = .true.
    if(.not.any(stop_type(:))) then
      if(myid.eq.0) print*, 'ERROR: stopping criterion not chosen.'
      passed = .false.
    endif
    !
    return 
  end subroutine chk_stop_type
  !
  subroutine chk_dims(ng,dims,ipencil,passed)
    !
    implicit none
    !
    integer, intent(in ), dimension(3) :: ng
    integer, intent(in ), dimension(2) :: dims
    integer, intent(in )               :: ipencil
    logical, intent(out)               :: passed
    !
    logical :: passed_loc
    !
    passed = .true.
    passed_loc = all(mod(ng(:),2).eq.0)
    if(myid.eq.0.and.(.not.passed_loc)) &
      print*, 'ERROR: itot, jtot and ktot should be even.'
    passed = passed.and.passed_loc
    passed_loc = all(mod(ng(1:2),dims(1:2)).eq.0)
    if(myid.eq.0.and.(.not.passed_loc)) &
      print*, 'ERROR: itot and jtot should be divisable by dims(1) and dims(2), respectively.'
    passed = passed.and.passed_loc
    passed_loc = (mod(ng(2),dims(1)).eq.0).and.(mod(ng(3),dims(2)).eq.0)
    if(myid.eq.0.and.(.not.passed_loc)) &
      print*, 'ERROR: jtot should be divisable by both dims(1) and dims(2), and &
                     &ktot should be divisable by dims(2)'
    passed = passed.and.passed_loc
#if defined(_OPENACC)
    if(ipencil /= 3 .and. dims(1) /= 1) then
      if(myid == 0) print*, 'dims(1) =/ 1 not supported on a GPU run with -D_DECOMP_X or -D_DECOMP_Y.'
      if(myid == 0) print*, 'Aborting...'
      passed = .false.
      passed = passed.and.passed_loc
    endif
#endif
    !
    return
  end subroutine chk_dims
  !
  subroutine chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed)
    !
    implicit none
    !
    character(len=1), intent(in ), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in ), dimension(0:1,3  ) :: cbcpre
    real(rp)        , intent(in ), dimension(0:1,3,3) :: bcvel
    real(rp)        , intent(in ), dimension(0:1,3  ) :: bcpre
    logical         , intent(out)                     :: passed
    !
    character(len=2) :: bc01v,bc01p
    integer :: ivel,idir
    logical :: passed_loc
    !
    passed = .true.
    !
    ! check validity of pressure and velocity BCs
    !
    passed_loc = .true.
    do ivel = 1,3
      do idir=1,3
        bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
        passed_loc = passed_loc.and.( (bc01v.eq.'PP').or. &
                                      (bc01v.eq.'ND').or. &
                                      (bc01v.eq.'DN').or. &
                                      (bc01v.eq.'NN').or. &
                                      (bc01v.eq.'DD') )
      enddo
    enddo
    if(myid.eq.0.and.(.not.passed_loc)) print*, 'ERROR: velocity BCs not valid.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01p.eq.'PP').or. &
                                    (bc01p.eq.'ND').or. &
                                    (bc01p.eq.'DN').or. &
                                    (bc01p.eq.'NN').or. &
                                    (bc01p.eq.'DD') )
    enddo
    if(myid.eq.0.and.(.not.passed_loc)) print*, 'ERROR: pressure BCs not valid.' 
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      ivel = idir
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01v.eq.'PP'.and.bc01p.eq.'PP').or. &
                                    (bc01v.eq.'ND'.and.bc01p.eq.'DN').or. &
                                    (bc01v.eq.'DN'.and.bc01p.eq.'ND').or. &
                                    (bc01v.eq.'DD'.and.bc01p.eq.'NN').or. &
                                    (bc01v.eq.'NN'.and.bc01p.eq.'DD') )
    enddo
    if(myid.eq.0.and.(.not.passed_loc)) print*, 'ERROR: velocity and pressure BCs not compatible.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,2
      passed_loc = passed_loc.and.((bcpre(0,idir).eq.0.).and.(bcpre(1,idir).eq.0.))
    enddo
    if(myid.eq.0.and.(.not.passed_loc)) &
      print*, 'ERROR: pressure BCs in directions x and y must be homogeneous (value = 0.).'
    passed = passed.and.passed_loc
    !
    return 
  end subroutine chk_bc
  !
  subroutine chk_outflow(cbcpre,is_outflow,passed)
    !
    implicit none
    !
    logical         , intent(in), dimension(0:1,3) :: is_outflow
    character(len=1), intent(in), dimension(0:1,3) :: cbcpre
    logical         , intent(out) :: passed
    integer :: idir,ibound
    passed = .true.
    !
    ! 1) check for compatibility between pressure BCs and outflow BC
    !
    do idir=1,3
      do ibound = 0,1
        passed = passed.and. &
                 (cbcpre(ibound,idir).eq.'D'.and.(is_outflow(ibound,idir))) .or. &
                 (.not.is_outflow(ibound,idir))
      enddo
    enddo
    if(myid.eq.0.and.(.not.passed)) &
      print*, 'ERROR: Dirichlet pressure BC should be an outflow direction; check the BC or is_outflow in dns.in.'
    !
    return 
  end subroutine chk_outflow
  !
  subroutine chk_forcing(cbcpre,is_forced,passed)
    !
    implicit none
    !
    character(len=1), intent(in ), dimension(0:1,3) :: cbcpre
    logical         , intent(in ), dimension(3)     :: is_forced
    logical         , intent(out)                   :: passed
    !
    integer :: idir
    passed = .true.
    !
    ! 1) check for compatibility between pressure BCs and forcing BC
    !
    do idir=1,3
      if(is_forced(idir)) then
        passed = passed.and.(cbcpre(0,idir)//cbcpre(1,idir).eq.'PP')
      endif
    enddo
    if(myid.eq.0.and.(.not.passed)) &
    print*, 'ERROR: Flow cannot be forced in a non-periodic direction; check the BCs and is_forced in dns.in.'
    !
    ! 2) ensure that the thermophysical properties of the two phases are the same
    !    if bulk-velocity forcing is employed 
    !
    if(any(is_forced(:)).and.(rho1.ne.rho2)) then
      passed = .false.
    endif
    if(myid.eq.0.and.(.not.passed)) then
      print*, 'ERROR: Bulk-velocity forcing currently implemented only for matched density' 
      print*, 'As temporary alternative, use constant pressure gradient forcing (bforce(:) in dns.in)' 
    endif
    !
    return 
  end subroutine chk_forcing
  !
  subroutine chk_solvers(n,dims,dims_xyz,dli,nh_d,nh_u,nh_p,halo_u,halo_up,halo_p, &
                         dzci,dzfi,cbcvel,cbcpre,bcvel,bcpre,is_outflow,passed)
    !
    implicit none
    !
    integer         , intent(in ), dimension(3)                :: n,dims
    integer         , intent(in ), dimension(3,3)              :: dims_xyz
    real(rp)        , intent(in ), dimension(3)                :: dli
    integer         , intent(in )                              :: nh_d,nh_u,nh_p
    integer         , intent(in ), dimension(3)                :: halo_u,halo_up,halo_p
    real(rp)        , intent(in ), dimension(1-nh_d:n(3)+nh_d) :: dzci,dzfi
    character(len=1), intent(in ), dimension(0:1,3,3)          :: cbcvel
    character(len=1), intent(in ), dimension(0:1,3)            :: cbcpre
    real(rp)        , intent(in ), dimension(0:1,3,3)          :: bcvel
    real(rp)        , intent(in ), dimension(0:1,3)            :: bcpre
    logical         , intent(in ), dimension(0:1,3)            :: is_outflow
    logical         , intent(out)                              :: passed
    !
    real(rp), dimension(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) :: u,v,w
    real(rp), dimension(1-nh_p:n(1)+nh_p,1-nh_p:n(2)+nh_p,1-nh_p:n(3)+nh_p) :: p,up,vp,wp
#if defined(_OPENACC)
    attributes(managed) :: u,v,w,p,up,vp,wp,dzfi,dzci
    integer :: istat
#endif
    real(rp), dimension(1-nh_d:n(3)+nh_d) :: dzc,dzf
    type(C_PTR), dimension(2,2) :: arrplan
    real(rp), dimension(n(1),n(2)) :: lambdaxy
    real(rp) :: normfft
    real(rp), dimension(n(3)) :: a,b,c,bb
    real(rp), dimension(n(2),n(3),0:1) :: rhsbx
    real(rp), dimension(n(1),n(3),0:1) :: rhsby
    real(rp), dimension(n(1),n(2),0:1) :: rhsbz
    logical , dimension(0:1,3)         :: no_outflow
    real(rp), dimension(3) :: dl
    real(rp) :: dt,dti,alpha
    real(rp) :: divtot,divmax,resmax
    logical  :: passed_loc
#if defined(_OPENACC)
    attributes(managed) :: a,b,bb,c,lambdaxy,rhsbx,rhsby,rhsbz
#endif
    passed = .true.
    !
    ! initialize velocity below with some random noise
    !
    up(:,:,:) = 0._rp
    vp(:,:,:) = 0._rp
    wp(:,:,:) = 0._rp
    call add_noise(n(1),n(2),n(3),dims,123,0.5_rp,up(1:n(1),1:n(2),1:n(3)))
    call add_noise(n(1),n(2),n(3),dims,456,0.5_rp,vp(1:n(1),1:n(2),1:n(3)))
    call add_noise(n(1),n(2),n(3),dims,789,0.5_rp,wp(1:n(1),1:n(2),1:n(3)))
    !
    ! test pressure correction
    !
    call initsolver(n,dims,dims,dli,nh_d,dzci,dzfi,cbcpre,bcpre(:,:),(/'c','c','c'/), &
                    lambdaxy,a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
    dl  = dli**(-1)
    dzc = dzci**(-1)
    dzf = dzfi**(-1)
    dt  = acos(-1.) ! value is irrelevant
    dti = dt**(-1)
    no_outflow(:,:) = .false.
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_p,halo_up,no_outflow,dl,dzc,dzf,up,vp,wp)
#if !defined(_USE_VOF)
    call fillps(n,dli,dzfi,dti,up,vp,wp,p)
#endif
    call updt_rhs_b(n(1),n(2),n(3),(/'c','c','c'/),cbcpre,nh_p,rhsbx,rhsby,rhsbz,p(1:n(1),1:n(2),1:n(3)))
#if defined(_OPENACC)
    call solver(n_z,dims_xyz(:,3),arrplan,normfft,lambdaxy,a,b,c,cbcpre(:,:),(/'c','c','c'/),p)
#else
    call solver_cpu(n,arrplan,normfft,lambdaxy,a,b,c,cbcpre(:,3),(/'c','c','c'/),p)
#endif
#if !defined(_USE_VOF)
    call correc(n,dli,dzci,dt,p,up,vp,wp,u,v,w)
#endif
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf,u,v,w)
    call chkdiv(n(1),n(2),n(3),dli(1),dli(2),dli(3),nh_d,nh_u,dzfi,u,v,w,divtot,divmax)
    passed_loc = divmax.lt.small
    if(myid.eq.0.and.(.not.passed_loc)) &
    print*, 'ERROR: Pressure correction: Divergence is too large.'
    passed = passed.and.passed_loc
    call fftend(arrplan)
    !
    return
  end subroutine chk_solvers
  !
  subroutine abortit
    !
    implicit none
    !
    if(myid.eq.0) print*, ''
    if(myid.eq.0) print*, '*** Simulation aborted due to errors in the input file ***'
    if(myid.eq.0) print*, '    check dns.in'
    call decomp_2d_finalize
    call MPI_FINALIZE(ierr)
    call exit
    !
    return
  end subroutine abortit
  !
end module mod_sanity
