!
! SPDX-License-Identifier: MIT
!
!---------------------------------------------------------------------------------
!
! FluTAS -- Fluid Transport Accelerated Solver                            
!                                                                                 
!  a. Menu Title: FluTAS_VoF_inc;                                                 
!  b. Feautures of FluTAS_VoF_inc:                                                
!      --> two-fluid incompressible and adiabatic solver with VoF;
!      --> allow for density and viscosity mismatch;                              
!      --> momentum equation advanced with Adams-Bashforth (explicit diffusion);  
!      --> pressure equation solved with a FFT-based direct solver.               
! 
!---------------------------------------------------------------------------------
!
program flutas
  !
  ! module declaration 
  !  note: --> import what you really neeed 
  !
  use iso_c_binding , only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound     , only: boundp,updt_rhs_b,bounduvw
  use mod_chkdiv    , only: chkdiv
  use mod_chkdt     , only: chkdt
  use mod_common_mpi, only: myid,ierr,comm_cart,n_z,ijk_start,ipencil
  use mod_correc    , only: correc
  use mod_debug     , only: chkmean
  use mod_fft       , only: fftini,fftend
  use mod_fillps    , only: fillps
  use mod_initflow  , only: initflow
  use mod_initgrid  , only: initgrid
  use mod_initmpi   , only: initmpi,halo_gen
#if defined(_OPENACC)
  use mod_initmpi   , only: alloc_buf
#endif
  use mod_initsolver, only: initsolver
  use mod_load      , only: load, load_scalar
  use mod_rk        , only: rk
  use mod_output    , only: out0d,out1d,out1d_2,out2d,out3d,write_visu_2d,write_visu_3d
#if defined(_TURB_FORCING)
  use mod_output    , only: budget
#endif
  use mod_param     , only: lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,small,is_wallturb, &
                            cbcvel,bcvel,cbcpre,bcpre, &
                            icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                            nstep,time_max,tw_max,stop_type,restart, &
                            cfl,    &
                            constant_dt,dt_input, &
                            inivel, &
                            itot,jtot,ktot,dims_in, &
                            nthreadsmax, &
                            gr, &
                            is_outflow,no_outflow,is_forced,bforce, &
                            rho1,rho2,rho0,mu1,mu2,cbcvof,bcvof,late_init,i_late_init, &
                            rho0, &
                            n,ng,l,dl,dli, &
                            read_input, &
                            dpdl
  !
  use mod_sanity    , only: test_sanity
#if defined(_OPENACC)
  use mod_solver    , only: solver
#else
  use mod_solver_cpu, only: solver_cpu
#endif
  use mod_types
  use mod_vof       , only: initvof,advvof,update_vof,update_property,alloc_vof_var
  use profiler
!@cuf use mod_common_mpi, only: mydev
!@cuf use cudafor
  !
  !$ use omp_lib
  !
  implicit none
  !
  ! Variables declaration
  !  note: --> first decleare arrays, then the other variables;
  !        --> order of declaration: type, dimension, allocatable;
  !        --> remember to deallocate at the end.
  !
  real(rp), dimension(:,:,:)  , allocatable :: u,v,w,p,up,vp,wp,pp
  real(rp), dimension(:,:,:)  , allocatable :: pold
  !
  real(rp), dimension(:,:,:)  , allocatable :: psi,mu,rho,kappa   
  real(rp), dimension(:,:,:,:), allocatable :: cur_t                 
  real(rp), dimension(:,:,:,:), allocatable :: nor                   
  real(rp), dimension(:,:,:)  , allocatable :: d_thinc
  !
  real(rp), dimension(:), allocatable :: dzc, dzf, zc, zf, dzci, dzfi
  real(rp), dimension(:), allocatable :: dzc_g, dzf_g, zc_g, zf_g, dzci_g, dzfi_g
  !
  real(rp), dimension(:,:,:), pointer :: dudtrko,dvdtrko,dwdtrko
  !
  type(C_PTR), dimension(2,2) :: arrplanp
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:)   :: ap,bp,cp
  real(rp) :: normfftp
  ! 
  real(rp), allocatable, dimension(:,:,:) :: rhsbp_x, rhsbp_y, rhsbp_z
  !
  real(rp), dimension(3) :: f
  real(rp) :: dt,dti,dtmax,time,dtrk,dtrki,divtot,divmax
  real(rp) :: dto 
#if defined(_TIMING)
  real(rp) :: dt12,dt12av,dt12min,dt12max
#endif
  !
  integer, dimension(3)   :: halo_u,halo_up,halo_p,halo_d,halo_v,halo_t
  integer, dimension(3)   :: dims
  integer, dimension(3,3) :: dims_xyz
  integer  :: nh_d,nh_u,nh_up,nh_p,nh_v

  integer  :: i,j,k,im,ip,jm,jp,km,kp
  integer  :: irk,istep
  character(len=100) :: datadir
  character(len=1)   :: action_load
  logical  :: is_first_vel,is_data
  !
  real(rp) :: meanvel,meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: dpdl_c
  real(rp), dimension(10) :: var
  ! 
  character(len=9) :: fldnum
  real(rp) :: f1d,f2d,f3d
  real(rp) :: twi,tw
  integer :: kk
  logical :: is_done,kill
  real(rp) :: rho0i
  !

  !@cuf integer :: istat
  !@cuf integer(kind=cuda_count_kind) :: freeMem, totMem
  !@cuf attributes(managed) :: pold, kappa, mu, rho, psi, d_thinc, cur_t, nor
  !@cuf attributes(managed) :: u, v, w, p, up, vp, wp, pp 
  !@cuf attributes(managed) :: dzc, dzf, zc, zf, dzci, dzfi
  !@cuf attributes(managed) :: lambdaxyp, ap, bp, cp, rhsbp_x, rhsbp_y, rhsbp_z
  !@cuf attributes(managed) :: dudtrko, dvdtrko, dwdtrko
  !if we don't use dropcheck.f90 we can comment the next line  
  real(rp) :: xd,yd,zd,ut,vt,wt,zcd,ycd,xcd,vol
  real(rp) :: vof_min,vof_max,vol_p1
  real(rp) :: aux
  integer  :: ii,jj
  integer :: n1, n2, n3
  integer :: ng1, ng2, ng3
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  !
  call profiler_init()
  !
  ! read parameter file
  !
  call read_input(myid)
  !
  ! create data folder and subfolders for post-processing, if they do not exist.
  !
  inquire(file='data/',exist=is_data)
  if(.not.is_data.and.myid.eq.0) call execute_command_line('mkdir -p data')
  datadir = 'data/'
  !
  ! initialize MPI/OpenMP
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(.false.,ng,cbcpre,dims_in,dims_xyz,dims,n)
  !
  n1 = n(1)
  n2 = n(2)
  n3 = n(3)
  ng1 = ng(1)
  ng2 = ng(2)
  ng3 = ng(3)
  !
  twi = MPI_WTIME()
  !
  ! halo calculation (for the variables we generally 
  !                   do enlarge the stencil support)
  !
  nh_u = 1
  nh_d = nh_u 
  !
  ! allocate memory
  !
  allocate(p(0:n(1)+1,0:n(2)+1,0:n(3)+1) , &
           u(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           v(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           w(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           up(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           vp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           wp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           pp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(dudtrko(n(1),n(2),n(3)), &
           dvdtrko(n(1),n(2),n(3)), &
           dwdtrko(n(1),n(2),n(3)))
  allocate(pold(   0:n(1)+1,0:n(2)+1,0:n(3)+1)   , &
           psi(    0:n(1)+1,0:n(2)+1,0:n(3)+1)   , &
           kappa(  0:n(1)+1,0:n(2)+1,0:n(3)+1)   , &
           mu(     0:n(1)+1,0:n(2)+1,0:n(3)+1)   , &
           rho(    0:n(1)+1,0:n(2)+1,0:n(3)+1)   , &
           cur_t(  0:n(1)+1,0:n(2)+1,0:n(3)+1,6) , &
           nor(    0:n(1)+1,0:n(2)+1,0:n(3)+1,3) , &
           d_thinc(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(lambdaxyp(n_z(1),n_z(2)))
  allocate(ap(n_z(3)),bp(n_z(3)),cp(n_z(3)))
  allocate(dzc( 1-nh_d:n(3)+nh_d), &
           dzf( 1-nh_d:n(3)+nh_d), &
           zc(  1-nh_d:n(3)+nh_d), &
           zf(  1-nh_d:n(3)+nh_d), &
           dzci(1-nh_d:n(3)+nh_d), &
           dzfi(1-nh_d:n(3)+nh_d))
  allocate(dzc_g( 1-nh_d:ng(3)+nh_d), &
           dzf_g( 1-nh_d:ng(3)+nh_d), &
           zc_g(  1-nh_d:ng(3)+nh_d), &
           zf_g(  1-nh_d:ng(3)+nh_d), &
           dzci_g(1-nh_d:ng(3)+nh_d), &
           dzfi_g(1-nh_d:ng(3)+nh_d))
  allocate(rhsbp_x(n(2),n(3),0:1), &
           rhsbp_y(n(1),n(3),0:1), &
           rhsbp_z(n(1),n(2),0:1))
  !
  ! prefetching of the variables (TODO: remember to add the one of x-pencil!)
  !
  !@cuf istat = cudaMemAdvise(u, size(u), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(v, size(v), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(w, size(w), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(p, size(p), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(up, size(up), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(vp, size(vp), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(wp, size(wp), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(pp, size(pp), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dudtrko, size(dudtrko), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dvdtrko, size(dvdtrko), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dwdtrko, size(dwdtrko), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(rhsbp_x, size(rhsbp_x), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(rhsbp_y, size(rhsbp_y), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(rhsbp_z, size(rhsbp_z), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(lambdaxyp, size(lambdaxyp), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(zc, size(zc), cudaMemAdviseSetReadMostly, 0)
  !@cuf istat = cudaMemAdvise(zf, size(zf), cudaMemAdviseSetReadMostly, 0)
  !@cuf istat = cudaMemAdvise(dzc, size(dzc), cudaMemAdviseSetReadMostly, 0)
  !@cuf istat = cudaMemAdvise(dzf, size(dzf), cudaMemAdviseSetReadMostly, 0)
  !@cuf istat = cudaMemAdvise(dzci, size(dzci), cudaMemAdviseSetReadMostly, 0)
  !@cuf istat = cudaMemAdvise(dzfi, size(dzfi), cudaMemAdviseSetReadMostly, 0)
  !
#if defined(_USE_VOF)
  !@cuf istat = cudaMemAdvise(mu, size(mu), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(pold, size(pold), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(rho, size(rho), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(kappa, size(kappa), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(psi, size(psi), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(d_thinc, size(d_thinc), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(nor, size(nor), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(cur_t, size(cur_t), cudaMemAdviseSetPreferredLocation, mydev)
#endif
  !
#endif
  !
  if(myid.eq.0) print*, '*******************************'
  if(myid.eq.0) print*, '*** Beginning of simulation ***'
  if(myid.eq.0) print*, '*******************************'
#if defined(_OPENACC)
  if(myid.eq.0) then
    print*, ' GPU accelerated version, grid size:', n(1)*dims(1), n(2)*dims(2), n(3)*dims(3)
  endif
#endif
  !
  ! halo calculation (for the variables we generally 
  !                   do not enlarge the stencil support)
  !
  nh_p  = abs(lbound(p  ,1))+1
  nh_up = abs(lbound(up ,1))+1
  nh_v  = abs(lbound(psi,1))+1
  !
#if defined(_OPENACC)
  !
  ! Allocate buffers for halo communications (GPU-only)
  !
  call alloc_buf(n,nh_d)
  !
#else
  ! 
  ! halo generation using MPI derivate datatypes (CPU-only)
  !
  call halo_gen(n,nh_u ,halo_u )
  call halo_gen(n,nh_p ,halo_p )
  call halo_gen(n,nh_up,halo_up)
  call halo_gen(n,nh_v ,halo_v )
  call halo_gen(n,nh_d ,halo_d )
  !
#endif
  !
  ! initialize the grid (using global variables along z)
  !
  call initgrid(inivel,ng(3),gr,lz,nh_d,dzc_g,dzf_g,zc_g,zf_g) 
  !
  if(myid.eq.0) then
    open(99,file=trim(datadir)//'grid.bin',access='direct',recl=4*ng(3)*sizeof(1._rp))
    write(99,rec=1) dzc_g(1:ng(3)),dzf_g(1:ng(3)),zc_g(1:ng(3)),zf_g(1:ng(3))
    close(99)
    open(99,file=trim(datadir)//'grid.out')
    do kk=1-nh_d,ng(3)+nh_d
      write(99,'(5E15.7)') 1._rp*kk,zf_g(kk),zc_g(kk),dzf_g(kk),dzc_g(kk)
    enddo
    close(99)
    open(99,file=trim(datadir)//'geometry.out')
      write(99,*) ng(1),ng(2),ng(3)
      write(99,*) l(1),l(2),l(3)
    close(99)
  endif
  !
  do k=1-nh_d,ng3+nh_d
    dzfi_g(k) = 1._rp/dzf_g(k)
    dzci_g(k) = 1._rp/dzc_g(k)
  enddo
  !
  ! compute the spacing along z in local coordinates
  !
  do k=1-nh_d,n3+nh_d
    kk      = k + ijk_start(3)
    zc(k)   = zc_g(kk)
    zf(k)   = zf_g(kk)
    dzf(k)  = dzf_g(kk)
    dzc(k)  = dzc_g(kk)
    dzfi(k) = 1._rp/dzf(k)
    dzci(k) = 1._rp/dzc(k)
  enddo
  !
  !@cuf istat = cudaMemPrefetchAsync(dzci, size(dzci), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(dzfi, size(dzfi), mydev, 0)
  !
  ! test input files before proceeding with the calculation
  !
  call test_sanity(ng,n,dims_xyz(:,3),ipencil,nh_d,nh_u,nh_p,halo_d,halo_u,halo_p,stop_type, &
                   cbcvel,cbcpre,bcvel,bcpre,is_outflow,is_forced,dli,dzci_g,dzfi_g)
  !
  if(.not.restart) then
    !
    istep = 0
    time  = 0._rp
    is_first_vel = .true.
    !
    if(.not.late_init) then 
      !
      ! Initialize VoF 
      !
      call initvof(n,dli,psi)
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,psi)
      call alloc_vof_var(n(1),n(2),n(3))
    else
      !$acc kernels
      psi(:,:,:) = 0._rp 
      !$acc end kernels
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,psi)
    endif
    !
    ! [NOTE] zc is managed, lz is scalar. mmmhhh
    call initflow(inivel,n(1),n(2),n(3),dims,nh_d,nh_u,nh_p,zc/lz,dzc/lz,dzf/lz,u,v,w,p)
    !
    ! set to zeros the rhs of momentum equation 
    ! (only for the first time-step, not for the restarting)
    ! 
    !$acc kernels
    do k=1,n3
      do j=1,n2
        do i=1,n1
          dudtrko(i,j,k) = 0._rp
          dvdtrko(i,j,k) = 0._rp
          dwdtrko(i,j,k) = 0._rp
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    if(myid.eq.0) print*, '*** Initial condition succesfully set ***'
    !
  else
    !
    action_load = 'r'
    call load(action_load,trim(datadir)//'fldu.bin',n,      u(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldv.bin',n,      v(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldw.bin',n,      w(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'flddu.bin',n,dudtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'flddv.bin',n,dvdtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'flddw.bin',n,dwdtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldp.bin',n,      p(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldpsi.bin',n,    psi(1:n(1),1:n(2),1:n(3)))
    call load_scalar(action_load,trim(datadir)//'scalar.out',time,istep,dto)
    !
    if(myid.eq.0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
    !
    is_first_vel = .false.
    !
    call alloc_vof_var(n(1),n(2),n(3)) ! to put here (important for the restarting)
  endif
  !
  !@cuf istat = cudaMemPrefetchAsync(dudtrko, size(dudtrko), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(dvdtrko, size(dvdtrko), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(dwdtrko, size(dwdtrko), mydev, 0)
  !
  ! set boundary conditions on the initial/loaded fields
  !
  call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
  call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,psi)
  !
  ! for the first time-step and the restarting we use a 0th order extrapolation 
  ! in time-splitting of the pressure equation
  !
  !$acc kernels 
  do k=1,n3
    do j=1,n2
      do i=1,n1
        pold(i,j,k) = p(i,j,k)
      end do
    end do
  end do
  !$acc end kernels
  !
  call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,pold)
  !
  ! update the quantities derived from VoF using the lastest available VoF field
  !
  call update_vof(n,dli,nh_d,dzc,dzf,nh_v,halo_v,psi,nor,cur_t,kappa,d_thinc)
  call update_property(n,(/rho1,rho2/),psi,rho)
  call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)
  call update_property(n,(/mu1 ,mu2 /),psi,mu ) 
  call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,mu)

  !
  ! post-process and write initial condition
#if !defined(__BENCHMARK_NO_IO)
  !
  call profiler_start("OUT:initial", tag = .true., tag_color = COLOR_WHITE)
  !
  ! Prefetching back pre IO
  !@cuf istat = cudaMemPrefetchAsync(u, size(u), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(v, size(v), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(w, size(w), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(p, size(p), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(psi, size(psi), cudaCpuDeviceId, 0)
  !
  write(fldnum,'(i9.9)') istep
  include 'out1d.h90'
  include 'out2d.h90'
  include 'out3d.h90'
  !
  ! Prefetching post IO
  !@cuf istat = cudaMemPrefetchAsync(u, size(u), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(v, size(v), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(w, size(w), mydev, 0)
  !
  call profiler_stop("OUT:initial")
#endif
  !
  ! compute an initial time-step
  !
  if(.not.constant_dt) then
    call chkdt(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
    dt = cfl*dtmax
  else
    if(myid.eq.0) print*, 'the simulation is run at constant time-step'
    dtmax = dt_input
    dt    = dtmax
  endif
  if(istep.eq.0) dto = dt
  if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ', dt
  dti  = 1._rp/dt
  kill = .false.
  !
  ! Prefetching post IO
  !@cuf istat = cudaMemPrefetchAsync(psi, size(psi), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(p, size(p), mydev, 0)
  !
  ! preliminary checks
  ! 
  vof_min = minval(psi(1:n1,1:n2,1:n3))
  call mpi_allreduce(MPI_IN_PLACE,vof_min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  vof_max = maxval(psi(1:n1,1:n2,1:n3))
  call mpi_allreduce(MPI_IN_PLACE,vof_max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  vol_p1  = sum(psi(1:n1,1:n2,1:n3))*dl(1)*dl(2)*dl(3)
  call mpi_allreduce(MPI_IN_PLACE,vol_p1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  var(:) = 0._rp
  var(1) = 1._rp*istep
  var(2) = dt
  var(3) = time
  var(4) = vof_min
  var(5) = vof_max
  var(6) = vol_p1
  call out0d(trim(datadir)//'vof_info.out',6,var)
  !
  ! initialize Poisson solver
  ! and deallocate global arrays (not needed anymore) 
  !
  call initsolver(n,dims,dims_xyz(:,3),dli,nh_d,dzci_g,dzfi_g,cbcpre,bcpre(:,:),(/'c','c','c'/),lambdaxyp, & 
                  ap,bp,cp,arrplanp,normfftp,rhsbp_x,rhsbp_y,rhsbp_z)
  !@cuf istat = cudaMemPrefetchAsync(rhsbp_x, size(rhsbp_x), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(rhsbp_y, size(rhsbp_y), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(rhsbp_z, size(rhsbp_z), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(lambdaxyp, size(lambdaxyp), mydev, 0)
  deallocate(dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g)
  !
  ! main loop
  !
  if(myid.eq.0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  !
  !@cuf istat = cudaMemGetInfo( freeMem, totMem )
  !@cuf if(myid.eq.0) print*, 'Used memory = ', totMem - freeMem
  !
  do while(.not.is_done)
    !
#if defined(_TIMING)
    dt12  = MPI_WTIME()
#endif
    !
    istep = istep + 1
    !
    call profiler_start("STEP", tag = .true., tag_color = COLOR_WHITE)
    !
    time  = time + dt
    !
    if(myid.eq.0) print*, 'Timestep #', istep, 'Time = ', time
    !
    if(late_init.and.(istep.eq.i_late_init)) then
      !
      call initvof(n,dli,psi)
      call alloc_vof_var(n(1),n(2),n(3))
      call update_vof(n,dli,nh_d,dzc,dzf,nh_v,halo_v,psi,nor,cur_t,kappa,d_thinc)
      !
      call update_property(n,(/mu1 ,mu2 /),psi,mu ) 
      call update_property(n,(/rho1,rho2/),psi,rho)
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,mu )
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)
      call chkdt(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
      !
    endif
    !
    ! 0. save quantities from the previous time-step
    !
#if defined(_OPENACC)
    !$acc parallel loop collapse(3) present(p, pp, pold)
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(p,pold,pp)
#endif
    do k=1,n3
      do j=1,n2
        do i=1,n1
          pp(i,j,k)   = (1._rp+(dt/dto))*p(i,j,k) - (dt/dto)*pold(i,j,k)
          pold(i,j,k) = p(i,j,k)
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end parallel loop 
#else
    !$OMP END PARALLEL DO
#endif
    call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,pp  )
    call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,pold)
    !
    ! 1. VoF advection and properties update --> vof^(n+1) 
    !
    if((.not.late_init).or. &                     ! Normal run (initialized at i=0)
       (late_init.and.istep.ge.i_late_init)) then ! in case we want first to run without advecting the VoF
      !
      call profiler_start("VOF", tag = .true., tag_color = COLOR_YELLOW)
      !
      call advvof(n,dli,dt,halo_v,nh_d,nh_u,dzc,dzf,u,v,w,psi,nor,cur_t,kappa,d_thinc)
      !
      call update_property(n,(/rho1,rho2/),psi,rho)           
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)
      call update_property(n,(/mu1 ,mu2 /),psi,mu )           
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,mu )
      !
      call profiler_stop("VOF")
      !
    endif
    !
#if !defined(_VOF_DBG)
    !
    ! 2. two-fluid Navier-Stokes --> (u,v,w,p)^(n+1)
    !
    rho0i = 1._rp/rho0
    dpdl_c(1:3) = 0._rp
    !
    !
    call profiler_start("RK", tag = .true., tag_color = COLOR_RED)
    !
    call rk(is_first_vel,n(1),n(2),n(3),dims,dxi,dyi,dzi,nh_d,nh_u,dt,dto,rho0i,dzci,dzfi,u,v,w,p,pp, &
            kappa,psi,mu,rho,dudtrko,dvdtrko,dwdtrko,up,vp,wp,f)
    !
    call profiler_stop("RK")
    !
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_up,halo_up,no_outflow,dl,dzc,dzf,up,vp,wp)
    dpdl_c(1:3) = dpdl_c(1:3) + f(1:3) + dpdl(1:3)*dt
    !
    call fillps(n(1),n(2),n(3),nh_d,nh_up,dxi,dyi,dzi,dzci,dzfi,dti,rho0,up,vp,wp,p)
    !
    call updt_rhs_b(n(1),n(2),n(3),(/'c','c','c'/),cbcpre,nh_p,rhsbp_x,rhsbp_y,rhsbp_z,p)
    !
    call profiler_start("SOLVER", tag = .true., tag_color = COLOR_GREEN)
    !
#if defined(_OPENACC)
    call solver(n_z,dims_xyz(:,3),arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,:),(/'c','c','c'/),p)
#else
    call solver_cpu(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),p)
#endif
    !
    call profiler_stop("SOLVER")
    !
    call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
    !
    call profiler_start("CORREC")
    !
    call correc(n(1),n(2),n(3),nh_d,nh_u,dxi,dyi,dzi,dzci,dt,rho0,p,up,vp,wp,rho,u,v,w)
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf,u,v,w)
    !
#if defined(_OPENACC)
    !$acc parallel loop collapse(3)
#else
    !$OMP WORKSHARE
#endif
    do k=1,n3
      do j=1,n2
        do i=1,n1
          p(i,j,k) = pold(i,j,k) + p(i,j,k)
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end parallel loop
#else
    !$OMP END WORKSHARE
#endif
    call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
    !
    call profiler_stop("CORREC")
    ! 
    !
#endif ! [NOTE] This endif is misaligned compared to HT, please check
    ! 
    write(fldnum,'(i9.9)') istep
    !
    if(.not.late_init.and.istep.ge.i_late_init) then
      include 'dropcheck.h90'
    endif
    !
#if defined(_TURB_FORCING)
    call budget(datadir,n(1),n(2),n(3),ng(1),ng(2),ng(3),rho1,rho2,dx,dy,dz,nh_u, &
                u,v,w,rho,psi,time,istep)
#endif
    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep.ge.nstep   ) is_done = is_done.or..true.
    endif
    if(stop_type(2)) then ! maximum simulation time reached
      if(time .ge.time_max) is_done = is_done.or..true.
    endif
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600._rp
      if(tw   .ge.tw_max  ) is_done = is_done.or..true.
    endif
    !
    ! check time-step size and velocity divergence
    !
    dto = dt  
    if(mod(istep,icheck).eq.0) then
      !
      if(myid.eq.0) print*, 'Checking stability and divergence...'
      !
      if(.not.constant_dt) then
        if(myid.eq.0) print*, 'updating the time-step size ...'
        call chkdt(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
        dt = cfl*dtmax
      else
        dtmax = dt_input
        dt    = dtmax
      endif
      !
      if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ',dt
      if(dtmax.lt.small) then
        if(myid.eq.0) print*, 'ERROR: timestep is too small.'
        if(myid.eq.0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      endif
      dti = 1._rp/dt
      !
      if(myid.eq.0) print*, 'checking the velocity divergence ...'
      call chkdiv(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzfi,u,v,w,divtot,divmax)
      if(divmax.gt.small.or.divtot.ne.divtot) then
        if(myid.eq.0) print*, 'ERROR: maximum divergence is too large.'
        if(myid.eq.0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      endif
      !
    endif
    !
    call profiler_stop("STEP")
    !
    ! output routines below
    !
#if !defined(__BENCHMARK_NO_IO)
    if(mod(istep,iout0d).eq.0) then
      !
      call profiler_start("OUT:iout0d", tag = .true., tag_color = COLOR_WHITE)
      !
      var(:) = 0._rp
      var(1) = 1._rp*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
      !
      vof_min = minval(psi(1:n(1),1:n(2),1:n(3)))
      call mpi_allreduce(MPI_IN_PLACE,vof_min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      vof_max = maxval(psi(1:n(1),1:n(2),1:n(3)))
      call mpi_allreduce(MPI_IN_PLACE,vof_max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      vol_p1  = sum(psi(1:n(1),1:n(2),1:n(3)))*dl(1)*dl(2)*dl(3)
      call mpi_allreduce(MPI_IN_PLACE,vol_p1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !
      var(:) = 0._rp
      var(1) = 1._rp*istep
      var(2) = dt
      var(3) = time
      var(4) = vof_min
      var(5) = vof_max
      var(6) = vol_p1
      call out0d(trim(datadir)//'vof_info.out',6,var)
      !
      if(any(is_forced(:)).or.(abs(sum(dpdl_c)).ge.small)) then ! control strategy: constant flow-rate
        call chkmean(n(1),n(2),n(3),dims,nh_d,nh_up,1._rp/(dzfi*lz),u,meanvelu)
        call chkmean(n(1),n(2),n(3),dims,nh_d,nh_up,1._rp/(dzfi*lz),v,meanvelv)
        call chkmean(n(1),n(2),n(3),dims,nh_d,nh_up,1._rp/(dzfi*lz),w,meanvelw)
        var(:)   = 0._rp
        var(1)   = time
        var(2)   = 1._rp*istep
        var(3:5) = -dpdl_c(1:3)*dti*rho2
        var(6)   = sqrt(maxval(abs(var(3:5)))*lz*0.5_rp)
        var(7)   = rho2*var(6)*0.5_rp*lz/mu2 ! Re_tau
        var(8)   = meanvelu
        var(9)   = meanvelv
        var(10)  = meanvelw
        call out0d(trim(datadir)//'forcing.out',10,var)
      endif 
      !
      call profiler_stop("OUT:iout0d")
      !
    endif
#endif
    !
    if(mod(istep,iout1d).eq.0) then
      call profiler_start("OUT:iout1d", tag = .true., tag_color = COLOR_WHITE)
      ! [TODO] Prefetching GPU->CPU (and viceversa)
      include 'out1d.h90'
      call profiler_stop("OUT:iout1d")
    endif
    if(mod(istep,iout2d).eq.0) then
      call profiler_start("OUT:iout2d", tag = .true., tag_color = COLOR_WHITE)
      ! [TODO] Prefetching GPU->CPU (and viceversa)
      include 'out2d.h90'
      call profiler_stop("OUT:iout2d")
    endif
    if(mod(istep,iout3d).eq.0) then
      call profiler_start("OUT:iout3d", tag = .true., tag_color = COLOR_WHITE)
      ! [TODO] Prefetching GPU->CPU (and viceversa)
      include 'out3d.h90'
      call profiler_stop("OUT:iout3d")
    endif
    !
    if ( (mod(istep,isave) .eq. 0) .or. is_done .and. (.not.kill) ) then
      !
      call profiler_start("OUT:isave", tag = .true., tag_color = COLOR_WHITE)
      !
      action_load = 'w'
      !
      !@cuf istat = cudaMemPrefetchAsync(u, size(u), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(v, size(v), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(w, size(w), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(dudtrko, size(dudtrko), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(dvdtrko, size(dvdtrko), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(dwdtrko, size(dwdtrko), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(p, size(p), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(psi, size(psi), cudaCpuDeviceId, 0)
      !
      inquire(file='data/fldu.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldu.bin   data/fldu_old.bin')
      call load(action_load,trim(datadir)//'fldu.bin',n,      u(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(u, size(u), mydev, 0)
      !
      inquire(file='data/fldv.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldv.bin   data/fldv_old.bin')
      call load(action_load,trim(datadir)//'fldv.bin',n,      v(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(v, size(v), mydev, 0)
      !
      inquire(file='data/fldw.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldw.bin   data/fldw_old.bin')
      call load(action_load,trim(datadir)//'fldw.bin',n,      w(1:n(1),1:n(2),1:n(3)))   
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(w, size(w), mydev, 0)
      !
      inquire(file='data/flddu.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddu.bin  data/flddu_old.bin')
      call load(action_load,trim(datadir)//'flddu.bin',n,dudtrko(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(dudtrko, size(dudtrko), mydev, 0)
      !
      inquire(file='data/flddv.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddv.bin  data/flddv_old.bin')
      call load(action_load,trim(datadir)//'flddv.bin',n,dvdtrko(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(dvdtrko, size(dvdtrko), mydev, 0)
      !
      inquire(file='data/flddw.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddw.bin  data/flddw_old.bin')
      call load(action_load,trim(datadir)//'flddw.bin',n,dwdtrko(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(dwdtrko, size(dwdtrko), mydev, 0)
      !
      inquire(file='data/flddp.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldp.bin   data/fldp_old.bin')
      call load(action_load,trim(datadir)//'fldp.bin',n,      p(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(p, size(p), mydev, 0)
      !
      inquire(file='data/fldpsi.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldpsi.bin data/fldpsi_old.bin')
      call load(action_load,trim(datadir)//'fldpsi.bin',n,    psi(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(psi, size(psi), mydev, 0)
      !
      inquire(file='data/scalar.out', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/scalar.out data/scalar_old.out')
      call load_scalar(action_load,trim(datadir)//'scalar.out',time,istep,dto)
      !
      if(myid.eq.0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
      !
      call profiler_stop("OUT:isave")
      !
    endif
    !
#if defined(_TIMING)
    dt12 = MPI_WTIME()-dt12
    call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    !
    var(:) = 0._rp
    var(1) = 1._rp*istep
    var(2) = time
    var(3) = dt12av/(1._rp*product(dims))
    var(4) = dt12min
    var(5) = dt12max
    call out0d(trim(datadir)//'performance.out',5,var)
#endif
    !
  enddo
  !
  ! clear ffts
  !
  call fftend(arrplanp)
  !
  ! deallocate memory
  !
  deallocate(u,v,w,p,up,vp,wp,pp, &
             pold, &
             psi,kappa,mu,rho, &     
             cur_t, &                   
             nor, &                 
             d_thinc, &
             dzc,dzf,zc,zf,dzci,dzfi, &
             dudtrko,dvdtrko,dwdtrko)
  !
  deallocate(lambdaxyp)
  deallocate(ap,bp,cp)
  deallocate(rhsbp_x,rhsbp_y,rhsbp_z)
  !
  if(myid.eq.0.and.(.not.kill)) print*, '*** Fim ***'
  !
  call profiler_report()
  !
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
  call exit
  !
end program flutas
