!
! SPDX-License-Identifier: MIT
!
!--------------------------------------------------------------------------------------------------------------------
! FluTAS -- Fluid Transport Accelerated Solver                            
!                                                                                       
!  a. Menu Title: FluTAS_single_phase;                                                         
!  b. Feautures of FluTAS_single_phase:                                                        
!      --> single-phase solver, optionally with heat transfer;
!      --> Oberbeck-Boussinesq approximation; 
!      --> momentum advanced with Adams-Bashforth (explicit diffusion);                 
!      --> pressure equation solved with FFT-based direct solver.                                       
!--------------------------------------------------------------------------------------------------------------------
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
  use mod_chkdt     , only: chkdt_sp
  use mod_common_mpi, only: myid,ierr,comm_cart,n_z,ijk_start,ipencil
  use mod_correc    , only: correc
  use mod_debug     , only: cmpt_mean
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
  use mod_rk        , only: rk,cmpt_time_factors
  use mod_output    , only: out0d,out1d,out2d,out3d,write_visu_2d,write_visu_3d
  use mod_param     , only: lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,small, &
                            cbcvel,bcvel,cbcpre,bcpre, &
                            icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                            nstep,time_max,tw_max,stop_type,restart, &
                            cfl,     &
                            constant_dt,dt_input, &
                            inivel,  &
                            itot,jtot,ktot,dims_in, &
                            nthreadsmax,gr, &
                            is_outflow,no_outflow,is_forced,rho_sp,mu_sp, &
#if defined(_DO_POSTPROC)
                            do_avg,time_deltai,&
                            avg_dir,&
                            do_favre,&
                            do_wall,wall_deltai,&
#endif
#if defined(_HEAT_TRANSFER)
                            initmp, &
                            cbctmp,bctmp, &
                            kappa_sp,cp_sp, &
#endif
                            n,ng,l,dl,dli, &
                            bulk_ftype,rkcoeff, &
                            time_scheme,space_scheme_mom,n_stage, &
                            read_input
  !
  use mod_source    , only: bulk_forcing_src,grav_sp_src,pres_sp_src
#if defined(_TURB_FORCING)
  use mod_source    , only: forc_src 
#endif
  use mod_sanity    , only: test_sanity
#if defined(_HEAT_TRANSFER)
  use mod_initflow  , only: inittmp
  use mod_rks       , only: rk_sca
#endif
#if defined(_OPENACC)
  use mod_solver_gpu, only: solver_gpu
#else
  use mod_solver_cpu, only: solver_cpu
#endif
  use mod_types
#if defined(_DO_POSTPROC)
  use mod_post      , only: time_sp_avg,wall_avg,compute_vorticity,mixed_variables
#endif
  use profiler
  !@cuf use mod_common_mpi, only: mydev
  !@cuf use cudafor
  !
  !$ use omp_lib
  !
  implicit none
  !
  ! Variables declaration
  !  note: --> first declare arrays, then the other variables;
  !        --> order of declaration: type, dimension, allocatable;
  !
  real(rp), dimension(:,:,:), allocatable :: u,v,w,p
  real(rp), dimension(:,:,:), allocatable :: pold
  real(rp), dimension(:,:,:), allocatable :: mu,rho
  real(rp), dimension(:,:,:), pointer :: dudtrko,dvdtrko,dwdtrko
  real(rp), dimension(:), allocatable :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi
  real(rp), dimension(:), allocatable :: dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g
  !
#if defined(_HEAT_TRANSFER)
  real(rp), dimension(:,:,:)  , allocatable :: tmp
  real(rp), dimension(:,:,:)  , allocatable :: dtmpdtrk
  real(rp), dimension(:,:,:)  , allocatable :: dtmpdtrkold
  real(rp), dimension(:,:,:)  , allocatable :: ka,cpp
#endif
  !
  type(C_PTR), dimension(2,2) :: arrplanp
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:)   :: ap,bp,cp
  real(rp) :: normfftp
  ! 
  real(rp), allocatable, dimension(:,:,:) :: rhsbp_x, rhsbp_y, rhsbp_z
  !
  real(rp), dimension(3) :: f
  real(rp) :: dt,dto,dti,dtmax,time,dtrk,dtrki,divtot,divmax, &
              f_t1,f_t2,f_t12,f_t12_i!,f_t12_o
#if defined(_TIMING)
  real(rp) :: dt12,dt12av,dt12min,dt12max
#endif
  !
  integer, dimension(3)   :: halo_u,halo_p,halo_d
  integer, dimension(3)   :: dims
  integer, dimension(3,3) :: dims_xyz
  integer  :: nh_d,nh_u,nh_p
#if defined(_HEAT_TRANSFER)
  integer, dimension(3) :: halo_t
  integer  :: nh_t
#endif
  !
  integer  :: i,j,k,im,ip,jm,jp,km,kp
  integer  :: irk_ss,istep
  character(len=1) :: action_load
  logical  :: is_data
  !
  real(rp) :: meanvel,meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: dpdl_c
  real(rp), dimension(20) :: var
  ! 
  character(len=9)   :: fldnum
  character(len=100) :: datadir,restart_dir
  real(rp) :: f1d,f2d,f3d
  real(rp) :: twi,tw
  integer :: kk
  logical :: is_done,kill
  !
  !@cuf integer :: istat
  !@cuf integer(kind=cuda_count_kind) :: freeMem, totMem
  !@cuf attributes(managed) :: u, v, w, p, mu, rho, pold 
  !@cuf attributes(managed) :: dudtrko, dvdtrko, dwdtrko
  !@cuf attributes(managed) :: dzc  , dzf  , dzci, dzfi, zc, zf, lambdaxyp, ap, bp, cp, rhsbp_x, rhsbp_y, rhsbp_z
  !@cuf attributes(managed) :: zc_g, zf_g
#if defined(_HEAT_TRANSFER) 
  !@cuf attributes(managed) :: tmp, ka, cpp, dtmpdtrk, dtmpdtrkold
#endif
  !
#if defined(_DO_POSTPROC)
  real(rp), allocatable, dimension(:,:,:) :: uv,vw,wu 
  real(rp), allocatable, dimension(:,:,:) :: vorx 
  real(rp), allocatable, dimension(:) :: u_avg_g,v_avg_g,w_avg_g,u_sqr_g,v_sqr_g,w_sqr_g 
  real(rp), allocatable, dimension(:) :: uv_avg_g,vw_avg_g,wu_avg_g,uv_sqr_g,vw_sqr_g,wu_sqr_g 
  real(rp), allocatable, dimension(:) :: vorx_avg_g,vorx_sqr_g 
  real(rp) :: u_vol_avg_g,u_vol_sqr_g,v_vol_avg_g,v_vol_sqr_g,w_vol_avg_g,w_vol_sqr_g 
  real(rp) :: uv_vol_avg_g,uv_vol_sqr_g,vw_vol_avg_g,vw_vol_sqr_g,wu_vol_avg_g,wu_vol_sqr_g 
  real(rp) :: vorx_vol_avg_g,vorx_vol_sqr_g 
  real(rp) :: vory_vol_avg_g,vory_vol_sqr_g 
  real(rp) :: vorz_vol_avg_g,vorz_vol_sqr_g 
  integer  :: i_av 
  real(rp), allocatable, dimension(:,:,:) :: utmp,vtmp,wtmp 
  real(rp), allocatable, dimension(:) :: tmp_avg_g,tmp_sqr_g 
  real(rp), allocatable, dimension(:) :: utmp_avg_g,vtmp_avg_g,wtmp_avg_g,utmp_sqr_g,vtmp_sqr_g,wtmp_sqr_g  
  real(rp) :: tmp_vol_avg_g,tmp_vol_sqr_g 
  real(rp) :: utmp_vol_avg_g,utmp_vol_sqr_g,vtmp_vol_avg_g,vtmp_vol_sqr_g,wtmp_vol_avg_g,wtmp_vol_sqr_g 
  !@cuf attributes(managed) :: vorx,uv,vw,wu,utmp,vtmp,wtmp
  !@cuf attributes(managed) :: u_avg_g,v_avg_g,w_avg_g,u_sqr_g,v_sqr_g,w_sqr_g 
  !@cuf attributes(managed) :: uv_avg_g,vw_avg_g,wu_avg_g,uv_sqr_g,vw_sqr_g,wu_sqr_g 
  !@cuf attributes(managed) :: vorx_avg_g,vorx_sqr_g 
  !@cuf attributes(managed) :: tmp_avg_g,tmp_sqr_g 
  !@cuf attributes(managed) :: utmp_avg_g,vtmp_avg_g,wtmp_avg_g,utmp_sqr_g,vtmp_sqr_g,wtmp_sqr_g  
#endif
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
  inquire(file='data/restart_dir/',exist=is_data)
  if(.not.is_data.and.myid.eq.0) call execute_command_line('mkdir -p data/restart_dir')
  restart_dir = 'data/restart_dir/'
  !
#if defined(_DO_POSTPROC)
  inquire(file='data/post/',exist=is_data)
  if(.not.is_data.and.myid.eq.0) call execute_command_line('mkdir -p data/post')
  inquire(file='data/post/time_averaging/',exist=is_data)
  if(do_avg.and..not.is_data.and.myid.eq.0) call execute_command_line('mkdir -p data/post/time_averaging')
  if(do_avg) i_av = 0
  inquire(file='data/post/wall/',exist=is_data)
  if(do_wall.and..not.is_data.and.myid.eq.0) call execute_command_line('mkdir -p data/post/wall')
#endif
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
  ! halo calculation
  !
  nh_p = 1
  if(    space_scheme_mom.eq.'cen') then
    nh_u = 1
  elseif(space_scheme_mom.eq.'fll') then
    nh_u = 3
  endif
#if defined(_HEAT_TRANSFER)
  nh_t = 3 ! we use WENO5 for scalar advection
#endif
  !
  nh_d = max(nh_u,nh_p) ! take the maximum of the previous ones
#if defined(_HEAT_TRANSFER)
  nh_d = max(nh_d,nh_t)
#endif
  !
  ! allocate memory
  !
  allocate(p(0:n(1)+1,0:n(2)+1,0:n(3)+1) , &
           u(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           v(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           w(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           pold(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(rho(0,0,0),mu(0,0,0))
#if defined(_HEAT_TRANSFER)
  allocate(tmp(1-nh_t:n(1)+nh_t,1-nh_t:n(2)+nh_t,1-nh_t:n(3)+nh_t), &
            ka(0,0,0),cpp(0,0,0), &
              dtmpdtrk(1:n(1),1:n(2),1:n(3)),&
           dtmpdtrkold(1:n(1),1:n(2),1:n(3)))
#endif
  allocate(dudtrko(n(1),n(2),n(3)), &
           dvdtrko(n(1),n(2),n(3)), &
           dwdtrko(n(1),n(2),n(3)))
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
  !@cuf istat = cudaMemAdvise(pold, size(pold), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(mu, size(mu), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(rho, size(rho), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dudtrko, size(dudtrko), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dvdtrko, size(dvdtrko), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dwdtrko, size(dwdtrko), cudaMemAdviseSetPreferredLocation, mydev)
#if defined(_HEAT_TRANSFER)
  !@cuf istat = cudaMemAdvise(tmp, size(tmp), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise( ka, size( ka), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(cpp, size(cpp), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dtmpdtrk, size(dtmpdtrk), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dtmpdtrkold, size(dtmpdtrkold), cudaMemAdviseSetPreferredLocation, mydev)
#endif
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
  !@cuf istat = cudaMemAdvise(zc_g , size(zc_g) , cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(zf_g , size(zf_g) , cudaMemAdviseSetPreferredLocation, mydev)
  !
#if defined(_DO_POSTPROC)
  !@cuf istat = cudaMemAdvise(vorx, size(vorx), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(  uv,   size(uv), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(  vw,   size(vw), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(  wu,   size(wu), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(utmp, size(utmp), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(vtmp, size(vtmp), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(wtmp, size(wtmp), cudaMemAdviseSetPreferredLocation, mydev)
#endif
  !
  if(myid.eq.0) print*, '***************************************************'
  if(myid.eq.0) print*, '*** Beginning of simulation (SINGLE-PHASE mode) ***'
  if(myid.eq.0) print*, '***************************************************'
  !
#if defined(_OPENACC)
  if(myid.eq.0) then
    print*, ' GPU accelerated version, grid size:', n(1)*dims(1), n(2)*dims(2), n(3)*dims(3)
  endif
#endif
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
  call halo_gen(n,nh_d ,halo_d )
#if defined(_HEAT_TRANSFER)
  call halo_gen(n,nh_t ,halo_t )
#endif
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
  !@cuf istat = cudaMemPrefetchAsync(zc_g , size( zc_g), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(zf_g , size( zf_g), mydev, 0)
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
    !
#if defined(_HEAT_TRANSFER)
    !
    ! Initialize temperature 
    !
    call inittmp(initmp,n(1),n(2),n(3),nh_t,dims,.false.,0.5d0,tmp)
    !@cuf istat = cudaMemPrefetchAsync(tmp, size(tmp), mydev, 0)
    call boundp(cbctmp,n,bctmp,nh_d,nh_t,halo_t,dl,dzc,dzf,tmp)
    !
#endif
    !
    call initflow(inivel,n(1),n(2),n(3),dims,nh_d,nh_u,nh_p,rho_sp,mu_sp,zc/lz,dzc/lz,dzf/lz,u,v,w,p)
    !
    !$acc kernels
    do k=1,n3
      do j=1,n2
        do i=1,n1
          ! 
          ! set to zeros the rhs of momentum equation 
          ! (only for the first time-step, not for the restarting)
          !
          dudtrko(i,j,k) = 0._rp
          dvdtrko(i,j,k) = 0._rp
          dwdtrko(i,j,k) = 0._rp
          !
#if defined(_HEAT_TRANSFER)
          !
          ! set to zeros the rhs of temperature equation 
          ! (only for the first time-step, not for the restarting)
          !
          dtmpdtrkold(i,j,k) = 0._rp
#endif
          !
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
    !
    call load(action_load,trim(restart_dir)//'fldu.bin',n,      u(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldv.bin',n,      v(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldw.bin',n,      w(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldp.bin',n,      p(1:n(1),1:n(2),1:n(3)))
    if(time_scheme.ne.'rk3') then ! RK3 is self-restarting
      call load(action_load,trim(restart_dir)//'flddu.bin',n,dudtrko(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(restart_dir)//'flddv.bin',n,dvdtrko(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(restart_dir)//'flddw.bin',n,dwdtrko(1:n(1),1:n(2),1:n(3)))
    endif
#if defined(_HEAT_TRANSFER)
    call load(action_load,trim(restart_dir)//'fldtmp.bin',n,    tmp(1:n(1),1:n(2),1:n(3)))
    if(time_scheme.ne.'rk3') then ! RK3 is self-restarting
      call load(action_load,trim(restart_dir)//'flddtmp.bin',n,dtmpdtrkold(1:n(1),1:n(2),1:n(3)))
    endif
#endif
    call load_scalar(action_load,trim(restart_dir)//'scalar.out',time,istep,dto)
    !
    if(myid.eq.0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
    !
  endif
  !
  !@cuf istat = cudaMemPrefetchAsync(dudtrko, size(dudtrko), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(dvdtrko, size(dvdtrko), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(dwdtrko, size(dwdtrko), mydev, 0)
#if defined(_HEAT_TRANSFER)
  !@cuf istat = cudaMemPrefetchAsync(dtmpdtrkold, size(dtmpdtrkold), mydev, 0)
#endif
  !
  ! set boundary conditions on the initial/loaded fields
  !
  call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
#if defined(_HEAT_TRANSFER)
  call boundp(cbctmp,n,bctmp,nh_d,nh_t,halo_t,dl,dzc,dzf,tmp)
#endif
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
  ! we fill these "dummy" arrays with unreasonable values
  ! to be sure they are not employed anywhere in the code
  ! (otherwise, it would crash)
  !
  rho(:,:,:) = - huge(1._rp)
  mu(:,:,:)  = - huge(1._rp)
#if defined(_HEAT_TRANSFER)
  ka(:,:,:)  = - huge(1._rp)
  cpp(:,:,:) = - huge(1._rp)
#endif
  !
  ! post-process and write initial condition
  !
#if !defined(_BENCHMARK_NO_IO)
  !
  call profiler_start("OUT:initial", tag = .true., tag_color = COLOR_WHITE)
  !
  ! Prefetching back pre IO
  !@cuf istat = cudaMemPrefetchAsync(u, size(u), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(v, size(v), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(w, size(w), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(p, size(p), cudaCpuDeviceId, 0)
#if defined(_HEAT_TRANSFER)
  !@cuf istat = cudaMemPrefetchAsync(tmp, size(tmp), cudaCpuDeviceId, 0)   
#endif
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
  !@cuf istat = cudaMemPrefetchAsync(p, size(p), mydev, 0)
#if defined(_HEAT_TRANSFER)
  !@cuf istat = cudaMemPrefetchAsync(tmp, size(tmp), mydev, 0)
#endif
  !
  call profiler_stop("OUT:initial")
#endif
  !
  ! compute an initial time-step
  !
  if(.not.constant_dt) then
    call chkdt_sp(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
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
  ! initialize Poisson solver
  ! and deallocate global arrays (not needed anymore) 
  !
  call initsolver(n,dims,dims_xyz(:,3),dli,nh_d,dzci_g,dzfi_g,cbcpre,bcpre(:,:),(/'c','c','c'/),lambdaxyp, & 
                  ap,bp,cp,arrplanp,normfftp,rhsbp_x,rhsbp_y,rhsbp_z)
  !@cuf istat = cudaMemPrefetchAsync(rhsbp_x, size(rhsbp_x), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(rhsbp_y, size(rhsbp_y), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(rhsbp_z, size(rhsbp_z), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(lambdaxyp, size(lambdaxyp), mydev, 0)
  deallocate(dzc_g,dzf_g,dzci_g,dzfi_g)
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
    if(any(is_forced(:))) dpdl_c(1:3) = 0._rp
    !
    do irk_ss = 1,n_stage ! n_stage=1 for AB2, n_stage=3 for RK3
      !
      ! 0. compute the coefficients for the time advancement
      !
      call cmpt_time_factors(time_scheme,restart,istep,irk_ss,rkcoeff,dt,dto,f_t1,f_t2,f_t12)
      f_t12_i = 1._rp/f_t12
      !
#if defined(_HEAT_TRANSFER)
      !
      ! 1. Temperature advection --> tmp^(n+1)
      !
      call profiler_start("TMP", tag = .true., tag_color = COLOR_YELLOW)
      !
      call rk_sca(f_t1,f_t2,n(1),n(2),n(3),dli(1),dli(2),dli(3),nh_d,nh_u,nh_t,rho,cpp,ka,dzci,dzfi, &
                  u,v,w,tmp,dtmpdtrk,dtmpdtrkold)
      !
      call boundp(cbctmp,n,bctmp,nh_d,nh_t,halo_t,dl,dzc,dzf,tmp)
      !
      call profiler_start("TMP", tag = .true., tag_color = COLOR_YELLOW)
      !
#endif
      !
      ! 2. Navier-Stokes solver --> (u,v,w,p)^(n+1)
      !
      call profiler_start("RK", tag = .true., tag_color = COLOR_RED)
      !
      call rk(space_scheme_mom,f_t1,f_t2,n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi, &
              u,v,w,mu,rho,dudtrko,dvdtrko,dwdtrko)
      !
      ! 2b. add the source terms
      !
      ! --> add gravity terms
      !
      call grav_sp_src(n(1),n(2),n(3),f_t12,cbcpre,dxi,dyi,dzi,nh_d,dzfi, &
#if defined(_HEAT_TRANSFER) && defined(_BOUSSINESQ)
                       nh_t,tmp, & 
#endif
                       nh_u,u,v,w)  
      !
      ! --> add pressure gradient
      !
      call pres_sp_src(n(1),n(2),n(3),f_t12,dxi,dyi,dzi,nh_d,nh_u,dzci,1._rp/rho_sp,p,u,v,w)
      !
#if defined(_TURB_FORCING)
      !
      ! --> add the forcing to sustain turbulence (triperiodic cases)
      !
      call forc_src(n(1),n(2),n(3),nh_d,nh_u,f_t12,dx,dy,dz,zc,u,v,w)
#endif
      !
      ! --> add the bulk velocity forcing 
      !     note: compute and add the bulk velocity forcing at the end so that the 
      !           f(1:3) accounts for all the terms in the prediction
      !
      if(any(is_forced(:))) then
        call bulk_forcing_src(bulk_ftype,is_forced,n(1),n(2),n(3),dx,dy,dz,f_t12, &
                              nh_d,nh_u,dzc,dzf,rho,u,v,w,f)
        dpdl_c(1:3) = dpdl_c(1:3) + f(1:3) 
      endif
      !
      call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,no_outflow,dl,dzc,dzf,u,v,w)
      !
      call profiler_stop("RK")
      !
      ! 2c. construct the Poisson equation
      !
#if defined(_OPENACC)
      !$acc parallel loop collapse(3) present(p, pold)
#else
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,j,k) &
      !$OMP SHARED(p,pold)
#endif
      do k=1,n3
        do j=1,n2
          do i=1,n1
            pold(i,j,k) = p(i,j,k)
          enddo
        enddo
      enddo
#if defined(_OPENACC)
      !$acc end parallel loop 
#else
      !$OMP END PARALLEL DO
#endif
      call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,pold)
      !
      call fillps(n(1),n(2),n(3),nh_d,nh_u,dxi,dyi,dzi,dzfi,f_t12_i,rho_sp,u,v,w,p)
      !
      call updt_rhs_b(n(1),n(2),n(3),(/'c','c','c'/),cbcpre,nh_p,rhsbp_x,rhsbp_y,rhsbp_z,p)
      !
      call profiler_start("SOLVER", tag = .true., tag_color = COLOR_GREEN)
      !
#if defined(_OPENACC)
      call solver_gpu(n_z,dims_xyz(:,3),arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,:),(/'c','c','c'/),p)
#else
      call solver_cpu(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),p)
#endif
      call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
      !
      call profiler_stop("SOLVER")
      !
      call profiler_start("CORREC")
      !
      ! 2d. correct the velocity and update the pressure
      !
      call correc(n(1),n(2),n(3),nh_d,nh_u,dxi,dyi,dzi,dzci,f_t12,rho_sp,p,u,v,w,rho)
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
      !
      call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
      !
      call profiler_stop("CORREC")
      !
    enddo
    !
    write(fldnum,'(i9.9)') istep
    !
    ! 3. post-processing
    !
#if defined(_DO_POSTPROC)
    !
    call profiler_start("POSTPROC", tag = .true., tag_color = COLOR_WHITE)
    !
    if((mod(istep,wall_deltai).eq.0).and.do_wall) then
      call wall_avg(avg_dir,n(1),n(2),n(3),ng(1),ng(2),ng(3),dli(1),dli(2),dli(3),nh_t,ka,tmp,time)
    endif
    !
    if((mod(istep,time_deltai).eq.0).and.do_avg) then
      if(.not.allocated(vorx)) then ! allocate if this was not done before
        allocate(vorx( n(1),n(2),n(3)) , &
                 uv( n(1),n(2),n(3))   , &
                 vw( n(1),n(2),n(3))   , &
                 wu( n(1),n(2),n(3))   , &
                 utmp( n(1),n(2),n(3)) , &
                 vtmp( n(1),n(2),n(3)) , &
                 wtmp( n(1),n(2),n(3)) )
        if(    avg_dir.eq.1) then
         include 'allocation_x.h90'
        elseif(avg_dir.eq.2) then 
         include 'allocation_y.h90'
        else
         include 'allocation_z.h90'
        endif
      endif
      i_av = i_av + 1
      call compute_vorticity(n(1),n(2),n(3),dli(1),dli(2),dli(3),nh_u,v,w,vorx)
      call mixed_variables(n(1),n(2),n(3),dli(1),dli(2),dli(3),nh_u,nh_t,u,v,w,tmp, &
                           utmp,vtmp,wtmp,uv,vw,wu)
     !include 'time_averaging_gas.h90'
     !include 'time_averaging_heat_transfer.h90'
     if(mod(istep,iout1d).eq.0) i_av = 0
    endif
    !
    call profiler_stop("POSTPROC")
    !
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
        call chkdt_sp(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
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
      !
    endif
    !
    call profiler_stop("STEP")
    !
    ! output routines below
    !
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
      if(any(is_forced(:))) then 
        call cmpt_mean(n(1),n(2),n(3),nh_d,nh_u,dx,dy,dzf,lx,ly,lz,u,meanvelu)
        call cmpt_mean(n(1),n(2),n(3),nh_d,nh_u,dx,dy,dzf,lx,ly,lz,v,meanvelv)
        call cmpt_mean(n(1),n(2),n(3),nh_d,nh_u,dx,dy,dzf,lx,ly,lz,w,meanvelw)
        var(:)   = 0._rp
        var(1)   = time
        var(2)   = 1._rp*istep
        var(3:5) = -dpdl_c(1:3)*dti*rho_sp
        var(6)   = sqrt(maxval(abs(var(3:5)))*lz*0.5_rp) ! u_tau
        var(7)   = rho_sp*var(6)*0.5_rp*lz/mu_sp ! Re_tau
        var(8)   = meanvelu
        var(9)   = meanvelv
        var(10)  = meanvelw
        var(11)  = rho_sp*var(8)*lz/mu_sp ! Re_bulk
        call out0d(trim(datadir)//'forcing.out',11,var)
      endif 
      !
      call profiler_stop("OUT:iout0d")
      !
    endif
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
#if !defined(_BENCHMARK_NO_IO)
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
#if defined(_HEAT_TRANSFER)
      !@cuf istat = cudaMemPrefetchAsync(tmp, size(tmp), cudaCpuDeviceId, 0)   
      !@cuf istat = cudaMemPrefetchAsync(dtmpdtrkold, size(dtmpdtrkold), cudaCpuDeviceId, 0)  
#endif
      !
      inquire(file='data/fldu.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldu.bin   data/fldu_old.bin')
      call load(action_load,trim(restart_dir)//'fldu.bin',n,      u(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(u, size(u), mydev, 0)
      !
      inquire(file='data/fldv.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldv.bin   data/fldv_old.bin')
      call load(action_load,trim(restart_dir)//'fldv.bin',n,      v(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(v, size(v), mydev, 0)
      !
      inquire(file='data/fldw.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldw.bin   data/fldw_old.bin')
      call load(action_load,trim(restart_dir)//'fldw.bin',n,      w(1:n(1),1:n(2),1:n(3)))   
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(w, size(w), mydev, 0)
      !
      inquire(file='data/flddp.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldp.bin   data/fldp_old.bin')
      call load(action_load,trim(restart_dir)//'fldp.bin',n,      p(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(p, size(p), mydev, 0)
      !
      if(time_scheme.ne.'rk3') then ! RK3 is self-restarting
        inquire(file='data/flddu.bin', exist=is_data)
        if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddu.bin  data/flddu_old.bin')
        call load(action_load,trim(restart_dir)//'flddu.bin',n,dudtrko(1:n(1),1:n(2),1:n(3)))
        !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(dudtrko, size(dudtrko), mydev, 0)
        !
        inquire(file='data/flddv.bin', exist=is_data)
        if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddv.bin  data/flddv_old.bin')
        call load(action_load,trim(restart_dir)//'flddv.bin',n,dvdtrko(1:n(1),1:n(2),1:n(3)))
        !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(dvdtrko, size(dvdtrko), mydev, 0)
        !
        inquire(file='data/flddw.bin', exist=is_data)
        if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddw.bin  data/flddw_old.bin')
        call load(action_load,trim(restart_dir)//'flddw.bin',n,dwdtrko(1:n(1),1:n(2),1:n(3)))
        !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(dwdtrko, size(dwdtrko), mydev, 0)
      endif
      !
#if defined(_HEAT_TRANSFER)
      inquire(file='data/fldtmp.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldtmp.bin data/fldtmp_old.bin')
      call load(action_load,trim(restart_dir)//'fldtmp.bin',n,    tmp(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(tmp, size(tmp), mydev, 0)
      !
      if(time_scheme.ne.'rk3') then ! RK3 is self-restarting
        inquire(file='data/flddtmp.bin', exist=is_data)
        if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddtmp.bin data/flddtmp_old.bin')
        call load(action_load,trim(restart_dir)//'flddtmp.bin',n,dtmpdtrkold(1:n(1),1:n(2),1:n(3)))
        !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(dtmpdtrkold, size(dtmpdtrkold), mydev, 0)
      endif
#endif
      !
      inquire(file='data/scalar.out', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/scalar.out data/scalar_old.out')
      call load_scalar(action_load,trim(restart_dir)//'scalar.out',time,istep,dto)
      !
      if(myid.eq.0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
      !
      call profiler_stop("OUT:isave")
      !
    endif
#endif
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
  if(myid.eq.0.and.(.not.kill)) print*, '*** Fim ***'
  !
  call profiler_report()
  !
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
  call exit
  !
end program flutas
