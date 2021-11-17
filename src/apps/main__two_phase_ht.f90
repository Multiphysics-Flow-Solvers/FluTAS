!
! SPDX-License-Identifier: MIT
!
!--------------------------------------------------------------------------------------------------------------------
! FluTAS -- Fluid Transport Accelerated Solver                            
!                                                                                       
!  a. Menu Title: FluTAS_VoF_ht;                                                         
!  b. Feautures of FluTAS_VoF_ht:                                                        
!      --> two-fluid solver with VoF with heat transfer;
!      --> both Oberbeck-Boussinesq approximations and low Mach in the gas phase; 
!      --> allow for density and viscosity mismatch;                                    
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
  use mod_param     , only: lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,small, &
                            cbcvel,bcvel,cbcpre,bcpre, &
                            icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                            nstep,time_max,tw_max,stop_type,restart, &
                            cfl,     &
                            constant_dt,dt_input, &
                            inivel,  &
                            itot,jtot,ktot,dims_in, &
                            nthreadsmax, &
                            gr, &
                            is_outflow,no_outflow,is_forced,bforce, &
                            rho1,rho2,rho0,mu1,mu2,cbcvof,bcvof,late_init,i_late_init, &
#if defined(_DO_POSTPROC)
                            do_avg,time_deltai,&
                            avg_dir,&
                            do_favre,&
                            do_wall,wall_deltai,&
#endif
                            initmp, &
                            cbctmp,bctmp, &
                            cp1,cv1,kappa1, &
                            cp2,cv2,kappa2, &
                            tl0,tg0, &
                            rho0, &
                            n,ng,l,dl,dli, &
                            read_input, dpdl
  !
  use mod_sanity    , only: test_sanity
  use mod_initflow  , only: inittmp
  use mod_rks       , only: ab2_tmp
#if defined(_OPENACC)
  use mod_solver    , only: solver
#else
  use mod_solver_cpu, only: solver_cpu
#endif
  use mod_types
  use mod_vof       , only: initvof,advvof,update_vof,update_property,alloc_vof_var
#if defined(_DO_POSTPROC)
  use mod_post      , only: time_avg,wall_avg,compute_vorticity,mixed_variables
#endif
#if defined(_OPENACC)
  use cudafor
  use mod_common_mpi, only: mydev
#endif
#if defined(_USE_NVTX)
  use nvtx
#endif
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
  real(rp), dimension(:,:,:)  , allocatable :: tmp
  real(rp), dimension(:,:,:)  , allocatable :: dtmpdtrk
  real(rp), dimension(:,:,:)  , allocatable :: dtmpdtrkold
  real(rp), dimension(:,:,:)  , allocatable :: ka,cpp
  real(rp), dimension(:), allocatable :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi
  real(rp), dimension(:), allocatable :: dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g
  !
  real(rp), dimension(:,:,:), pointer :: dudtrko,dvdtrko,dwdtrko
  !
  type(C_PTR), dimension(2,2) :: arrplanp
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:)   :: ap,bp,cp
  real(rp) :: normfftp
  ! 
  type rhs_bound
#if defined(_OPENACC)
   real(rp), allocatable, dimension(:,:,:), managed :: x,y,z
#else
   real(rp), allocatable, dimension(:,:,:) :: x,y,z
#endif
  end type rhs_bound 
  type(rhs_bound), allocatable :: rhsbp
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
  integer  :: nh_d,nh_u,nh_up,nh_p,nh_v,nh_t

  integer  :: i,j,k,im,ip,jm,jp,km,kp
  integer  :: irk,istep
  character(len=1) :: action_load
  logical  :: is_first_vel,is_data, is_first_tmp
  !
  real(rp) :: meanvel,meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: dpdl_c
  real(rp), dimension(10) :: var
  ! 
  character(len=9)   :: fldnum
  character(len=100) :: datadir
  real(rp) :: f1d,f2d,f3d
  real(rp) :: twi,tw
  integer :: kk
  logical :: is_done,kill
  real(rp) :: rho0i
#if defined(_OPENACC)
  attributes(managed) :: pold, kappa, mu, rho, psi, d_thinc, cur_t, nor
  attributes(managed) :: u, v, w, p, up, vp, wp, pp 
  attributes(managed) :: tmp, ka, cpp, dtmpdtrk, dtmpdtrkold
  integer :: istat
  integer(kind=cuda_count_kind) :: freeMem, totMem
  integer(rp) :: totEle
  attributes(managed) :: dzc, dzf, zc, zf, dzci, dzfi, lambdaxyp, ap, bp, cp, rhsbp
  attributes(managed) :: dudtrko, dvdtrko, dwdtrko
#endif
  !if we don't use dropcheck.f90 we can comment the next line  
  !real(rp) :: xd,yd,zd,ut,vt,wt,zcd,ycd,xcd,vol
  real(rp) :: vof_min,vof_max,vol_p1
  real(rp) :: aux
  integer :: ii,jj!,kk
#if defined(_DO_POSTPROC)
  real(rp), allocatable, dimension(:,:,:) :: uv,vw,wu 
  real(rp), allocatable, dimension(:,:,:) :: vorx 
  real(rp), allocatable, dimension(:) :: u_avg_g,v_avg_g,w_avg_g,u_sqr_g,v_sqr_g,w_sqr_g 
  real(rp), allocatable, dimension(:) :: u_avg_l,v_avg_l,w_avg_l,u_sqr_l,v_sqr_l,w_sqr_l 
  real(rp), allocatable, dimension(:) :: uv_avg_g,vw_avg_g,wu_avg_g,uv_sqr_g,vw_sqr_g,wu_sqr_g 
  real(rp), allocatable, dimension(:) :: uv_avg_l,vw_avg_l,wu_avg_l,uv_sqr_l,vw_sqr_l,wu_sqr_l 
  real(rp), allocatable, dimension(:) :: vorx_avg_g,vorx_sqr_g 
  real(rp), allocatable, dimension(:) :: vorx_avg_l,vorx_sqr_l 
  real(rp), allocatable, dimension(:) :: void_avg,void_sqr 
  real(rp) :: u_vol_avg_g,u_vol_sqr_g,v_vol_avg_g,v_vol_sqr_g,w_vol_avg_g,w_vol_sqr_g 
  real(rp) :: u_vol_avg_l,u_vol_sqr_l,v_vol_avg_l,v_vol_sqr_l,w_vol_avg_l,w_vol_sqr_l 
  real(rp) :: uv_vol_avg_g,uv_vol_sqr_g,vw_vol_avg_g,vw_vol_sqr_g,wu_vol_avg_g,wu_vol_sqr_g 
  real(rp) :: uv_vol_avg_l,uv_vol_sqr_l,vw_vol_avg_l,vw_vol_sqr_l,wu_vol_avg_l,wu_vol_sqr_l 
  real(rp) :: vorx_vol_avg_g,vorx_vol_sqr_g 
  real(rp) :: vorx_vol_avg_l,vorx_vol_sqr_l 
  real(rp) :: vory_vol_avg_g,vory_vol_sqr_g 
  real(rp) :: vory_vol_avg_l,vory_vol_sqr_l 
  real(rp) :: vorz_vol_avg_g,vorz_vol_sqr_g 
  real(rp) :: vorz_vol_avg_l,vorz_vol_sqr_l 
  real(rp) :: void_vol_avg,void_vol_sqr 
  integer  :: i_av 
  real(rp), allocatable, dimension(:,:,:) :: utmp,vtmp,wtmp 
  real(rp), allocatable, dimension(:) :: tmp_avg_g,tmp_sqr_g 
  real(rp), allocatable, dimension(:) :: tmp_avg_l,tmp_sqr_l 
  real(rp), allocatable, dimension(:) :: utmp_avg_g,vtmp_avg_g,wtmp_avg_g,utmp_sqr_g,vtmp_sqr_g,wtmp_sqr_g  
  real(rp), allocatable, dimension(:) :: utmp_avg_l,vtmp_avg_l,wtmp_avg_l,utmp_sqr_l,vtmp_sqr_l,wtmp_sqr_l  
  real(rp) :: tmp_vol_avg_g,tmp_vol_sqr_g 
  real(rp) :: tmp_vol_avg_l,tmp_vol_sqr_l 
  real(rp) :: utmp_vol_avg_g,utmp_vol_sqr_g,vtmp_vol_avg_g,vtmp_vol_sqr_g,wtmp_vol_avg_g,wtmp_vol_sqr_g 
  real(rp) :: utmp_vol_avg_l,utmp_vol_sqr_l,vtmp_vol_avg_l,vtmp_vol_sqr_l,wtmp_vol_avg_l,wtmp_vol_sqr_l 
#if defined(_OPENACC)
  attributes(managed) :: vorx,uv,vw,wu,utmp,vtmp,wtmp
  attributes(managed) :: u_avg_g,v_avg_g,w_avg_g,u_sqr_g,v_sqr_g,w_sqr_g 
  attributes(managed) :: u_avg_l,v_avg_l,w_avg_l,u_sqr_l,v_sqr_l,w_sqr_l 
  attributes(managed) :: uv_avg_g,vw_avg_g,wu_avg_g,uv_sqr_g,vw_sqr_g,wu_sqr_g 
  attributes(managed) :: uv_avg_l,vw_avg_l,wu_avg_l,uv_sqr_l,vw_sqr_l,wu_sqr_l 
  attributes(managed) :: vorx_avg_g,vorx_sqr_g 
  attributes(managed) :: vorx_avg_l,vorx_sqr_l 
  attributes(managed) :: void_avg,void_sqr 
  attributes(managed) :: tmp_avg_g,tmp_sqr_g 
  attributes(managed) :: tmp_avg_l,tmp_sqr_l 
  attributes(managed) :: utmp_avg_g,vtmp_avg_g,wtmp_avg_g,utmp_sqr_g,vtmp_sqr_g,wtmp_sqr_g  
  attributes(managed) :: utmp_avg_l,vtmp_avg_l,wtmp_avg_l,utmp_sqr_l,vtmp_sqr_l,wtmp_sqr_l  
#endif
#endif
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  !
  ! read parameter file
  !
#if defined(_USE_NVTX)
  call nvtxStartRange("readInput")
#endif
  call read_input(myid)
#if defined(_USE_NVTX)
  call nvtxEndRange
#endif
  !
  ! create data folder and subfolders for post-processing, if they do not exist.
  !
  inquire(file='data/',exist=is_data)
  if(.not.is_data.and.myid.eq.0) call execute_command_line('mkdir -p data')
  datadir = 'data/'
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
#if defined(_USE_NVTX)
  call nvtxStartRange("initmpi")
#endif
  call initmpi(.false.,ng,cbcpre,dims_in,dims_xyz,dims,n)
#if defined(_USE_NVTX)
  call nvtxEndRange
#endif
  !
  twi = MPI_WTIME()
  !
  ! halo calculation (for the variables we generally 
  !                   do enlarge the stencil support)
  !
  nh_u = 3
  nh_t = 3
  nh_d = nh_u ! take the maximum of the previous (e.g., max(nh_u,nh_t)) 
  !
  ! allocate memory
  !
#if defined(_USE_NVTX)
  call nvtxStartRange("variable_allocations")
#endif
  allocate(p(0:n(1)+1,0:n(2)+1,0:n(3)+1) , &
           u(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           v(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           w(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           up(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           vp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           wp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           pp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(tmp(1-nh_t:n(1)+nh_t,1-nh_t:n(2)+nh_t,1-nh_t:n(3)+nh_t), &
            ka(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           cpp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
              dtmpdtrk(1:n(1),1:n(2),1:n(3)),&
           dtmpdtrkold(1:n(1),1:n(2),1:n(3)))
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
  allocate(rhsbp)
  allocate(rhsbp%x(n(2),n(3),0:1), &
           rhsbp%y(n(1),n(3),0:1), &
           rhsbp%z(n(1),n(2),0:1))
  !
  ! prefetching of the variables (TODO: remember to add the one of x-pencil!)
  !
#if defined(_OPENACC)
  istat = cudaMemAdvise(u, size(u), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(v, size(v), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(w, size(w), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(p, size(p), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(up, size(up), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(vp, size(vp), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(wp, size(wp), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(pp, size(pp), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(dudtrko, size(dudtrko), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(dvdtrko, size(dvdtrko), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(dwdtrko, size(dwdtrko), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemPrefetchAsync(u, size(u), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(v, size(v), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(w, size(w), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(p, size(p), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(up, size(up), mydev, 0)
  istat = cudaMemPrefetchAsync(vp, size(vp), mydev, 0)
  istat = cudaMemPrefetchAsync(wp, size(wp), mydev, 0)
  istat = cudaMemPrefetchAsync(pp, size(pp), mydev, 0)
  istat = cudaMemPrefetchAsync(dudtrko, size(dudtrko), mydev, 0)
  istat = cudaMemPrefetchAsync(dvdtrko, size(dvdtrko), mydev, 0)
  istat = cudaMemPrefetchAsync(dwdtrko, size(dwdtrko), mydev, 0)
  !
  istat = cudaMemAdvise(tmp, size(tmp), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise( ka, size( ka), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(cpp, size(cpp), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(   dtmpdtrk,    size(dtmpdtrk), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(dtmpdtrkold, size(dtmpdtrkold), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemPrefetchAsync(tmp, size(tmp), mydev, 0)
  istat = cudaMemPrefetchAsync(ka, size(ka), mydev, 0)
  istat = cudaMemPrefetchAsync(cpp, size(cpp), mydev, 0)
  istat = cudaMemPrefetchAsync(   dtmpdtrk,    size(dtmpdtrk), mydev, 0)
  istat = cudaMemPrefetchAsync(dtmpdtrkold, size(dtmpdtrkold), mydev, 0)
  !
  istat = cudaMemAdvise(rhsbp%x, size(rhsbp%x), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(rhsbp%y, size(rhsbp%y), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(rhsbp%z, size(rhsbp%z), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemPrefetchAsync(rhsbp%x, size(rhsbp%x), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(rhsbp%y, size(rhsbp%y), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(rhsbp%z, size(rhsbp%z), cudaCpuDeviceId, 0)
  !
  istat = cudaMemAdvise(zc, size(zc), cudaMemAdviseSetReadMostly, 0)
  istat = cudaMemAdvise(zf, size(zf), cudaMemAdviseSetReadMostly, 0)
  istat = cudaMemAdvise(dzc, size(dzc), cudaMemAdviseSetReadMostly, 0)
  istat = cudaMemAdvise(dzf, size(dzf), cudaMemAdviseSetReadMostly, 0)
  istat = cudaMemAdvise(dzci, size(dzci), cudaMemAdviseSetReadMostly, 0)
  istat = cudaMemAdvise(dzfi, size(dzfi), cudaMemAdviseSetReadMostly, 0)
  !
#if defined(_USE_VOF)
  istat = cudaMemAdvise(mu, size(mu), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(pold, size(pold), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(rho, size(rho), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(kappa, size(kappa), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(psi, size(psi), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(d_thinc, size(d_thinc), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(nor, size(nor), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(cur_t, size(cur_t), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemPrefetchAsync(pold, size(pold), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(mu, size(mu), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(rho, size(rho), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(kappa, size(kappa), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(psi, size(psi), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(d_thinc, size(d_thinc), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(nor, size(nor), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(cur_t, size(cur_t), cudaCpuDeviceId, 0)
#endif
  !
#if defined(_DO_POSTPROC)
  istat = cudaMemAdvise(vorx, size(vorx), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(  uv,   size(uv), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(  vw,   size(vw), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(  wu,   size(wu), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(utmp, size(utmp), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(vtmp, size(vtmp), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemAdvise(wtmp, size(wtmp), cudaMemAdviseSetPreferredLocation, mydev)
  istat = cudaMemPrefetchAsync(vorx, size(vorx), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(  uv,   size(uv), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(  vw,   size(vw), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(  wu,   size(wu), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(utmp, size(utmp), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(vtmp, size(vtmp), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(wtmp, size(wtmp), cudaCpuDeviceId, 0)
#endif
  !
#endif
#if defined(_USE_NVTX)
  call nvtxEndRange
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
  ! halo generation
  !
  call halo_gen(n,nh_u ,halo_u )
  call halo_gen(n,nh_t ,halo_t )
  call halo_gen(n,nh_p ,halo_p )
  call halo_gen(n,nh_up,halo_up)
  call halo_gen(n,nh_v ,halo_v )
  call halo_gen(n,nh_d ,halo_d )
  !
#if defined(_OPENACC)
  !
  ! Allocate buffers for GPU halo communications
  !
  call alloc_buf(n,nh_d)
  !
#endif
  !
  ! initialize the grid (using global variables along z)
  !
#if defined(_OPENACC)
  istat = cudaMemPrefetchAsync(zc, size(zc), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(zf, size(zf), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(dzc, size(dzc), cudaCpuDeviceId, 0)
  istat = cudaMemPrefetchAsync(dzf, size(dzf), cudaCpuDeviceId, 0)
#endif
#if defined(_USE_NVTX)
  call nvtxStartRange("initgrid")
#endif
  call initgrid(inivel,ng(3),gr,lz,nh_d,dzc_g,dzf_g,zc_g,zf_g) 
#if defined(_USE_NVTX)
  call nvtxEndRange
#endif
  !
#if defined(_OPENACC)
  istat = cudaMemPrefetchAsync(zc, size(zc), mydev, 0)
  istat = cudaMemPrefetchAsync(zf, size(zf), mydev, 0)
  istat = cudaMemPrefetchAsync(dzc, size(dzc), mydev, 0)
  istat = cudaMemPrefetchAsync(dzf, size(dzf), mydev, 0)
#endif
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
  !$acc kernels
  do k=1-nh_d,ng(3)+nh_d
    dzfi_g(k) = 1._rp/dzf_g(k)
    dzci_g(k) = 1._rp/dzc_g(k)
  enddo
  !$acc end kernels
  !@cuf istat=cudaDeviceSynchronize()
  !
#if defined(_OPENACC)
  istat = cudaMemPrefetchAsync(dzci, size(dzci), mydev, 0)
  istat = cudaMemPrefetchAsync(dzfi, size(dzfi), mydev, 0)
  istat = cudaMemPrefetchAsync(zc, size(zc), mydev, 0)
#endif
  !
  ! compute the spacing along z in local coordinates
  !
  !$acc kernels
  do k=1-nh_d,n(3)+nh_d
    kk      = k + ijk_start(3)
    zc(k)   = zc_g(kk)
    zf(k)   = zf_g(kk)
    dzf(k)  = dzf_g(kk)
    dzc(k)  = dzc_g(kk)
    dzfi(k) = 1._rp/dzf(k)
    dzci(k) = 1._rp/dzc(k)
  enddo
  !$acc end kernels
  !@cuf istat=cudaDeviceSynchronize()
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
    is_first_tmp = .true.
    is_first_vel = .true.
    !
#if defined(_OPENACC)
    istat = cudaMemPrefetchAsync(dzci, size(dzci), cudaCpuDeviceId, 0)
    istat = cudaMemPrefetchAsync(dzfi, size(dzfi), cudaCpuDeviceId, 0)
#endif
    !
    if(.not.late_init) then 
      !
      ! Initialize VoF 
      !
#if defined(_USE_NVTX)
      call nvtxStartRange('initvof')
#endif
      call initvof(n,dli,psi)
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,psi)
      call alloc_vof_var(n(1),n(2),n(3))
#if defined(_USE_NVTX)
      call nvtxEndRange
#endif
      !
    endif
    !
    ! Initialize temperature 
    !
#if defined(_USE_NVTX)
    call nvtxStartRange('inittmp')
#endif
    call inittmp(initmp,n(1),n(2),n(3),nh_v,nh_t,dims,.false.,0.5d0,psi,tmp)
    call boundp(cbctmp,n,bctmp,nh_d,nh_t,halo_t,dl,dzc,dzf,tmp)
#if defined(_USE_NVTX)
    call nvtxEndRange
#endif
    !
#if defined(_USE_NVTX)
    call nvtxEndRange
#endif
    call initflow(inivel,n(1),n(2),n(3),dims,nh_d,nh_u,nh_p,zc/lz,dzc/lz,dzf/lz,u,v,w,p)
    !
    ! set to zeros the rhs of momentum equation 
    ! (only for the first time-step, not for the restarting)
    ! 
    !$acc kernels
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dudtrko(i,j,k) = 0._rp
          dvdtrko(i,j,k) = 0._rp
          dwdtrko(i,j,k) = 0._rp
        enddo
      enddo
    enddo
    !$acc end kernels
    !@cuf istat=cudaDeviceSynchronize()
    !
    ! set to zeros the rhs of temperature equation 
    ! (only for the first time-step, not for the restarting)
    ! 
    !$acc kernels
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dtmpdtrkold(i,j,k) = 0._rp
        enddo
      enddo
    enddo
    !$acc end kernels
    !@cuf istat=cudaDeviceSynchronize()
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
    call load(action_load,trim(datadir)//'fldtmp.bin',n,    tmp(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'flddtmp.bin',n,dtmpdtrkold(1:n(1),1:n(2),1:n(3)))
    call load_scalar(action_load,trim(datadir)//'scalar.out',time,istep,dto)
    !
    if(myid.eq.0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
    !
    is_first_vel = .false.
    is_first_tmp = .false.
    !
    call alloc_vof_var(n(1),n(2),n(3)) ! to put here (important for the restarting)
    !
  endif
#if defined(_OPENACC)
  istat = cudaMemPrefetchAsync(u, size(u), mydev, 0)
  istat = cudaMemPrefetchAsync(v, size(v), mydev, 0)
  istat = cudaMemPrefetchAsync(w, size(w), mydev, 0)
  istat = cudaMemPrefetchAsync(p, size(p), mydev, 0)
  istat = cudaMemPrefetchAsync(pold, size(pold), mydev, 0)
#endif
  !
  ! set boundary conditions on the initial/loaded fields
  !
  call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
  call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,psi)
  call boundp(cbctmp,n,bctmp,nh_d,nh_t,halo_t,dl,dzc,dzf,tmp)
  !
  ! for the first time-step and the restarting we use a 0th order extrapolation 
  ! in time-splitting of the pressure equation
  !
  !$acc kernels 
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        pold(i,j,k) = p(i,j,k)
      end do
    end do
  end do
  !$acc end kernels
  !@cuf istat=cudaDeviceSynchronize()
  call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,pold)
  !
  ! update the quantities derived from VoF using the lastest available VoF field
  !
  call update_vof(n,dli,nh_d,dzc,dzf,nh_v,halo_v,psi,nor,cur_t,kappa,d_thinc)
  call update_property(n,(/rho1,rho2/),psi,rho)
  call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)
  call update_property(n,(/mu1 ,mu2 /),psi,mu ) 
  call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,mu)
  call update_property(n,(/kappa1,kappa2/),psi,ka)           
  call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,ka )
  call update_property(n,(/cp1,cp2/),psi,cpp)           
  call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,cpp)
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i9.9)') istep
  include 'out1d.h90'
  include 'out2d.h90'
  include 'out3d.h90'
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
  ! preliminary checks
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
  ! initialize Poisson solver
  ! and deallocate global arrays (not needed anymore) 
  !
  call initsolver(n,dims,dims_xyz(:,3),dli,nh_d,dzci_g,dzfi_g,cbcpre,bcpre(:,:),(/'c','c','c'/),lambdaxyp, & 
                  ap,bp,cp,arrplanp,normfftp,rhsbp%x,rhsbp%y,rhsbp%z)
  deallocate(dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g)
  !
  ! main loop
  !
  if(myid.eq.0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  !
  do while(.not.is_done)
    !
#if defined(_TIMING)
    dt12  = MPI_WTIME()
#endif
    !
    istep = istep + 1
    time  = time + dt
    !
    if(myid.eq.0) print*, 'Timestep #', istep, 'Time = ', time
    !
    if(late_init.and.(istep.eq.i_late_init)) then
      !
      call initvof(n,dli,psi)
      call update_vof(n,dli,nh_d,dzc,dzf,nh_v,halo_v,psi,nor,cur_t,kappa,d_thinc)
      call alloc_vof_var(n(1),n(2),n(3))
      !
      call update_property(n,(/mu1 ,mu2 /),psi,mu ) 
      call update_property(n,(/rho1,rho2/),psi,rho)
      call update_property(n,(/kappa1,kappa2/),psi,ka)           
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,mu )
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,ka )
      call chkdt(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
      !
    endif
    !
    ! 0. save quantities from the previous time-step
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(p,pold,pp)
#endif
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          pp(i,j,k)   = (1._rp+(dt/dto))*p(i,j,k) - (dt/dto)*pold(i,j,k)
          pold(i,j,k) = p(i,j,k)
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end kernels
    !@cuf istat=cudaDeviceSynchronize()
#else
    !$OMP END PARALLEL DO
#endif
    call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,pp  )
    call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,pold)
    !
    ! 1. VoF advection and properties update --> vof^(n+1) 
    !
    if((.not.late_init).or.        &              ! Normal run (initialized at i=0)
       (late_init.and.istep.ge.i_late_init)) then ! in case we want first to run without advecting the VoF

      !
#if defined(_USE_NVTX)
      call nvtxStartRange("FullVof")
#endif
      call advvof(n,dli,dt,halo_v,nh_d,nh_u,dzc,dzf,u,v,w,psi,nor,cur_t,kappa,d_thinc)
      !
      call update_property(n,(/mu1 ,mu2 /),psi,mu )           
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,mu )
      call update_property(n,(/rho1,rho2/),psi,rho)           
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)
      call update_property(n,(/kappa1,kappa2/),psi,ka)           
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,ka )
      call update_property(n,(/cp1,cp2/),psi,cpp)           
      call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,cpp)
#if defined(_USE_NVTX)
      call nvtxEndRange
#endif
      !
    endif
    !
    ! 2. Temperature advection and thermodynamic pressure --> tmp^(n+1)
    !
#if defined(_USE_NVTX)
    call nvtxStartRange("ab2_tmp")
#endif
    call ab2_tmp(is_first_tmp,n(1),n(2),n(3),dli(1),dli(2),dli(3),nh_d,nh_t,dt,dto,rho,cpp,ka,u,v,w, &
                 tmp,dtmpdtrk,dtmpdtrkold)
#if defined(_USE_NVTX)
    call nvtxEndRange
#endif
    call boundp(cbctmp,n,bctmp,nh_d,nh_t,halo_t,dl,dzc,dzf,tmp)
    !
#if !defined(_VOF_DBG)
    !
    ! 3. two-fluid Navier-Stokes --> (u,v,w,p)^(n+1)
    !
    rho0i = 1._rp/rho0
    dpdl_c(1:3) = 0._rp
    !
#if defined(_USE_NVTX)
    call nvtxStartRange("rk")
#endif
    call rk(is_first_vel,n(1),n(2),n(3),dims,dxi,dyi,dzi,nh_d,nh_u,dt,dto,rho0i,dzci,dzfi,u,v,w,p,pp, &
            kappa,psi,mu,rho, & 
#if defined(_BOUSSINESQ)
            nh_t,tmp,& 
#endif
            dudtrko,dvdtrko,dwdtrko,up,vp,wp,f)
#if defined(_USE_NVTX)
    call nvtxEndRange
#endif
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_up,halo_up,no_outflow,dl,dzc,dzf,up,vp,wp)
    dpdl_c(1:3) = dpdl_c(1:3) + f(1:3) + dpdl(1:3)*dt
    !
#if defined(_USE_NVTX)
    call nvtxStartRange("fillps")
#endif
    call fillps(n(1),n(2),n(3),nh_d,nh_up,dxi,dyi,dzi,dzci,dzfi,dti,rho0,up,vp,wp,p)
    !
#if defined(_USE_NVTX)
    call nvtxEndRange
#endif
    !
    call updt_rhs_b(n(1),n(2),n(3),(/'c','c','c'/),cbcpre,nh_p,rhsbp%x,rhsbp%y,rhsbp%z,p)
#if defined(_USE_NVTX)
    call nvtxStartRange("solver")
#endif
    !
#if defined(_OPENACC)
    call solver(n_z,dims_xyz(:,3),arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,:),(/'c','c','c'/),p)
#else
    call solver_cpu(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),p)
#endif
#if defined(_USE_NVTX)
    call nvtxEndRange
#endif
    call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
#if defined(_USE_NVTX)
    call nvtxStartRange("CorrecEcc")
#endif
    !
    call correc(n(1),n(2),n(3),nh_d,nh_u,dxi,dyi,dzi,dzci,dt,rho0,p,up,vp,wp,rho,u,v,w)
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf,u,v,w)
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP WORKSHARE
#endif
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          p(i,j,k) = pold(i,j,k) + p(i,j,k)
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end kernels
    !@cuf istat=cudaDeviceSynchronize()
#else
    !$OMP END WORKSHARE
#endif
    call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
    !
#if defined(_USE_NVTX)
    call nvtxEndRange
#endif
    !
#endif 
    ! 
    write(fldnum,'(i9.9)') istep
    !
    !include 'dropcheck.h90'
    !
    ! post-processing
    !
#if defined(_DO_POSTPROC)
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
     include 'time_averaging_liquid.h90'
     include 'time_averaging_gas.h90'
     include 'time_averaging_heat_transfer.h90'
    endif
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
      !
    endif
    !
    ! output routines below
    !
    if(mod(istep,iout0d).eq.0) then
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
    endif
    !
    if(mod(istep,iout1d).eq.0) then
     include 'out1d.h90'
    endif
    if(mod(istep,iout2d).eq.0) then
     include 'out2d.h90'
    endif
    if(mod(istep,iout3d).eq.0) then
     include 'out3d.h90'
    endif
    !
    if(mod(istep,isave).eq.0.or.(is_done.and..not.kill)) then
      !
      action_load = 'w'
      !
      inquire(file='data/fldu.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldu.bin   data/fldu_old.bin')
      inquire(file='data/fldv.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldv.bin   data/fldv_old.bin')
      inquire(file='data/fldw.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldw.bin   data/fldw_old.bin')
      inquire(file='data/flddu.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddu.bin  data/flddu_old.bin')
      inquire(file='data/flddv.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddv.bin  data/flddv_old.bin')
      inquire(file='data/flddw.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddw.bin  data/flddw_old.bin')
      inquire(file='data/flddp.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldp.bin   data/fldp_old.bin')
      inquire(file='data/fldpsi.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldpsi.bin data/fldpsi_old.bin')
      inquire(file='data/fldtmp.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldtmp.bin data/fldtmp_old.bin')
      inquire(file='data/flddtmp.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddtmp.bin data/flddtmp_old.bin')
      inquire(file='data/scalar.out', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/scalar.out data/scalar_old.out')
      !
      call load(action_load,trim(datadir)//'fldu.bin',n,      u(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldv.bin',n,      v(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldw.bin',n,      w(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'flddu.bin',n,dudtrko(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'flddv.bin',n,dvdtrko(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'flddw.bin',n,dwdtrko(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldp.bin',n,      p(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldpsi.bin',n,    psi(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldtmp.bin',n,    tmp(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'flddtmp.bin',n,dtmpdtrkold(1:n(1),1:n(2),1:n(3)))
      call load_scalar(action_load,trim(datadir)//'scalar.out',time,istep,dto)
      !
      if(myid.eq.0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
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
             tmp,ka,cpp,dtmpdtrk,dtmpdtrkold, &     
             cur_t, &                   
             nor, &                 
             d_thinc, &
             dzc,dzf,zc,zf,dzci,dzfi, &
             dudtrko,dvdtrko,dwdtrko)
#if defined(_DO_POSTPROC)
  deallocate(vorx,uv,vw,wu,utmp,vtmp,wtmp, &
             u_avg_g,v_avg_g,w_avg_g,      &
             u_avg_l,v_avg_l,w_avg_l,      &
             u_sqr_g,v_sqr_g,w_sqr_g,      &
             u_sqr_l,v_sqr_l,w_sqr_l,      &
             tmp_avg_g,tmp_sqr_g,          &
             tmp_avg_l,tmp_sqr_l,          &
             utmp_avg_g,utmp_sqr_g,        &
             utmp_avg_l,utmp_sqr_l,        &
             vtmp_avg_g,vtmp_sqr_g,        &
             vtmp_avg_l,vtmp_sqr_l,        &
             wtmp_avg_g,wtmp_sqr_g,        &
             wtmp_avg_l,wtmp_sqr_l,        &
             uv_avg_g,vw_avg_g,wu_avg_g,   &                
             uv_avg_l,vw_avg_l,wu_avg_l,   &             
             uv_sqr_g,vw_sqr_g,wu_sqr_g,   &             
             uv_sqr_l,vw_sqr_l,wu_sqr_l,   &             
             void_avg,void_sqr,            &
             vorx_avg_g,vorx_sqr_g,        &              
             vorx_avg_l,vorx_sqr_l)             
#endif
  !
  deallocate(lambdaxyp)
  deallocate(ap,bp,cp)
  deallocate(rhsbp%x,rhsbp%y,rhsbp%z)
  !
  if(myid.eq.0.and.(.not.kill)) print*, '*** Fim ***'
  !
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
  call exit
  !
end program flutas
