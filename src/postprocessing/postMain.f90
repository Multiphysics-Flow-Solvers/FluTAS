program flutas
  use iso_c_binding  , only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound      , only: boundp,updt_rhs_b,bounduvw
  use mod_common_mpi , only: myid,ierr, coord, comm_cart
  use mod_initmpi    , only: initmpi
  use mod_load       , only: load
  use mod_initgrid   , only: initgrid
  use mod_output     , only: out0d,out1d,out1d_2,out2d,out3d!,out3d_c

  use mod_param      , only: itot,jtot,ktot,lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,uref,lref,rey,visc,small, &
  cbcvel,bcvel,cbcpre,bcpre, &
  icheck,iout0d,iout1d,iout2d,iout3d,isave, &
  nstep,time_max,tw_max,stop_type,restart, &
  rkcoeff, &
  datadir, &
  cfl,     &
  inivel,  &
  imax,jmax,dims, &
  nthreadsmax, &
  gr, &
  is_outflow,no_outflow,is_forced,bforce, &
  rho1,rho2,mu1,mu2,cbcvof,bcvof, lateInit, i_lateInit, &
  doSpectra, doFlux, doSF, sf_deltai, spectra_deltai, flux_deltai, nbin,doHist, hist_deltai, doBalance, &
  balance_deltai, doRR, RR_deltai, doTagging, tagging_deltai,&
  doTime,time_deltai,&
  doWall,wall_deltai,&
  n,ng,l,dl,dli, &
  read_input
  
  use mod_sanity     , only: test_sanity
  use mod_solver     , only: solver
  use mod_types
  use mod_vof        , only: initvof,advvof,update_vof,update_property             
  use mod_spectra, only: post_globalT
  use mod_post,    only: histograms,energyBalance_mf, &
                         time_avg,wall_avg,compute_vorticity,mixed_variables
  use mod_tagging, only: dropletTagging
  !$ use omp_lib
  implicit none
!
  real(rp), dimension(:,:,:), allocatable   :: pold
  real(rp), dimension(:,:,:), allocatable   :: psi,kappa,mu,rho      
  real(rp), dimension(:,:,:,:), allocatable :: cur                   
  real(rp), dimension(:,:,:,:), allocatable :: nor                   
  real(rp), dimension(:,:,:), allocatable   :: d_thinc
  real(rp) :: dto 
  character(len=1) :: action_load
!
  real(rp), allocatable, dimension(:,:,:) :: u,v,w,p
!
  integer  :: i,j,k,im,ip,jm,jp,km,kp
! 
  real(rp) :: ristep
  real(rp) :: dt,dti,dtmax,time,dtrk,dtrki,divtot,divmax
  integer  :: irk,istep
  real(rp), allocatable, dimension(:) :: dzc,dzf,zc,zf,dzci,dzfi,dzflzi,dzclzi,zclzi
! 
  real(rp) :: meanvel,meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: dpdl
  !real(rp), allocatable, dimension(:) :: var
  real(rp), dimension(10) :: var
! 
  real(rp) :: dt12,dt12av,dt12min,dt12max
! 
  real(rp) :: f1d,f2d,f3d
  real(rp) :: twi,tw
  character(len=7) :: fldnum
  integer :: kk
  logical :: is_done,kill,is_data
! 
  real(rp)   , allocatable, dimension(:,:,:) :: ux,vx,wx,vofx ,uy,vy,wy,vofy, uz,vz,wz,vofz 
  real(rp)   , allocatable, dimension(:,:,:) :: uv,vw,wu 
  real(rp)   , allocatable, dimension(:,:,:) :: vorx 
  real(rp)   , allocatable, dimension(:) :: u_avg,v_avg,w_avg,u_sqr,v_sqr,w_sqr 
  real(rp)   , allocatable, dimension(:) :: uv_avg,vw_avg,wu_avg,uv_sqr,vw_sqr,wu_sqr 
  real(rp)   , allocatable, dimension(:) :: vorx_avg,vorx_sqr 
  real(rp)   , allocatable, dimension(:) :: void_avg,void_sqr 
  real(rp)   :: u_vol_avg,u_vol_sqr,v_vol_avg,v_vol_sqr,w_vol_avg,w_vol_sqr 
  real(rp)   :: uv_vol_avg,uv_vol_sqr,vw_vol_avg,vw_vol_sqr,wu_vol_avg,wu_vol_sqr 
  real(rp)   :: vorx_vol_avg,vorx_vol_sqr 
  real(rp)   :: vory_vol_avg,vory_vol_sqr 
  real(rp)   :: vorz_vol_avg,vorz_vol_sqr 
  real(rp)   :: void_vol_avg,void_vol_sqr 
  integer    :: reasonU, reasonV, reasonW, reasonP, reasonVof
  integer, dimension(4) ::nfiles
  character(len=100) :: ufile, vfile, wfile, pfile, voffile,dummy
  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  !
  ! read parameter file
  !
  call read_input(myid)
  !
  !
  ! Create data folder, if it does not exist.
  !
  inquire(file='data/',exist=is_data)
  if (is_data .eqv. .false.)  call system('mkdir data')
  inquire(file='data/post/',exist=is_data)
  if ((is_data .eqv. .false.).and.(myid.eq.0))  call system('mkdir data/post')
  inquire(file='data/post/spectra/',exist=is_data)
  if ((doSpectra.eqv..true.).and.(is_data .eqv. .false.).and.(myid.eq.0))   call system('mkdir data/post/spectra')
  inquire(file='data/post/flux/',exist=is_data)
  if ((doFlux.eqv..true.).and.(is_data .eqv. .false.).and.(myid.eq.0))   call system('mkdir data/post/flux')
  inquire(file='data/post/hist/',exist=is_data)
  if ((doHist.eqv..true.).and.(is_data .eqv. .false.).and.(myid.eq.0))   call system('mkdir data/post/hist')
  inquire(file='data/post/balance/',exist=is_data)
  if ((doBalance.eqv..true.).and.(is_data .eqv. .false.).and.(myid.eq.0))   then 
    call system('mkdir data/post/balance')
    open(92,file='data/post/balance/balance.out',position='append')
    write(92,'(A)') '# istep, 1-E_E, 2-ens_l, 3-ens_g, 4-eps_l, 5-eps_g,'&
    //'6-KE_l, 7-KE_g, 8-KEp_l, 9-KEp_g, 10-Tnu_l, 11-Tnu_g, 12-Tp_l, 13-Tp_g,   14-Psi_nu'
    close(92)
  endif
  inquire(file='data/post/sf/',exist=is_data)
  if ((doSF.eqv..true.).and.(is_data .eqv. .false.).and.(myid.eq.0)) call system('mkdir data/post/sf')
  inquire(file='data/post/RR/',exist=is_data)
  if ((doRR.eqv..true.).and.(is_data .eqv. .false.).and.(myid.eq.0)) call system('mkdir data/post/RR')
  inquire(file='data/post/tagging/',exist=is_data)
  if ((doTagging.eqv..true.).and.(is_data .eqv. .false.).and.(myid.eq.0)) call system('mkdir data/post/tagging')
  inquire(file='data/post/time_averaging/',exist=is_data)
  if ((doTime.eqv..true.).and.(is_data .eqv. .false.).and.(myid.eq.0)) call system('mkdir data/post/time_averaging')
  inquire(file='data/post/wall/',exist=is_data)
  if ((doWall.eqv..true.).and.(is_data .eqv. .false.).and.(myid.eq.0)) call system('mkdir data/post/wall')
  ! initialize MPI/OpenMP
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(ng,cbcpre)
  n  = (/imax,jmax,ktot/) ! now set in initmpi
  twi = MPI_WTIME()
  !
  ! allocate memory
  !
  !
  allocate(p( 0:n(1)+1,0:n(2)+1,0:n(3)+1) , &
   u( 0:n(1)+1,0:n(2)+1,0:n(3)+1) , &
   v( 0:n(1)+1,0:n(2)+1,0:n(3)+1) , &
   w( 0:n(1)+1,0:n(2)+1,0:n(3)+1)) 
  allocate(psi(    0:imax+1,0:jmax+1,0:ktot+1)   , &
    kappa(  0:imax+1,0:jmax+1,0:ktot+1)   , &
    mu(     0:imax+1,0:jmax+1,0:ktot+1)   , &
    rho(    0:imax+1,0:jmax+1,0:ktot+1)   , &
    cur(    0:imax+1,0:jmax+1,0:ktot+1,7) , &
    nor(    0:imax+1,0:jmax+1,0:ktot+1,3) , &
    d_thinc(0:imax+1,0:jmax+1,0:ktot+1)   )
  allocate(dzc(   0:n(3)+1), &
   dzf(   0:n(3)+1), &
   zc(    0:n(3)+1), &
   zf(    0:n(3)+1), &
   dzci(  0:n(3)+1), &
   dzfi(  0:n(3)+1), &
   dzclzi(0:n(3)+1), &
   dzflzi(0:n(3)+1), &
   zclzi( 0:n(3)+1))

   allocate(ux( ng(1),ng(2)/dims(1),ng(3)/dims(2)) , &
   vx( ng(1),ng(2)/dims(1),ng(3)/dims(2)) , &
   wx( ng(1),ng(2)/dims(1),ng(3)/dims(2)) , &
 vofx( ng(1),ng(2)/dims(1),ng(3)/dims(2)) , &
   uy( ng(1)/dims(1),ng(2),ng(3)/dims(2)) , &
   vy( ng(1)/dims(1),ng(2),ng(3)/dims(2)) , &
   wy( ng(1)/dims(1),ng(2),ng(3)/dims(2)) , &
 vofy( ng(1)/dims(1),ng(2),ng(3)/dims(2)) , &
   uz( n(1),n(2),n(3)) , &
   vz( n(1),n(2),n(3)) , &
   wz( n(1),n(2),n(3)) , &
 vofz( n(1),n(2),n(3)) , &
 vorx( n(1),n(2),n(3)) , &
   uv( n(1),n(2),n(3)) , &
   vw( n(1),n(2),n(3)) , &
   wu( n(1),n(2),n(3)) , &
   u_avg (n(3)), &
   v_avg (n(3)), &
   w_avg (n(3)), &
   u_sqr (n(3)), &
   v_sqr (n(3)), &
   w_sqr (n(3)), &
   uv_avg (n(3)), &
   vw_avg (n(3)), &
   wu_avg (n(3)), &
   uv_sqr (n(3)), &
   vw_sqr (n(3)), &
   wu_sqr (n(3)), &
   void_avg (n(3)), &
   void_sqr (n(3)), &
   vorx_avg (n(3)), &
   vorx_sqr (n(3)))

  if(myid.eq.0) print*, '***********************************'
  if(myid.eq.0) print*, '*** Beginning of postprocessing ***'
  if(myid.eq.0) print*, '***********************************'
  if(myid.eq.0) print*, ''
  !

  call initgrid(inivel,n(3),gr,lz,dzc,dzf,zc,zf)

  do k=0,ktot+1
   dzci(k)=1./dzc(k)
   dzfi(k)=1./dzf(k)
   dzflzi(k)=dzf(k)/lz
   dzclzi(k)=dzc(k)/lz
   zclzi(k)=zc(k)/lz
 enddo

! Listing files for u,v,w,p,vof to process
if (myid.eq.0) call system('rm nfiles')
if (myid.eq.0) call system('ls data/vex*>ufiles; cat ufiles|wc -l>>nfiles')
if (myid.eq.0) call system('ls data/vey*>vfiles; cat vfiles|wc -l>>nfiles')
if (myid.eq.0) call system('ls data/vez*>wfiles; cat wfiles|wc -l>>nfiles')
if (myid.eq.0) call system('ls data/pre*>pfiles; cat pfiles|wc -l>>nfiles')
if (myid.eq.0) call system('ls data/vof*>voffiles;  cat voffiles|wc -l>>nfiles')
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

open(30,file='nfiles'  ,action="read")
do i=1,4
  read(30,*) nfiles(i)
enddo
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

if (sum(nfiles)/4.ne.nfiles(1)) then
  if (myid.eq.0) print*, "Error in number of files"
  call exit
endif
 
open(31,FILE='ufiles'  ,action="read")
open(32,FILE='vfiles'  ,action="read")
open(33,FILE='wfiles'  ,action="read")
open(34,FILE='pfiles'  ,action="read")
open(35,FILE='voffiles',action="read")
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

spectra_deltai = 1
flux_deltai = 1
hist_deltai = 1
balance_deltai = 1
sf_deltai = 1
RR_deltai = 1
do i=1,nfiles(1)  
  read(31,'(A)') ufile
  read(32,'(A)') vfile
  read(33,'(A)') wfile
  read(34,'(A)') pfile
  read(35,'(A)') voffile
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  action_load = 'r'
  call load(action_load,trim(ufile),n,      u(1:n(1),1:n(2),1:n(3)))
  call load(action_load,trim(vfile),n,      v(1:n(1),1:n(2),1:n(3)))
  call load(action_load,trim(wfile),n,      w(1:n(1),1:n(2),1:n(3)))
  call load(action_load,trim(pfile),n,      p(1:n(1),1:n(2),1:n(3)))
  call load(action_load,trim(voffile),n,  psi(1:n(1),1:n(2),1:n(3)))
  if(myid.eq.0) print*, '*** Checkpoint loaded at timestep = ', ufile(14:20)
  call bounduvw(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
  ! post-process and write initial condition
  !
  call update_vof(n,dli,dzc,dzf,psi,nor,cur,d_thinc)
  kappa(:,:,:) = cur(:,:,:,7)
  !#if defined(_OPENACC)
  call update_property(n,(/rho1,rho2/),psi,rho)
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,rho) 
  call update_property(n,(/mu1 ,mu2 /),psi,mu ) 
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu ) 
  read(ufile(14:20),*) istep
  if(myid.eq.0) print*, 'vof updated'
  !
  ! main loop
  !
  if((doSpectra.eqv..true.).or. &
     (doFlux.eqv..true.)   .or. &
     (doSF  .eqv..true.)   .or. &
     (doRR  .eqv..true.))   then 
    if(myid.eq.0) print*, 'Doing trasnpose base postprocessing'
    call post_globalT(ng, dims, n, u(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                                   v(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                                   w(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                                 psi(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                                   ux,vx,wx,uy,vy,wy,uz,vz,wz,vofx,vofy,vofz, istep)
  endif
  if((doHist.eqv..true.)) then
    if(myid.eq.0) print*, 'Doing histogrmas '
    call histograms(u(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                    v(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                    w(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                    n, nbin, 0,istep, .true.)
  endif
  if((doBalance.eqv..true.)) then
    if(myid.eq.0) print*, 'Doing Energy balance'
    call energyBalance_mf(u(0:n(1)+1,0:n(2)+1,0:n(3)+1), & 
                          v(0:n(1)+1,0:n(2)+1,0:n(3)+1), & 
                          w(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                          p,psi,kappa,rho,mu,n,time, dzc(0:n(3)+1), dzf(0:n(3)+1))
  endif
  if((doTagging.eqv..true.)) then
    if(myid.eq.0) print*, 'Doing droplet tagging'
    call dropletTagging(n,psi,u(0:n(1)+1,0:n(2)+1,0:n(2)+1), &
                              v(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                              w(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                        istep)
  endif
  if((doWall.eqv..true.)) then
    call wall_avg(n,ng,dli,mu,   u(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                                 w(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                  time)
  endif
  if((doTime.eqv..true.)) then
    include 'postprocessing/time_averaging.h90'
  endif
  if(myid.eq.0) print*, "Elapsed time ", MPI_WTIME()-twi
  twi = MPI_WTIME()
enddo
  if(myid.eq.0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
  call exit
end program flutas
