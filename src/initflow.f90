!
! SPDX-License-Identifier: MIT
!
module mod_initflow
  !
  use mpi
  use decomp_2d
  use mod_common_mpi, only: ierr,myid,ijk_start
  use mod_sanity    , only: flutas_error
  use mod_types
#if defined(_OPENACC)
  use cudafor
  use openacc_curand
#endif
  !
  implicit none
  !
  private
  public  :: initflow, add_noise
#if defined(_HEAT_TRANSFER)
  public  :: inittmp
#endif
  !
  contains
  !
  subroutine initflow(inivel,nx,ny,nz,dims,is_noise_vel,noise_vel,nh_d,nh_u,nh_p,rho_p,mu_p,zclzi,dzclzi,dzflzi,u,v,w,p)
    !
    ! computes initial conditions for the velocity and pressure field
    !
    use mod_param, only: pi,dx,dy,dz, &
                         is_forced,bvel_x,bvel_y,bvel_z, &
                         dpdl_x,dpdl_y,dpdl_z,bulk_ftype,bcvel, &
#if defined(_TURB_FORCING)
                         u0_t,abc_x,abc_y,abc_z,k0_t,add_noise_abc, &
#endif
                         lx,ly,lz,is_wallturb,wallturb_type
    !
    implicit none
    !
    character(len=3), intent(in )                                     :: inivel
    integer         , intent(in )                                     :: nx,ny,nz
    integer         , intent(in ), dimension(3)                       :: dims
    logical         , intent(in )                                     :: is_noise_vel
    real(rp)        , intent(in )                                     :: noise_vel
    integer         , intent(in )                                     :: nh_d,nh_u,nh_p
    real(rp)        , intent(in )                                     :: rho_p,mu_p ! when we call initflow, we decide which phase
    real(rp)        , intent(in ), dimension(1-nh_d:)                 :: zclzi,dzclzi,dzflzi
    real(rp)        , intent(out), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp)        , intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    !
    !@cuf attributes(managed) :: u, v, w, p, zclzi, dzclzi, dzflzi
    !
#if defined(_TURB_FORCING)
    real(rp) :: add_sin,add_cos
#endif
    real(rp), allocatable, dimension(:) :: u1d
    real(rp) :: q
    real(rp) :: xc,yc,zc,xf,yf,zf,norm_u
    integer  :: i,j,k,ijk_start_x,ijk_start_y,ijk_start_z
    logical  :: is_mean
    !
    allocate(u1d(nz))
    !
    ijk_start_x = ijk_start(1)
    ijk_start_y = ijk_start(2)
    ijk_start_z = ijk_start(3)
    !
    q       = 0.5_rp
    is_mean = .false.
    norm_u  = 1._rp ! default case: initial velocity magnitude equal to 1
    !
    ! we compute the appropriate norm_u depending on 'inivel'
    !
    if(any(is_forced(:))) then
      !
      ! right norm_u for most of wall-bounded flows except for Couette
      !
      select case(bulk_ftype)
      case('cfr')
        if(    is_forced(1)) then
          norm_u = bvel_x
        elseif(is_forced(2)) then
          norm_u = bvel_y
        elseif(is_forced(3)) then 
          norm_u = bvel_z
        endif
      case('cpg')
        norm_u = ((-(dpdl_x+dpdl_y+dpdl_z)*lz/2._rp/rho_p)**0.5*lz/2._rp/mu_p/0.09) 
        norm_u = norm_u**(1._rp/0.88)*mu_p/lz ! compute norm_u to be consistent with the imposed dpdl_*
      end select
      !
    else
      !
      ! other relevant cases:
      !   1. Triperiodic cases: norm_u is not employed since we use u0_t. So we skip.
      !   2. Gravity-driven flow (e.g. Rayleigh-Benard convection). We typically start with 
      !      initial velocity equal to 0, i.e. case('zer') and trigger the scalar field to promote transition. So we skip. 
      !
    endif
    !
    ! we compute the initial velocity field
    !
    select case(inivel)
    case('cou')
      norm_u = abs( bcvel(0,3,1)-bcvel(1,3,1) )   ! assumes that the z-direction is the wall-normal one 
      call couette(   q,nz,nh_d,zclzi,norm_u,u1d)
    case('poi')
      norm_u = norm_u
      call poiseuille(q,nz,nh_d,zclzi,norm_u,u1d)
      is_mean = .true.
    case('zer')
      norm_u = 0._rp
      u1d(:) = norm_u
    case('uni')
      norm_u = 1._rp
      u1d(:) = norm_u
    case('log')
      norm_u = norm_u
      call log_profile(q,nz,nh_d,zclzi,mu_p,rho_p,lz,norm_u,u1d)
      is_mean = .true.
    case('hcp')
      norm_u = norm_u
      deallocate(u1d)
      allocate(u1d(2*nz))
      call poiseuille(q,2*nz,nh_d,zclzi,norm_u,u1d)
      is_mean = .true.
#if defined(_TURB_FORCING)
    case('tgv')
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            !
            zc = (k+ijk_start_z-0.5_rp)*dz/lz*2.0_rp*pi
            yc = (j+ijk_start_y-0.5_rp)*dy/ly*2.0_rp*pi
            yf = (j+ijk_start_y-0.0_rp)*dy/ly*2.0_rp*pi
            xc = (i+ijk_start_x-0.5_rp)*dx/lx*2.0_rp*pi
            xf = (i+ijk_start_x-0.0_rp)*dx/lx*2.0_rp*pi
            !
            u(i,j,k) =  u0_t*sin(k0_t*xf)*cos(k0_t*yc)*cos(k0_t*zc)
            v(i,j,k) = -u0_t*cos(k0_t*xc)*sin(k0_t*yf)*cos(k0_t*zc)
            w(i,j,k) = 0.0_rp
            p(i,j,k) = 0.0_rp!(cos(2._rp*xc)+cos(2._rp*yc))*(cos(2._rp*zc)+2._rp)/16._rp
            !
          enddo
        enddo
      enddo
      !$acc end kernels
    case('abc')
      if(add_noise_abc) then
        add_sin = 0.5_rp ! to promote transition when Re_abc is close to the critical value
        add_cos = 1.5_rp ! to promote transition when Re_abc is close to the critical value
      else
        add_sin = 0.0_rp
        add_cos = 0.0_rp
      endif
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            !
            zc = (k+ijk_start_z-0.5_rp)*dz/lz*2.0_rp*pi
            yc = (j+ijk_start_y-0.5_rp)*dy/ly*2.0_rp*pi
            yf = (j+ijk_start_y-0.0_rp)*dy/ly*2.0_rp*pi
            xc = (i+ijk_start_x-0.5_rp)*dx/lx*2.0_rp*pi
            xf = (i+ijk_start_x-0.0_rp)*dx/lx*2.0_rp*pi
            !
            u(i,j,k) =  u0_t*(abc_x*sin(k0_t*zc+add_sin) + abc_z*cos(k0_t*yc+add_cos))
            v(i,j,k) =  u0_t*(abc_y*sin(k0_t*xc+add_sin) + abc_x*cos(k0_t*zc+add_cos))
            w(i,j,k) =  u0_t*(abc_z*sin(k0_t*yc+add_sin) + abc_y*cos(k0_t*xc+add_cos))
            p(i,j,k) = 0._rp!(cos(2._rp*xc)+cos(2._rp*yc))*(cos(2._rp*zc)+2._rp)/16._rp
            !
          enddo
        enddo
      enddo
      !$acc end kernels
#endif
    case default
      call flutas_error('Error: invalid name of the initial velocity field. Simulation aborted. Check dns.in')
    end select
    !
    ! we set the initial velocity field
    !   note: in 2D, we solve the governing equations along the y and z direction. The x-velocity is set to 0 
    !
    if((inivel.ne.'tgv').and.(inivel.ne.'abc')) then
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
#if defined(_TWOD)
            u(i,j,k) = 0._rp
            v(i,j,k) = u1d(k)
#else
            u(i,j,k) = u1d(k)
            v(i,j,k) = 0._rp
#endif
            w(i,j,k) = 0._rp
            p(i,j,k) = 0._rp
          enddo
        enddo
      enddo
      !$acc end kernels
    endif
    !
    ! set the mean velocity
    !
    if(is_mean) then
#if defined(_TWOD)
      call set_mean(nx,ny,nz,nh_d,lx,ly,lz,dx,dy,dz,dzclzi,norm_u,v(1:nx,1:ny,1:nz))
#else
      call set_mean(nx,ny,nz,nh_d,lx,ly,lz,dx,dy,dz,dzclzi,norm_u,u(1:nx,1:ny,1:nz))
#endif
    endif
    !
    ! we superimpose a disturbance to the velocity components (v,w) 
    ! This option is very convenient in triggering transition for wall-bounded turbulent flows
    !
    if(is_wallturb) then
      !
      select case(wallturb_type)
      case('hkv')
        !
        ! initialize a streamwise vortex pair for a fast transition
        ! to turbulence in a pressure-driven channel:
        !        psi(x,y,z)  = f(z)*g(x,y), with
        !        f(z)        = (1-z**2)**2, and
        !        g(x,y)      = y*exp[-(16x**2-4y**2)]
        ! (x,y,z) --> (streamwise, spanwise, wall-normal) directions
        !
        ! see Henningson and Kim, JFM 1991
        !
        !!$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              !
              zc = 2._rp*zclzi(k) - 1._rp ! z rescaled to be between -1 and +1
              zf = 2._rp*(zclzi(k) + 0.5_rp*dzflzi(k)) - 1._rp
              yc = ((ijk_start_y+j-0.5_rp)*dy-0.5_rp*ly)*2._rp/lz
              yf = ((ijk_start_y+j-0.0_rp)*dy-0.5_rp*ly)*2._rp/lz
              xc = ((ijk_start_x+i-0.5_rp)*dx-0.5_rp*lx)*2._rp/lz
              xf = ((ijk_start_x+i-0.0_rp)*dx-0.5_rp*lx)*2._rp/lz
              !
              !u(i,j,k) = u1d(k) ! we do not touch the x-velocity
              v(i,j,k) =  1._rp * fz(zc)*dgxy(yf,xc)*norm_u*1.5_rp
              w(i,j,k) = -1._rp * gxy(yc,xc)*dfz(zf)*norm_u*1.5_rp
              p(i,j,k) = 0._rp
              !
            enddo
          enddo
        enddo
        !!$acc end kernels
      case('tgv')
        !
        ! trigger wall-bounded turbulence using Taylor-Green vortices
        !
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              !
              zc = (zclzi(k)                 )*2._rp*pi
              zf = (zclzi(k)+0.5_rp*dzclzi(k))*2._rp*pi
              yc = (j+ijk_start_y-0.5_rp)*dy/ly*2._rp*pi
              yf = (j+ijk_start_y-0.0_rp)*dy/ly*2._rp*pi
              xc = (i+ijk_start_x-0.5_rp)*dx/lx*2._rp*pi
              xf = (i+ijk_start_x-0.0_rp)*dx/lx*2._rp*pi
              !
              !u(i,j,k) = u1d(k) ! we do not touch the x-velocity
              v(i,j,k) =  sin(xc)*cos(yf)*cos(zc)*norm_u
              w(i,j,k) = -cos(xc)*sin(yc)*cos(zf)*norm_u
              p(i,j,k) = 0._rp!(cos(2.*xc)+cos(2.*yc))*(cos(2.*zc)+2.)/16.
              !
            enddo
          enddo
        enddo
        !$acc end kernels
      case default
        call flutas_error('Wrong setting specified for is_wallturb. Option available: hkv or tgv')
      end select
      !
    endif
    deallocate(u1d)
    !
    ! add noise
    !
    if(is_noise_vel) then
      call add_noise(nx,ny,nz,dims,123,noise_vel,u(1:nx,1:ny,1:nz))
      call add_noise(nx,ny,nz,dims,456,noise_vel,v(1:nx,1:ny,1:nz))
      call add_noise(nx,ny,nz,dims,789,noise_vel,w(1:nx,1:ny,1:nz))
      !
#if defined(_TWOD)
      !
      ! in the next, we ensure that if noise is applied, the initial condition is kept
      ! two-dimensional 
      ! 
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u(i,j,k) = 0._rp
            v(i,j,k) = v(nx,j,k)
            w(i,j,k) = w(nx,j,k)
          enddo
        enddo
      enddo
#endif
    endif
    ! 
    return
  end subroutine initflow
  !
#if defined(_HEAT_TRANSFER)
  subroutine inittmp(initmp,nx,ny,nz,nh_t,dims,is_noise_tmp,noise_tmp, &
#if defined(_USE_VOF)
                     nh_v,psi, &
#endif
                     tmp)
    !
    ! computes initial conditions for the temperature field
    !
    use mod_param, only: tmp0
#if defined(_USE_VOF)
    use mod_param, only: tl0,tg0
#endif
    !
    implicit none
    !
    character(len=3), intent(in )                                     :: initmp
    integer         , intent(in )                                     :: nx,ny,nz
    integer         , intent(in )                                     :: nh_t
    integer         , intent(in ), dimension(3)                       :: dims
    logical         , intent(in )                                     :: is_noise_tmp
    real(rp)        , intent(in )                                     :: noise_tmp
#if defined(_USE_VOF)
    integer         , intent(in )                                     :: nh_v
    real(rp)        , intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: psi
    !@cuf attributes(managed) :: psi
#endif
    real(rp)        , intent(out), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: tmp
    !
    !@cuf attributes(managed) :: tmp
    integer :: i,j,k
    !
    select case(initmp)
    case('uni')
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            tmp(i,j,k) = tmp0
          enddo
        enddo
      enddo
      !$acc end kernels
      !
#if defined(_USE_VOF)
    case('sin')
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            tmp(i,j,k) = tl0*psi(i,j,k)+tg0*(1._rp-psi(i,j,k))
          enddo
        enddo
      enddo
      !$acc end kernels
#endif
      !
    case default
      call flutas_error('Error: invalid name of the initial temperature field. Simulation aborted. Check heat_transfer.in')
    end select
    !
    if(is_noise_tmp) then
#if defined(_OPENACC)
      call add_noise_cuda(nx,ny,nz,dims,123456,noise_tmp,tmp(1:nx,1:ny,1:nz))
#else
      call add_noise(     nx,ny,nz,dims,123456,noise_tmp,tmp(1:nx,1:ny,1:nz))
#endif
      !
      ! in the next, we ensure that if noise is applied, the initial condition is kept
      ! two-dimensional 
      ! 
#if defined(_TWOD)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            tmp(i,j,k) = tmp(nx,j,k)
          enddo
        enddo
      enddo
#endif
    endif
    !
    return
  end subroutine inittmp
#endif   
  !
  subroutine add_noise(nx,ny,nz,dims,iseed,norm,p)
    !
    implicit none
    !
    integer , intent(in   )                      :: nx,ny,nz
    integer , intent(in   )                      :: iseed
    integer , intent(in   ), dimension(3)        :: dims
    real(rp), intent(in   )                      :: norm 
    real(rp), intent(inout), dimension(nx,ny,nz) :: p
    !
    integer(4), allocatable, dimension(:) :: seed
    real(rp) :: rn
    integer  :: nxg,nyg,nzg
    integer  :: i,j,k,ii,jj,kk
    integer  :: ijk_start_x,ijk_start_y,ijk_start_z
    !
    allocate(seed(64))
    seed(:) = iseed
    call random_seed( put = seed )
    !
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz*dims(3)
    !
    ijk_start_x = ijk_start(1)
    ijk_start_y = ijk_start(2)
    ijk_start_z = ijk_start(3)
    !
    do k=1,nzg
      do j=1,nyg
        do i=1,nxg
          !
          ii = i-ijk_start_x
          jj = j-ijk_start_y
          kk = k-ijk_start_z
          !
          call random_number(rn)
          if(ii.ge.1.and.ii.le.nx .and. &
             jj.ge.1.and.jj.le.ny .and. &
             kk.ge.1.and.kk.le.nz ) then
             p(ii,jj,kk) = p(ii,jj,kk) + 2.0_rp*(rn-0.5_rp)*norm
          endif
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine add_noise
  !
#if defined(_OPENACC)
  subroutine add_noise_cuda(nx,ny,nz,dims,iseed,norm,p)
    !
    implicit none
    !
    integer , intent(in   )                      :: nx,ny,nz
    integer , intent(in   )                      :: iseed
    integer , intent(in   ), dimension(3)        :: dims
    real(rp), intent(in   )                      :: norm 
    real(rp), intent(inout), dimension(nx,ny,nz) :: p
    !
    type(curandStateXORWOW) :: h
    real(rp) :: rn
    integer  :: nxg,nyg,nzg
    integer  :: i,j,k,ii,jj,kk
    integer  :: ijk_start_x,ijk_start_y,ijk_start_z
    !
    attributes(managed):: p
    !
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz*dims(3)
    !
    ijk_start_x = ijk_start(1)
    ijk_start_y = ijk_start(2)
    ijk_start_z = ijk_start(3)
    !
    !$acc parallel num_gangs(1) vector_length(1) private(h)
    call curand_init(iseed,0,0,h)
    !$acc loop seq
    do k=1,nzg
      do j=1,nyg
        do i=1,nxg
          !
          ii = i-ijk_start_x
          jj = j-ijk_start_y
          kk = k-ijk_start_z
          !
          rn = curand_uniform(h)
          if(ii.ge.1.and.ii.le.nx .and. &
             jj.ge.1.and.jj.le.ny .and. &
             kk.ge.1.and.kk.le.nz ) then
             p(ii,jj,kk) = p(ii,jj,kk) + 2._rp*(rn-0.5_rp)*norm
          endif
          !
        enddo
      enddo
    enddo
    !$acc end parallel
    !
    return
  end subroutine add_noise_cuda
#endif
  !
  subroutine set_mean(nx,ny,nz,nh_d,lx,ly,lz,dx,dy,dz,dzlzi,mean,p)
    !
    implicit none
    !
    integer , intent(in   )                      :: nx,ny,nz
    integer , intent(in   )                      :: nh_d
    real(rp), intent(in   )                      :: lx,ly,lz
    real(rp), intent(in   )                      :: dx,dy,dz
    real(rp), intent(in   ), dimension(1-nh_d:)  :: dzlzi 
    real(rp), intent(in   )                      :: mean
    real(rp), intent(inout), dimension(nx,ny,nz) :: p
    !
    real(rp) :: meanold
    integer  :: i,j,k
    !
    meanold = 0._rp
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,dzlzi) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP REDUCTION(+:meanold)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          meanold = meanold + p(i,j,k)*dx*dy*dzlzi(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    call mpi_allreduce(MPI_IN_PLACE,meanold,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    meanold = meanold/(1._rp*lx*ly) ! already divided by lz in the i,j,k loop
    !
    if(meanold.ne.0._rp) then
      !$OMP WORKSHARE
      p(1:nx,1:ny,1:nz) = p(1:nx,1:ny,1:nz)/meanold*mean
      !$OMP END WORKSHARE
    endif
    !
    return
  end subroutine set_mean
  !
  subroutine couette(q,n,nh_d,zc,norm,p)
    !
    ! plane couette profile normalized by the wall velocity difference
    !
    implicit none
    !
    real(rp), intent(in )                     :: q
    integer , intent(in )                     :: n
    integer , intent(in )                     :: nh_d
    real(rp), intent(in ), dimension(1-nh_d:) :: zc
    real(rp), intent(in )                     :: norm
    real(rp), intent(out), dimension(n)       :: p
    !
    real(rp) :: z
    integer  :: k
    !
    do k=1,n
      z    = zc(k)!1._rp*((k-1)+q)/(1._rp*n)
      p(k) = 0.5_rp*(1._rp-2._rp*z)/norm
    enddo
    !
    return
  end subroutine couette
  !
  subroutine poiseuille(q,n,nh_d,zc,norm,p)
    !
    implicit none
    !
    real(rp), intent(in )                     :: q
    integer , intent(in )                     :: n
    integer , intent(in )                     :: nh_d
    real(rp), intent(in ), dimension(1-nh_d:) :: zc
    real(rp), intent(in )                     :: norm
    real(rp), intent(out), dimension(n)       :: p
    !
    real(rp) :: z
    integer  :: k
    !
    ! plane poiseuille profile normalized by the bulk velocity
    !
    do k=1,n
      z    = zc(k)!1._rp*((k-1)+q)/(1._rp*n)
      p(k) = 6._rp*z*(1._rp-z)*norm
    enddo
    !
    return
  end subroutine poiseuille
  !
  subroutine log_profile(q,n,nh_d,zc,mu_p,rho_p,lref,uref,p)
    !
    implicit none
    !
    real(rp), intent(in )                     :: q
    integer , intent(in )                     :: n,nh_d
    real(rp), intent(in ), dimension(1-nh_d:) :: zc
    real(rp), intent(in )                     :: mu_p,rho_p,lref,uref
    real(rp), intent(out), dimension(n)       :: p
    !
    real(rp) :: z,reb,retau ! z/lz and bulk Reynolds number
    integer  :: k
    !
    reb = rho_p*uref*lref/mu_p
    retau = 0.09_rp*reb**(0.88) ! from Pope's book
    do k=1,n/2
      z    = zc(k)*2._rp*retau!1._rp*((k-1)+q)/(1._rp*n)*2.*retau
      p(k) = 2.5_rp*log(z) + 5.5_rp
      if(z.le.11.6_rp) p(k)=z
      p(n+1-k) = p(k)
    enddo
    !
    return
  end subroutine log_profile
  !
  ! functions to initialize the streamwise vortex pair
  ! (explained above)
  !
  function fz(zc)
  !!$acc routine(fz) seq
  real(rp), intent(in) :: zc
  real(rp) :: fz
    fz = ((1._rp-zc**2)**2)
  end function
  !
  function dfz(zc)
  !!$acc routine(dfz) seq
  real(rp), intent(in) :: zc
  real(rp) :: dfz
    dfz = -4._rp*zc*((1._rp-zc**2)**2)
  end function
  !
  function gxy(xc,yc)
  !!$acc routine(gxy) seq
  real(rp), intent(in) :: xc,yc
  real(rp) :: gxy
    gxy = yc*exp(-4._rp*(4._rp*xc**2+yc**2))
  end function
  !
  function dgxy(xc,yc)
  !!$acc routine(dgxy) seq
  real(rp), intent(in) :: xc,yc
  real(rp) :: dgxy
    dgxy = exp(-4._rp*(4._rp*xc**2+yc**2))*(1._rp-8._rp*yc**2)
  end function
  !
end module mod_initflow
