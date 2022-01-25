!
! SPDX-License-Identifier: MIT
!
module mod_initflow
  !
  use mpi
  use decomp_2d
  use mod_common_mpi, only: ierr,myid,ijk_start
  use mod_param     , only: pi,dx,dy,dz,ng, &
                            is_forced, bulk_velx,bulk_vely,bulk_velz,&
                            dpdl, small, &
#if defined(_TURB_FORCING)
                            u0_t, k0_t, abc, add_noise_abc, &
#endif
#if defined(_HEAT_TRANSFER)
                            tl0,tg0, &
#endif
                            lx,ly,lz,is_wallturb,mu2,rho2
  use mod_types
#if defined(_OPENACC)
  use cudafor
  use openacc_curand
#endif
  !
  implicit none
  !
  private
  public initflow, &
#if defined(_HEAT_TRANSFER)
         inittmp, &
#endif
         init_surf_tension,add_noise
  !
  contains
  !
  subroutine initflow(inivel,nx,ny,nz,dims,nh_d,nh_u,nh_p,zclzi,dzclzi,dzflzi,u,v,w,p)
    !
    ! computes initial conditions for the velocity and pressure field
    !
    implicit none
    !
    character(len=3), intent(in )                                     :: inivel
    integer         , intent(in )                                     :: nx,ny,nz
    integer         , intent(in ), dimension(3)                       :: dims
    integer         , intent(in )                                     :: nh_d,nh_u,nh_p
    real(rp)        , intent(in ), dimension(1-nh_d:)                 :: zclzi,dzclzi,dzflzi
    real(rp)        , intent(out), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp)        , intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    !@cuf attributes(managed) :: u, v, w, p, zclzi, dzclzi, dzflzi
    !
#if defined(_TURB_FORCING)
    real(rp) :: add_sin,add_cos
#endif
    real(rp), allocatable, dimension(:) :: u1d
    real(rp) :: q
    real(rp) :: xc,yc,zc,xf,yf,zf
    integer  :: i,j,k
    logical  :: is_noise,is_mean
    real(rp) :: norm = 1._rp
    !
    allocate(u1d(nz))
    ! TODO: prefetch/move on GPU async
    !
    is_noise = .false.
    is_mean  = .false.
    q = 0.5_rp
    if (is_forced(1))               then
      norm = bulk_velx
    elseif (is_forced(2))           then
      norm = bulk_vely
    elseif (is_forced(3))           then 
      norm = bulk_velz
    elseif (abs(sum(dpdl)).ge.small) then
      norm = ((-sum(dpdl)*lz/2._rp/rho2)**0.5*lz/2._rp/mu2/0.09)**(1._rp/0.88)* & !Re_tau
              mu2/lz
    else 
      norm = 1._rp
    endif
    !
    select case(inivel)
    case('cou')
      call couette(   q,nz,nh_d,zclzi,norm,u1d)
    case('poi')
      call poiseuille(q,nz,nh_d,zclzi,norm,u1d)
      is_mean=.true.
    case('zer')
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u(i,j,k) = 0.0_rp
            v(i,j,k) = 0.0_rp
            w(i,j,k) = 0.0_rp
          enddo
        enddo
      enddo
      !$acc end kernels
      u1d(:) = 0._rp
    case('uni')
      u1d(:) = 1._rp
    case('log')
      call log_profile(q,nz,nh_d,zclzi,mu2,rho2,lz,norm,u1d)
      is_noise = .true.
      is_mean  = .true.
    !case('hcl')
    !  deallocate(u1d)
    !  allocate(u1d(2*nz))
    !  call log_profile(q,2*nz,zclzi,visc,u1d)
    !  is_noise = .true.
    !  is_mean=.true.
    case('hcp')
      deallocate(u1d)
      allocate(u1d(2*nz))
      call poiseuille(q,2*nz,nh_d,zclzi,norm,u1d)
      is_mean = .true.
#if defined(_TURB_FORCING)
    case('tgv')
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            !
            zc = (k+ijk_start(3)-0.5_rp)*dz/lz*2.0_rp*pi
            yc = (j+ijk_start(2)-0.5_rp)*dy/ly*2.0_rp*pi
            yf = (j+ijk_start(2)-0.0_rp)*dy/ly*2.0_rp*pi
            xc = (i+ijk_start(1)-0.5_rp)*dx/lx*2.0_rp*pi
            xf = (i+ijk_start(1)-0.0_rp)*dx/lx*2.0_rp*pi
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
            zc = (k+ijk_start(3)-0.5_rp)*dz/lz*2.0_rp*pi
            yc = (j+ijk_start(2)-0.5_rp)*dy/ly*2.0_rp*pi
            yf = (j+ijk_start(2)-0.0_rp)*dy/ly*2.0_rp*pi
            xc = (i+ijk_start(1)-0.5_rp)*dx/lx*2.0_rp*pi
            xf = (i+ijk_start(1)-0.0_rp)*dx/lx*2.0_rp*pi
            !
            u(i,j,k) =  u0_t*(abc(1)*sin(k0_t*zc+add_sin) + abc(3)*cos(k0_t*yc+add_cos))
            v(i,j,k) =  u0_t*(abc(2)*sin(k0_t*xc+add_sin) + abc(1)*cos(k0_t*zc+add_cos))
            w(i,j,k) =  u0_t*(abc(3)*sin(k0_t*yc+add_sin) + abc(2)*cos(k0_t*xc+add_cos))
            p(i,j,k) = 0._rp!(cos(2._rp*xc)+cos(2._rp*yc))*(cos(2._rp*zc)+2._rp)/16._rp
            !
          enddo
        enddo
      enddo
      !$acc end kernels
#endif
    case default
      if(myid.eq.0) print*, 'ERROR: invalid name for initial velocity field'
      if(myid.eq.0) print*, ''
      if(myid.eq.0) print*, '*** Simulation abortited due to errors in the case file ***'
      if(myid.eq.0) print*, '    check the file dns.in'
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      call exit
    end select
    !
    if((inivel.ne.'tgv').and.(inivel.ne.'abc')) then
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u(i,j,k) = u1d(k)
            v(i,j,k) = 0._rp
            w(i,j,k) = 0._rp
            p(i,j,k) = 0._rp
          enddo
        enddo
      enddo
      !$acc end kernels
    endif
    !
    if(is_mean) then
      call set_mean(nx,ny,nz,nh_d,dims,dzclzi,norm,u(1:nx,1:ny,1:nz))
    endif
    if(is_wallturb.ne.0) then
      if(is_wallturb.eq.1) then
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
        do k=1,nz
          do j=1,ny
            do i=1,nx
              !
              zc = 2._rp*zclzi(k) - 1._rp ! z rescaled to be between -1 and +1
              zf = 2._rp*(zclzi(k) + 0.5_rp*dzflzi(k)) - 1._rp
              yc = ((ijk_start(2)+j-0.5_rp)*dy-0.5_rp*ly)*2._rp/lz
              yf = ((ijk_start(2)+j-0.0_rp)*dy-0.5_rp*ly)*2._rp/lz
              xc = ((ijk_start(1)+i-0.5_rp)*dx-0.5_rp*lx)*2._rp/lz
              xf = ((ijk_start(1)+i-0.0_rp)*dx-0.5_rp*lx)*2._rp/lz
              !
              u(i,j,k) = u1d(k)
              v(i,j,k) =  1._rp * fz(zc)*dgxy(yf,xc)*norm*1.5_rp
              w(i,j,k) = -1._rp * gxy(yc,xc)*dfz(zf)*norm*1.5_rp
              p(i,j,k) = 0._rp
              !
            enddo
          enddo
        enddo
      elseif(is_wallturb.eq.2) then
        !
        ! Initialize turbulence using Taylor-Green vortices
        !
        do k=1,nz
          zc = (zclzi(k)                 )*2._rp*pi
          zf = (zclzi(k)+0.5_rp*dzclzi(k))*2._rp*pi
          do j=1,ny
            yc = (j+ijk_start(2)-0.5_rp)*dy/ly*2._rp*pi
            yf = (j+ijk_start(2)-0.0_rp)*dy/ly*2._rp*pi
            do i=1,nx
              xc = (i+ijk_start(1)-0.5_rp)*dx/lx*2._rp*pi
              xf = (i+ijk_start(1)-0.0_rp)*dx/lx*2._rp*pi
              !u(i,j,k) = u1d(k)
              v(i,j,k) =  sin(xc)*cos(yf)*cos(zc)*norm
              w(i,j,k) = -cos(xc)*sin(yc)*cos(zf)*norm
              p(i,j,k) = 0._rp!(cos(2.*xc)+cos(2.*yc))*(cos(2.*zc)+2.)/16.
            enddo
          enddo
        enddo
      else
       if(myid.eq.0) print*, 'Wrong setting specified for is_wallturb. &
                              Please set 1 for Henningson and Kim vortices, &
                              2 for Taylor-Green vortices or 0 to disable it.'
       if(myid.eq.0) print*, 'Aborting...'
       call MPI_FINALIZE(ierr)
       call exit()
      endif
      !!$acc end kernels
    endif
    deallocate(u1d)
    !
    if(is_noise) then
      call add_noise(nx,ny,nz,dims,123,0.5_rp,u(1:nx,1:ny,1:nz))
      call add_noise(nx,ny,nz,dims,456,0.5_rp,v(1:nx,1:ny,1:nz))
      call add_noise(nx,ny,nz,dims,789,0.5_rp,w(1:nx,1:ny,1:nz))
      !
#if defined(_TWOD)
      !
      ! in the next, we ensure that if noise is applied, the initial condition is kept
      ! two-dimensional 
      ! 
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u(i,j,k) = u(nx/2,j,k)
            v(i,j,k) = v(nx/2,j,k)
            w(i,j,k) = w(nx/2,j,k)
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
  subroutine inittmp(initmp,nx,ny,nz,nh_v,nh_t,dims,is_noise,norm_t,vof,tmp)
    !
    ! computes initial conditions for the temperature field
    !
    implicit none
    !
    character(len=3), intent(in )                                     :: initmp
    integer         , intent(in )                                     :: nx,ny,nz
    integer         , intent(in )                                     :: nh_v,nh_t
    integer         , intent(in ), dimension(3)                       :: dims
    logical         , intent(in )                                     :: is_noise
    real(rp)        , intent(in )                                     :: norm_t
    real(rp)        , intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: vof
    real(rp)        , intent(out), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: tmp
#if defined(_OPENACC)
    attributes(managed):: tmp,vof
#endif
    !
    integer :: i,j,k
    !
    select case(initmp)
    case('uni')
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            tmp(i,j,k) = tg0
          enddo
        enddo
      enddo
      !$acc end kernels
      !
    case('sin')
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            tmp(i,j,k) = tg0*vof(i,j,k)+tl0*(1._rp-vof(i,j,k))
          enddo
        enddo
      enddo
      !$acc end kernels
      !
    case default
      if(myid.eq.0) print*, 'ERROR: invalid name for initial temperature field'
      if(myid.eq.0) print*, ''
      if(myid.eq.0) print*, '*** Simulation abortited due to ierrs in the case file ***'
      if(myid.eq.0) print*, '    check heat_transfer.in'
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      call exit
    end select
    !
    if(is_noise)then
#if defined(_OPENACC)
      call add_noise_cuda(nx,ny,nz,dims,123456,norm_t,tmp(1:nx,1:ny,1:nz))
#else
      call add_noise(     nx,ny,nz,dims,123456,norm_t,tmp(1:nx,1:ny,1:nz))
#endif
      !
      ! in the next, we ensure that if noise is applied, the initial condition is kept
      ! two-dimensional 
      ! 
#if defined(_TWOD)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            tmp(i,j,k) = tmp(nx/2,j,k)
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
  subroutine init_surf_tension(nx,ny,nz,dxi,dyi,dzi,nh_d,dzci,kappa,psi,rho, &
                               ssx_o,ssy_o,ssz_o,ssx,ssy,ssz)
    !
    ! computes initial conditions for the surface tension
    !
    use mod_param,  only: sigma
    !
    implicit none
    !
    integer , intent(in )                      :: nx,ny,nz
    real(rp), intent(in )                      :: dxi,dyi,dzi
    integer , intent(in )                      :: nh_d
    real(rp), intent(in ), dimension(1-nh_d:)  :: dzci
    real(rp), intent(in ), dimension(0:,0:,0:) :: kappa,psi,rho
    real(rp), intent(out), dimension(0:,0:,0:) :: ssx_o,ssy_o,ssz_o,ssx,ssy,ssz
    !
    real(rp) :: rhox,rhoy,rhoz,kappasx,kappasy,kappasz
    integer  :: i,j,k,ip,jp,kp
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,ip,jp,kp) &
    !$OMP PRIVATE(kappasx,kappasy,kappasz,rhox,rhoy,rhoz,sigma) &
    !$OMP SHARED(n,rho,psi,kappa,ssx_p,ssy_p,ssz_p,ssx,ssy,ssz)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i + 1
          jp = j + 1
          kp = k + 1
          !
          rhox    = 0.5_rp*(rho(ip,j,k)+rho(i,j,k))
          rhoy    = 0.5_rp*(rho(i,jp,k)+rho(i,j,k))
          rhoz    = 0.5_rp*(rho(i,j,kp)+rho(i,j,k))
          kappasx = 0.5_rp*(kappa(ip,j,k)+kappa(i,j,k))
          kappasy = 0.5_rp*(kappa(i,jp,k)+kappa(i,j,k))
          kappasz = 0.5_rp*(kappa(i,j,kp)+kappa(i,j,k))
          !
          ssx(i,j,k)   = dxi*sigma*kappasx*(psi(ip,j,k)-psi(i,j,k)) 
          ssx_o(i,j,k) = ssx(i,j,k)
          ssy(i,j,k)   = dyi*sigma*kappasy*(psi(i,jp,k)-psi(i,j,k)) 
          ssy_o(i,j,k) = ssy(i,j,k) 
          ssz(i,j,k)   = dzci(k)*sigma*kappasz*(psi(i,j,kp)-psi(i,j,k)) 
          ssz_o(i,j,k) = ssz(i,j,k) 
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine init_surf_tension
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
    !
    allocate(seed(64))
    seed(:) = iseed
    call random_seed( put = seed )
    !
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz*dims(3)
    !
    do k=1,nzg
      do j=1,nyg
        do i=1,nxg
          !
          ii = i-ijk_start(1)
          jj = j-ijk_start(2)
          kk = k-ijk_start(3)
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
    !
    attributes(managed):: p
    !
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz*dims(3)
    !
    !$acc parallel num_gangs(1) vector_length(1) private(h)
    call curand_init(iseed,0,0,h)
    !$acc loop seq
    do k=1,nzg
      do j=1,nyg
        do i=1,nxg
          !
          ii = i-ijk_start(1)
          jj = j-ijk_start(2)
          kk = k-ijk_start(3)
          !
          rn = curand_uniform(h)
          if(ii.ge.1.and.ii.le.nx .and. &
             jj.ge.1.and.jj.le.ny .and. &
             kk.ge.1.and.kk.le.nz ) then
             p(ii,jj,kk) = p(ii,jj,kk) + 2.*(rn-.5)*norm
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
  subroutine set_mean(nx,ny,nz,nh_d,dims,dzlzi,mean,p)
    !
    implicit none
    !
    integer , intent(in   )                      :: nx,ny,nz
    integer , intent(in   )                      :: nh_d
    integer , intent(in   ), dimension(3)        :: dims
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
          meanold = meanold + p(i,j,k)*dzlzi(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !@
    call mpi_allreduce(MPI_IN_PLACE,meanold,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    meanold = meanold/(1._rp*nx*dims(1)*ny*dims(2))
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
  subroutine log_profile(q,n,nh_d,zc,mu2,rho2,lref,uref,p)
    !
    implicit none
    !
    real(rp), intent(in )                     :: q
    integer , intent(in )                     :: n,nh_d
    real(rp), intent(in ), dimension(1-nh_d:) :: zc
    real(rp), intent(in )                     :: mu2,rho2,lref,uref
    real(rp), intent(out), dimension(n)       :: p
    !
    real(rp) :: z,reb,retau ! z/lz and bulk Reynolds number
    integer  :: k
    !
    reb = rho2*uref*lref/mu2
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
  real(rp), intent(in) :: zc
  real(rp) :: fz
    fz = ((1._rp-zc**2)**2)
  end function
  !
  function dfz(zc)
  real(rp), intent(in) :: zc
  real(rp) :: dfz
    dfz = -4._rp*zc*((1._rp-zc**2)**2)
  end function
  !
  function gxy(xc,yc)
  real(rp), intent(in) :: xc,yc
  real(rp) :: gxy
    gxy = yc*exp(-4._rp*(4._rp*xc**2+yc**2))
  end function
  !
  function dgxy(xc,yc)
  real(rp), intent(in) :: xc,yc
  real(rp) :: dgxy
    dgxy = exp(-4._rp*(4._rp*xc**2+yc**2))*(1._rp-8._rp*yc**2)
  end function
  !
end module mod_initflow
