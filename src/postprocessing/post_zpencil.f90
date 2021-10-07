module mod_post 
  !
  use mod_types
  use mod_common_mpi, only: myid, comm_Cart,ierr,ijk_start
  use mpi
  use mod_bound, only: boundp
!   use mod_param
!   use mod_common_mpi
!   use decomp_2d
!   use mod_phase_indicator
  implicit none
  private
  public  time_avg, wall_avg, compute_vorticity, mixed_variables
contains

subroutine time_avg(fname,n,ng,istep,istep_av,iout1d,idir,z,dzlzi,psi,p,pout1,pout2,pvol1,pvol2)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions
    !
    ! fname    -> name of the file
    ! n        -> size of the input array
    ! ng       -> total size of computational domain
    ! istep    -> current time step
    ! istep_av -> size of statistical sample
    ! iout1d   -> print file every iout1d time steps
    ! idir     -> direction of the profile (all other directions are averaged)
    ! z        -> z coordinate (grid is non-uniform in z)
    ! dzlzi    -> dz/lz weight of a grid cell for averaging over z
    ! p        -> 3D vof field (0 --> liquid, 1 --> gas)
    ! p        -> 3D input scalar field to be averaged
    ! pout1    ->  first order time statistics of plane averaged field (mean)
    ! pout2    -> second order time statistics of plane averaged field (rms)
    ! pvol     ->  first order time statistics of volume averaged field (mean)
    ! pout2    -> second order time statistics of volume averaged field (rms)
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: n,ng
    integer , intent(in) :: istep,istep_av,iout1d
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(0:) :: z,dzlzi
    real(rp), intent(in), dimension(1:,1:,1:) :: psi
    real(rp), intent(in), dimension(1:,1:,1:) :: p
    real(rp), intent(inout), dimension(1:) :: pout1,pout2
    real(rp), intent(inout) :: pvol1,pvol2
    real(rp), allocatable, dimension(:) :: p1d1,p1d2
    integer :: i,j,k,ii,jj
    integer :: iunit
    real(rp) :: factor
    !
    iunit  = 61
    factor = REAL(istep_av)
    select case(idir)
    case(3)
      allocate(p1d1(n(3)))
      allocate(p1d2(n(3)))
      do k=1,ng(3)
        p1d1(k) = 0.
        p1d2(k) = 0.
        do j=1,n(2)
          do i=1,n(1)
            p1d1(k) = p1d1(k) + (1.d0-psi(i,j,k))*p(i,j,k)
            p1d2(k) = p1d2(k) + (1.d0-psi(i,j,k))*p(i,j,k)**2
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d1(:)  = p1d1(:)/(1.*ng(1)*ng(2))
      p1d2(:)  = p1d2(:)/(1.*ng(1)*ng(2))
      pout1(:) = ((factor-1.)*pout1(:)+p1d1(:))/factor
      pout2(:) = ((factor-1.)*pout2(:)+p1d2(:))/factor
      pvol1    = sum(pout1)
      pvol2    = sum(pout2)
      if(mod(istep,iout1d).eq.0) then
        if(myid.eq.0) then
          open(unit=iunit,file=fname)
          do k=1,n(3)
            write(iunit,'(3E15.7)') z(k),pout1(k),pout2(k)
          enddo
          close(iunit)
        endif
      endif
    case(2)
      allocate(p1d1(ng(2)))
      allocate(p1d2(ng(2)))
      p1d1(:) = 0.
      p1d2(:) = 0.
      do j=1,n(2)
        jj = ijk_start(2)+j
        p1d1(jj) = 0.
        p1d2(jj) = 0.
        do k=1,n(3)
          do i=1,n(1)
            p1d1(jj) = p1d1(jj) + (1.d0-psi(i,j,k))*p(i,j,k)*dzlzi(k)
            p1d2(jj) = p1d2(jj) + (1.d0-psi(i,j,k))*p(i,j,k)*p(i,j,k)*dzlzi(k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),ng(2),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),ng(2),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d1(:) = p1d1(:)/(1.*ng(1))
      p1d2(:) = p1d2(:)/(1.*ng(1))
      pout1(:) = ((factor-1.)*pout1(:)+p1d1(:))/factor
      pout2(:) = ((factor-1.)*pout2(:)+p1d2(:))/factor
      pvol1    = sum(pout1)
      pvol2    = sum(pout2)
      if(mod(istep,iout1d).eq.0) then
        if(myid.eq.0) then
          open(unit=iunit,file=fname)
          do j=1,ng(2)
            write(iunit,'(3E15.7)') (1.*j-.5)/(1.*ng(2)),pout1(j),pout2(j)
          enddo
          close(iunit)
        endif
      endif
    case(1)
      allocate(p1d1(ng(1)))
      allocate(p1d2(ng(1)))
      p1d1(:) = 0.
      p1d2(:) = 0.
      do i=1,n(1)
        ii = ijk_start(1)+i
        p1d1(i) = 0.
        p1d2(i) = 0.
        do k=1,n(3)
          do j=1,n(2)
            p1d1(ii) = p1d1(ii) + (1.d0-psi(i,j,k))*p(i,j,k)*dzlzi(k)
            p1d2(ii) = p1d2(ii) + (1.d0-psi(i,j,k))*p(i,j,k)*p(i,j,k)*dzlzi(k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),ng(1),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),ng(1),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d1(:) = p1d1(:)/(1.*ng(2))
      p1d2(:) = p1d2(:)/(1.*ng(2))
      pout1(:) = ((factor-1.)*pout1(:)+p1d1(:))/factor
      pout2(:) = ((factor-1.)*pout2(:)+p1d2(:))/factor
      pvol1    = sum(pout1)
      pvol2    = sum(pout2)
      if(mod(istep,iout1d).eq.0) then
        if(myid.eq.0) then
          open(unit=iunit,file=fname)
          do i=1,ng(1)
            write(iunit,'(3E15.7)') (1.*i-.5)/(1.*n(1)),pout1(i),pout2(i)
          enddo
          close(iunit)
        endif
      endif
    end select
  return
end subroutine time_avg

subroutine wall_avg(n,ng,dli,mu,u,w, &
#if defined(_HEAT_TRANSFER)
                              ka,tmp,&
#endif
                                time )  
#if defined(_HEAT_TRANSFER)
  use mod_param, only: kappa1,kappa2,lref,deltaT
#endif
  implicit none
  integer , intent(in), dimension(3)                  :: n
  integer , intent(in), dimension(3)                  :: ng
  real(rp), intent(in), dimension(3)                  :: dli
  real(rp), intent(in), dimension(0:,0:,0:)           :: mu
  real(rp), intent(in), dimension(0:,0:,0:)           :: u,w
#if defined(_HEAT_TRANSFER)
  real(rp), intent(in), dimension(0:,0:,0:)           :: ka
  real(rp), intent(in), dimension(0:,0:,0:)           :: tmp
#endif
  real(rp), intent(in)                                :: time
  integer  :: i,j,k
  integer  :: ip,km,kp
  real(rp) :: dudzm,dudzp,dwdxm,dwdxp,muzm,muzp
  real(rp) :: tau_up,tau_down
#if defined(_HEAT_TRANSFER)
  real(rp) :: dtdz,kap
  real(rp) :: nusselt_up,nusselt_down
  real(rp) :: ka_ref
#endif

    tau_up = 0.d0
    k  = n(3)+1
    km = n(3)
    do j=1,n(2)
      do i=1,n(1)
        ip = i+1
        dudzm = (u(i ,j,k ) - u(i,j,km))*dli(3)
        dwdxm = (w(ip,j,km) - w(i,j,km))*dli(1)
        muzm  = 0.25d0*(mu(i,j,k)+mu(i,j,km)+mu(ip,j,km)+mu(ip,j,k))
        tau_up = tau_up - muzm*(dudzm+dwdxm)
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,tau_up,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    tau_up = tau_up/(1.d0*ng(1)*ng(2))
    !
    tau_down = 0.d0
    k  = 0
    kp = 1
    do j=1,n(2) 
      do i=1,n(1)
        ip = i+1
        dudzp = (u(i ,j,kp) - u(i,j,k))*dli(3)
        dwdxp = (w(ip,j,k ) - w(i,j,k))*dli(1)
        muzp  = 0.25d0*(mu(i,j,k)+mu(i,j,kp)+mu(ip,j,kp)+mu(ip,j,k))
        tau_down = tau_down + muzp*(dudzp+dwdxp)
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,tau_down,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    tau_down = tau_down/(1.d0*ng(1)*ng(2))
  
    if(myid .eq. 0 ) then
       open(93,file='data/post/wall/shear.out',position='append')
       write(93,'(3E15.7)') time,tau_down,tau_up
       close(93)
    endif
    !
    !
    !
#if defined(_HEAT_TRANSFER)
    ka_ref = kappa2
    nusselt_up = 0.d0
    k  = n(3)+1
    km = n(3)
    do j=1,n(2)
      do i=1,n(1)
        dtdz = (tmp(i ,j,k ) - tmp(i,j,km))*dli(3)
        kap  = 0.5d0*(ka(i,j,k)+ka(i,j,km))
        nusselt_up = nusselt_up + kap*dtdz
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,nusselt_up,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    nusselt_up = (nusselt_up/(1.d0*ng(1)*ng(2)))*(lref/ka_ref/deltaT)
    !
    nusselt_down = 0.d0
    k  = 0
    kp = 1
    do j=1,n(2) 
      do i=1,n(1)
        dtdz  = (tmp(i ,j,kp) - tmp(i,j,k))*dli(3)
        kap   = 0.5d0*(ka(i,j,k)+ka(i,j,kp))
        nusselt_down = nusselt_down - kap*dtdz
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,nusselt_down,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    nusselt_down = (nusselt_down/(1.d0*ng(1)*ng(2)))*(lref/ka_ref/deltaT)
    
    if(myid .eq. 0 ) then
       open(94,file='data/post/wall/nusselt.out',position='append')
       write(94,'(3E15.7)') time,nusselt_down,nusselt_up
       close(94)
    endif
#endif
  return
end subroutine wall_avg

subroutine compute_vorticity(n,dli,v,w,vor)
  implicit none
  integer , intent(in), dimension(3)             :: n
  real(rp), intent(in), dimension(3)             :: dli
  real(rp), intent(in), dimension(0:,0:,0:)      :: v,w
  real(rp), intent(out), dimension(1:,1:,1:)     :: vor
  integer :: i,j,k,im,jm,km,ip,jp,kp
  real(rp):: vp,vm,wp,wm
  !
  do k=1,n(3)
    kp = k+1
    km = k-1
    do j=1,n(2)
      jp = j+1
      jm = j-1
      do i=1,n(1)
        ip = i+1
        im = i-1
        !
        vp = 0.25*(v(i,j,k)+v(i,jm,k)+v(i,j,kp)+v(i,jm,kp))
        vm = 0.25*(v(i,j,k)+v(i,jm,k)+v(i,j,km)+v(i,jm,km))
        wp = 0.25*(w(i,j,k)+w(i,j,km)+w(i,jp,k)+w(i,jp,km))
        wm = 0.25*(w(i,j,k)+w(i,j,km)+w(i,jm,k)+w(i,jm,km))
        vor(i,j,k)=(wp-wm)*dli(2)-(vp-vm)*dli(3)
      enddo
    enddo
  enddo
  return
end subroutine compute_vorticity

subroutine mixed_variables(n,dli,u,v,w,   &
#if defined(_HEAT_TRANSFER)
                           tmp,           &
                           utmp,vtmp,wtmp, &
#endif
                           uv,vw,wu)
  implicit none
  integer , intent(in), dimension(3)             :: n
  real(rp), intent(in), dimension(3)             :: dli
  real(rp), intent(in), dimension(0:,0:,0:)      :: u,v,w
#if defined(_HEAT_TRANSFER)
  real(rp), intent(in), dimension(0:,0:,0:)      :: tmp
  real(rp), intent(out), dimension(1:,1:,1:)     :: utmp,vtmp,wtmp
#endif
  real(rp), intent(out), dimension(1:,1:,1:)     :: uv,vw,wu
  integer :: i,j,k,im,jm,km
  !
  do k=1,n(3)
    km = k-1
    do j=1,n(2)
      jm = j-1
      do i=1,n(1)
        im = i-1
        !
        uv(i,j,k) = 0.25*(u(i,j,k)+u(im,j,k))*(v(i,j,k)+v(i,jm,k))
        vw(i,j,k) = 0.25*(v(i,j,k)+v(i,jm,k))*(w(i,j,k)+w(i,j,km))
        wu(i,j,k) = 0.25*(w(i,j,k)+w(i,j,km))*(u(i,j,k)+u(im,j,k))
#if defined(_HEAT_TRANSFER)
        utmp(i,j,k) = 0.5*(u(i,j,k)+u(im,j,k))*tmp(i,j,k)
        vtmp(i,j,k) = 0.5*(v(i,j,k)+v(i,jm,k))*tmp(i,j,k)
        wtmp(i,j,k) = 0.5*(w(i,j,k)+w(i,j,km))*tmp(i,j,k)
#endif
      enddo
    enddo
  enddo
  return
end subroutine mixed_variables
  !
end module mod_post
