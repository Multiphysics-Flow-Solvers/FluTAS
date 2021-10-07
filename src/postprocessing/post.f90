module mod_post 
  !
  use mod_types
  use mod_common_mpi, only: myid,ierr,ijk_start, &
                            MPI_PROC_NULL,       &
                            top,bottom,          &
                            front,back,          &
                            left,right
  use mpi
  use mod_param,     only: lx,ly,lz
#if defined(_HEAT_TRANSFER)
  use mod_param,     only: kappa1,kappa2,deltaT
#endif
#if defined(_OPENACC)
  use cudafor
  use mod_common_mpi, only: mydev
#endif
  !
  implicit none
  !
  private
  public  :: wall_avg,time_avg,compute_vorticity,mixed_variables
  contains
  !
  subroutine wall_avg(idir,nx,ny,nz,ngx,ngy,ngz,dxi,dyi,dzi,nh_t,ka,tmp,time)
    !
    implicit none
    !
    integer , intent(in )                                      :: idir
    integer , intent(in )                                      :: nx,ny,nz
    integer , intent(in )                                      :: ngx,ngy,ngz
    real(rp), intent(in )                                      :: dxi,dyi,dzi
    integer , intent(in )                                      :: nh_t
    real(rp), intent(in ), dimension(     0:,     0:,     0:)  :: ka
    real(rp), intent(in ), dimension(1-nh_t:,1-nh_t:,1-nh_t:)  :: tmp
    real(rp), intent(in )                                      :: time
    !
    real(rp) :: dtdz,dtdy,dtdx,kap
    real(rp) :: lref,ka_ref
    real(rp) :: nusselt_up,nusselt_down
    integer  :: i,j,k,ip,jp,kp,im,jm,km
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: ka, tmp
#endif
    !
    ka_ref = kappa2
    !
    select case(idir)
    !
    ! Z-ORIENTATION
    !
    case(3)
      !
      lref = lz
      nusselt_up = 0._rp
      !
      if(top.eq.MPI_PROC_NULL) then ! z - top boundary 
        k  = nz+1
        km = nz
        !$acc kernels
        do j=1,ny
          do i=1,nx
            dtdz = (tmp(i,j,k)-tmp(i,j,km))*dzi
            kap  = 0.5_rp*(ka(i,j,k)+ka(i,j,km))
            nusselt_up = nusselt_up + kap*dtdz
          enddo
        enddo
        !$acc end kernels 
        !@cuf istat=cudaDeviceSynchronize()
      endif
      call mpi_allreduce(MPI_IN_PLACE,nusselt_up,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      nusselt_up = (nusselt_up/(1._rp*ngx*ngy))*(lref/ka_ref/deltaT)
      !
      nusselt_down = 0._rp
      !
      if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary 
        k  = 0
        kp = 1
        !$acc kernels
        do j=1,ny 
          do i=1,nx
            dtdz = (tmp(i,j,kp)-tmp(i,j,k))*dzi
            kap  = 0.5_rp*(ka(i,j,k)+ka(i,j,kp))
            nusselt_down = nusselt_down - kap*dtdz
          enddo
        enddo
        !$acc end kernels 
        !@cuf istat=cudaDeviceSynchronize()
      endif
      call mpi_allreduce(MPI_IN_PLACE,nusselt_down,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      nusselt_down = (nusselt_down/(1._rp*ngx*ngy))*(lref/ka_ref/deltaT)
      ! 
      if(myid.eq.0) then
        open(94,file='data/post/wall/nusselt.out',position='append')
        write(94,'(3E15.7)') time,nusselt_down,nusselt_up
        close(94)
      endif
    !
    ! Y-ORIENTATION
    !
    case(2)
      !
      lref = ly
      nusselt_up = 0._rp
      !
      if(back.eq.MPI_PROC_NULL) then ! y - top boundary 
        j  = ny+1
        jm = ny
        !$acc kernels
        do k=1,nz
          do i=1,nx
            dtdy = (tmp(i,j,k)-tmp(i,jm,k))*dyi
            kap  = 0.5_rp*(ka(i,j,k)+ka(i,jm,k))
            nusselt_up = nusselt_up + kap*dtdy
          enddo
        enddo
        !$acc end kernels 
        !@cuf istat=cudaDeviceSynchronize()
      endif
      call mpi_allreduce(MPI_IN_PLACE,nusselt_up,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      nusselt_up = (nusselt_up/(1._rp*ngx*ngz))*(lref/ka_ref/deltaT)
      !
      nusselt_down = 0._rp
      !
      if(front.eq.MPI_PROC_NULL) then ! y - bottom boundary 
        j  = 0
        jp = 1
        !$acc kernels
        do k=1,nz 
          do i=1,nx
            dtdy = (tmp(i,jp,k)-tmp(i,j,k))*dyi
            kap  = 0.5_rp*(ka(i,j,k)+ka(i,jp,k))
            nusselt_down = nusselt_down - kap*dtdy
          enddo
        enddo
        !$acc end kernels 
        !@cuf istat=cudaDeviceSynchronize()
      endif
      call mpi_allreduce(MPI_IN_PLACE,nusselt_down,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      nusselt_down = (nusselt_down/(1._rp*ngx*ngz))*(lref/ka_ref/deltaT)
      ! 
      if(myid.eq.0) then
        open(94,file='data/post/wall/nusselt.out',position='append')
        write(94,'(3E15.7)') time,nusselt_down,nusselt_up
        close(94)
      endif
    !
    ! X-ORIENTATION
    !
    case(1)
      !
      lref = lx
      nusselt_up = 0._rp
      !
      if(right.eq.MPI_PROC_NULL) then ! x - top boundary 
        i  = nx+1
        im = nx
        !$acc kernels
        do k=1,nz
          do j=1,ny
            dtdx = (tmp(i,j,k)-tmp(im,j,k))*dxi
            kap  = 0.5_rp*(ka(i,j,k)+ka(im,j,k))
            nusselt_up = nusselt_up + kap*dtdx
          enddo
        enddo
        !$acc end kernels 
        !@cuf istat=cudaDeviceSynchronize()
      endif
      call mpi_allreduce(MPI_IN_PLACE,nusselt_up,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      nusselt_up = (nusselt_up/(1._rp*ngy*ngz))*(lref/ka_ref/deltaT)
      !
      nusselt_down = 0._rp
      !
      if(left.eq.MPI_PROC_NULL) then ! x - bottom boundary 
        i  = 0
        ip = 1
        !$acc kernels
        do k=1,nz 
          do j=1,ny
            dtdx = (tmp(ip,j,k)-tmp(i,j,k))*dxi
            kap  = 0.5_rp*(ka(i,j,k)+ka(ip,j,k))
            nusselt_down = nusselt_down - kap*dtdx
          enddo
        enddo
        !$acc end kernels 
        !@cuf istat=cudaDeviceSynchronize()
      endif
      call mpi_allreduce(MPI_IN_PLACE,nusselt_down,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      nusselt_down = (nusselt_down/(1._rp*ngy*ngz))*(lref/ka_ref/deltaT)
      ! 
      if(myid.eq.0) then
        open(94,file='data/post/wall/nusselt.out',position='append')
        write(94,'(3E15.7)') time,nusselt_down,nusselt_up
        close(94)
      endif
    end select
    !
    return
  end subroutine wall_avg
  !
  subroutine time_avg(idir,do_avg,do_favre,fname,n,ng,istep,istep_av,iout1d, &
                      nh_d,nh_v,nh_p,psi,rho_p,p,pout1,pout2,pvol1,pvol2)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions
    !
    ! idir     -> the direction normal to the averaging plane
    ! do_avg   -> do or not averaging
    ! do_favre -> Favre or regular averaging
    ! fname    -> name of the file
    ! n        -> size of the input array
    ! ng       -> total size of computational domain
    ! istep    -> current time step
    ! istep_av -> size of statistical sample
    ! iout1d   -> print file every iout1d time steps
    ! idir     -> direction of the profile (all other directions are averaged)
    ! nh_      -> halos point
    ! psi      -> 3D vof field (0 --> liquid, 1 --> gas)
    ! rho_p    -> density of the phase (for compressible phase)
    ! p        -> 3D input scalar field to be averaged
    ! pout1    ->  first order time statistics of plane averaged field (mean)
    ! pout2    -> second order time statistics of plane averaged field (rms)
    ! pvol     ->  first order time statistics of volume averaged field (mean)
    ! pvol2    -> second order time statistics of volume averaged field (rms)
    !
    use mod_param, only: dl
    !
    implicit none
    !
    integer         , intent(in   )                                     :: idir
    logical         , intent(in   )                                     :: do_avg
    logical         , intent(in   )                                     :: do_favre
    character(len=*), intent(in   )                                     :: fname
    integer         , intent(in   ), dimension(3)                       :: n,ng
    integer         , intent(in   )                                     :: istep,istep_av,iout1d
    integer         , intent(in   )                                     :: nh_d,nh_v,nh_p
    real(rp)        , intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: psi
    real(rp)        , intent(in   ), dimension(0     :,0     :,0     :) :: rho_p
    real(rp)        , intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    real(rp)        , intent(inout), dimension(1:)                      :: pout1,pout2
    real(rp)        , intent(inout)                                     :: pvol1,pvol2
    !
    real(rp), allocatable, dimension(:) :: p1d1,p1d2,rhop
    real(rp) :: factor,factor2,p_dl12
    integer  :: i,j,k,mm
    integer  :: nx,ny,nz,ngx,ngy,ngz,ng_idir
    integer  :: start
    integer  :: iunit
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: p,psi,rho_p,p1d1,p1d2,rhop,pout1,pout2
#endif
    !
    p_dl12 = dl(1)*dl(2)
    !
    nx      = n(1)
    ny      = n(2)
    nz      = n(3)
    ngx     = ng(1)
    ngy     = ng(2)
    ngz     = ng(3)
    ng_idir = ng(idir)
    start   = ijk_start(idir) 
    !
    allocate(p1d1(ng(idir)),p1d2(ng(idir)),rhop(ng(idir)))
    !
    iunit   = 10
    factor  = 1._rp*istep_av
    factor2 = 1._rp*ng_idir/(1._rp*ngx*ngy*ngz)
    !
    if(istep_av.eq.1) then  
      do mm=1,ng_idir
        pout1(mm) = 0._rp
        pout2(mm) = 0._rp
      enddo
    endif
    !
    !$acc kernels
    do mm=1,ng_idir
      p1d1(mm)  = 0._rp
      p1d2(mm)  = 0._rp
      rhop(mm)  = 0._rp
    enddo
    !$acc end kernels 
    !
    !   
    if(do_favre) then
      !
      ! Density-based averaging (Favre)
      !
      select case(idir)
      case(1)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + i
              !
              rhop(mm) = rhop(mm) + p_dl12*rho_p(i,j,k)
              p1d1(mm) = p1d1(mm) + p_dl12*rho_p(i,j,k)*(1._rp-psi(i,j,k))*p(i,j,k)     
              p1d2(mm) = p1d2(mm) + p_dl12*rho_p(i,j,k)*(1._rp-psi(i,j,k))*p(i,j,k)**2 
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
        !@cuf istat=cudaDeviceSynchronize()
      case(2)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + j
              !
              rhop(mm) = rhop(mm) + p_dl12*rho_p(i,j,k)
              p1d1(mm) = p1d1(mm) + p_dl12*rho_p(i,j,k)*(1._rp-psi(i,j,k))*p(i,j,k)     
              p1d2(mm) = p1d2(mm) + p_dl12*rho_p(i,j,k)*(1._rp-psi(i,j,k))*p(i,j,k)**2 
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
        !@cuf istat=cudaDeviceSynchronize()
      case(3)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + k
              !
              rhop(mm) = rhop(mm) + p_dl12*rho_p(i,j,k)
              p1d1(mm) = p1d1(mm) + p_dl12*rho_p(i,j,k)*(1._rp-psi(i,j,k))*p(i,j,k)     
              p1d2(mm) = p1d2(mm) + p_dl12*rho_p(i,j,k)*(1._rp-psi(i,j,k))*p(i,j,k)**2 
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
        !@cuf istat=cudaDeviceSynchronize()
      end select
      !
      call mpi_allreduce(MPI_IN_PLACE,rhop(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !
      !$acc kernels
      do mm=1,ng_idir
        p1d1(mm) = p1d1(mm)/(rhop(mm))
        p1d2(mm) = p1d2(mm)/(rhop(mm))
      enddo
      !$acc end kernels 
      !@cuf istat=cudaDeviceSynchronize()
      !
    else
      !
      ! Volumetric-based averaging (no Favre)
      !
      select case(idir)
      case(1)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + i
              !
              p1d1(mm) = p1d1(mm) + (1._rp-psi(i,j,k))*p(i,j,k)
              p1d2(mm) = p1d2(mm) + (1._rp-psi(i,j,k))*p(i,j,k)**2
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
        !@cuf istat=cudaDeviceSynchronize()
      case(2)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + j
              !
              p1d1(mm) = p1d1(mm) + (1._rp-psi(i,j,k))*p(i,j,k)
              p1d2(mm) = p1d2(mm) + (1._rp-psi(i,j,k))*p(i,j,k)**2
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
        !@cuf istat=cudaDeviceSynchronize()
      case(3)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + k
              !
              p1d1(mm) = p1d1(mm) + (1._rp-psi(i,j,k))*p(i,j,k)
              p1d2(mm) = p1d2(mm) + (1._rp-psi(i,j,k))*p(i,j,k)**2
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
        !@cuf istat=cudaDeviceSynchronize()
      end select
      !
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !
      !$acc kernels
      do mm=1,ng_idir
        p1d1(mm) = p1d1(mm)*factor2
        p1d2(mm) = p1d2(mm)*factor2
      enddo
      !$acc end kernels 
      !@cuf istat=cudaDeviceSynchronize()
      !
    endif
    !
    ! decide or not to averaging 
    !
    if(.not.do_avg) then
      !$acc kernels
      do mm=1,ng_idir
        pout1(mm) = p1d1(mm)
        pout2(mm) = p1d2(mm)
      enddo
      !$acc end kernels 
      !@cuf istat=cudaDeviceSynchronize()
    else
      !$acc kernels
      do mm=1,ng_idir
        pout1(mm) = ((factor-1._rp)*pout1(mm)+p1d1(mm))/factor
        pout2(mm) = ((factor-1._rp)*pout2(mm)+p1d2(mm))/factor
      enddo
      !$acc end kernels 
      !@cuf istat=cudaDeviceSynchronize()
    endif
    pvol1 = sum(pout1)
    pvol2 = sum(pout2)
    !
    ! print 
    !  note: we put this condition on iout1d in order to ensure that the 
    !        averaging frequency (iout_av) is indenpendent 
    !        of the print frequency of the files (iout1d)
    !
    if(mod(istep,iout1d).eq.0) then
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do mm=1,ng_idir
          write(iunit,'(3E15.7)') (mm-0.5_rp)*dl(idir),pout1(mm),pout2(mm)
        enddo
        close(iunit)
      endif
    endif
    !
    return
  end subroutine time_avg
  !
  !
  subroutine compute_vorticity(nx,ny,nz,dxi,dyi,dzi,nh_u,v,w,vor)
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_u
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: v,w
    real(rp), intent(out), dimension(     1:,     1:,     1:) :: vor
    !
    real(rp) :: vp,vm,wp,wm
    integer  :: i,j,k,im,jm,km,ip,jp,kp
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: v, w
#endif
    !
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          kp = k+1
          km = k-1
          jp = j+1
          jm = j-1
          ip = i+1
          im = i-1
          !
          vp = 0.25*(v(i,j,k)+v(i,jm,k)+v(i,j,kp)+v(i,jm,kp))
          vm = 0.25*(v(i,j,k)+v(i,jm,k)+v(i,j,km)+v(i,jm,km))
          wp = 0.25*(w(i,j,k)+w(i,j,km)+w(i,jp,k)+w(i,jp,km))
          wm = 0.25*(w(i,j,k)+w(i,j,km)+w(i,jm,k)+w(i,jm,km))
          !
          vor(i,j,k) = (wp-wm)*dyi-(vp-vm)*dzi
          !
        enddo
      enddo
    enddo
    !$acc end kernels 
    !@cuf istat=cudaDeviceSynchronize()
    !
    return
  end subroutine compute_vorticity
  !
  subroutine mixed_variables(nx,ny,nz,dxi,dyi,dzi,nh_u,nh_s1, & 
                             u,v,w,s1, &
                             us1,vs1,ws1,uv,vw,wu)
    !
    implicit none
    !
    integer , intent(in )                                        :: nx,ny,nz
    real(rp), intent(in )                                        :: dxi,dyi,dzi
    integer , intent(in )                                        :: nh_u,nh_s1!,nh_s2
    real(rp), intent(in ), dimension(1-nh_u :,1-nh_u :,1-nh_u :) :: u,v,w
    real(rp), intent(in ), dimension(1-nh_s1:,1-nh_s1:,1-nh_s1:) :: s1 ! generic scalar
    real(rp), intent(out), dimension(     1:,     1:,     1:)    :: us1,vs1,ws1
    real(rp), intent(out), dimension(     1:,     1:,     1:)    :: uv ,vw ,wu
    !
    integer :: i,j,k,im,jm,km
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: u,v,w,s1
#endif
    !
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          km = k-1
          jm = j-1
          im = i-1
          !
          uv(i,j,k)  = 0.25_rp*(u(i,j,k)+u(im,j,k))*(v(i,j,k)+v(i,jm,k))
          vw(i,j,k)  = 0.25_rp*(v(i,j,k)+v(i,jm,k))*(w(i,j,k)+w(i,j,km))
          wu(i,j,k)  = 0.25_rp*(w(i,j,k)+w(i,j,km))*(u(i,j,k)+u(im,j,k))
          !
          us1(i,j,k) = 0.5_rp*(u(i,j,k)+u(im,j,k))*s1(i,j,k)
          vs1(i,j,k) = 0.5_rp*(v(i,j,k)+v(i,jm,k))*s1(i,j,k)
          ws1(i,j,k) = 0.5_rp*(w(i,j,k)+w(i,j,km))*s1(i,j,k)
          !
        enddo
      enddo
    enddo
    !$acc end kernels 
    !@cuf istat=cudaDeviceSynchronize()
    !
    return
  end subroutine mixed_variables
  ! 
end module mod_post
