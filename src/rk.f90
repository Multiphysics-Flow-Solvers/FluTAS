!
! SPDX-License-Identifier: MIT
!
module mod_rk
  !
  use mpi
  use mod_types
  use mod_common_mpi, only: ierr
  use mod_param,      only: cbcpre,lx,ly,lz,gacc,is_forced, &
                            bulk_velx,bulk_vely,bulk_velz
  use mod_debug,      only: chkmean
  use mod_mom,        only: momad_xyz
#if defined(_USE_VOF)
  use mod_source,     only: surft_src
#endif
  use mod_source,     only: grav_src,pres_src
#if defined(_TURB_FORCING)
  use mod_source,     only: forc_src
#endif
#if defined(_OPENACC)
  use cudafor
#endif
  !
  implicit none
  !
  private
  public  :: rk
  !
  contains
  !
  subroutine rk(is_first,nx,ny,nz,dims,dxi,dyi,dzi,nh_d,nh_u,dt,dto,rho0i,dzci,dzfi,u,v,w,p,pp, &
                kappa,psi,mu,rho, &
#if defined(_BOUSSINESQ)
                nh_t,tmp, & 
#endif
                dudtrko,dvdtrko,dwdtrko,up,vp,wp,f)
    !
    ! 2nd-order Adams-Bashforth scheme 
    ! for time integration of the momentum equations.
    !
    implicit none
    !
    logical , intent(inout)                                     :: is_first
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   ), dimension(3)                       :: dims
    real(rp), intent(in   )                                     :: dxi,dyi,dzi
    integer , intent(in   )                                     :: nh_d,nh_u
    real(rp), intent(in   )                                     :: dt,dto,rho0i
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzci,dzfi
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in   ), dimension(     0:,     0:,     0:) :: p,pp,kappa,psi,mu,rho
#if defined(_BOUSSINESQ)
    integer , intent(in   )                                     :: nh_t
    real(rp), intent(in   ), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: tmp
#endif
    real(rp), intent(inout), dimension(      :,      :,      :) :: dudtrko,dvdtrko,dwdtrko
    real(rp), intent(out  ), dimension(     0:,     0:,     0:) :: up,vp,wp
    real(rp), intent(out  ), dimension(3), optional             :: f
    !
    real(rp), dimension(nx,ny,nz) :: dudtrk,dvdtrk,dwdtrk
    real(rp) :: factor1,factor2,factor12
    real(rp) :: rho_av
    real(rp) :: meanvel
    integer  :: i,j,k,nh_up
#if defined(_OPENACC)
    attributes(managed) :: dzci,dzfi,u,v,w,p, pp
    attributes(managed) :: dudtrk,dvdtrk,dwdtrk,up,vp,wp
    attributes(managed) :: p,kappa,psi,mu,rho
    attributes(managed) :: dudtrko,dvdtrko,dwdtrko
#if defined(_BOUSSINESQ)
    attributes(managed) :: tmp
#endif
    integer             :: istat
#endif
    !
    nh_up = abs(lbound(up ,1))+1
    !
    factor1  = +( 1._rp+0.5_rp*(dt/dto) )*dt
    factor2  = -(       0.5_rp*(dt/dto) )*dt
    factor12 = factor1 + factor2
    !
    ! calculate average density (to subtract net gravitational force per unit mass)
    ! 
    rho_av = 0._rp
    if(((cbcpre(0,1)//cbcpre(1,1).eq.'PP').and.gacc(1).ne.0._rp).or. &
       ((cbcpre(0,2)//cbcpre(1,2).eq.'PP').and.gacc(2).ne.0._rp).or. &
       ((cbcpre(0,3)//cbcpre(1,3).eq.'PP').and.gacc(3).ne.0._rp)     ) then
#if defined(_OPENACC)
      !$acc parallel loop reduction(+:rho_av)
#endif
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rho_av = rho_av + rho(i,j,k)/(dxi*dyi*dzfi(k))
          enddo
        enddo
      enddo
#if defined(_OPENACC)
      !$acc end parallel loop 
      !@cuf istat=cudaDeviceSynchronize()
#endif
      call MPI_ALLREDUCE(MPI_IN_PLACE,rho_av,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      rho_av = rho_av/(lx*ly*lz)
    endif
    !
    ! compute the r.h.s. of momentum equation
    !
    call momad_xyz(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,mu,rho, &
                   dudtrk,dvdtrk,dwdtrk)  
    !
    if(is_first) then
#if defined(_OPENACC)
      !$acc kernels
#endif
      do k=1,nz
        do j=1,ny
          do i=1,nx
            dudtrko(i,j,k) = dudtrk(i,j,k)
            dvdtrko(i,j,k) = dvdtrk(i,j,k)
            dwdtrko(i,j,k) = dwdtrk(i,j,k)
          enddo
        enddo
      enddo
#if defined(_OPENACC)
      !$acc end kernels 
      !@cuf istat=cudaDeviceSynchronize()
#endif
      is_first = .false.
    endif
    !
    ! compute the velocity prediction 
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor1,factor2,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk,dudtrko,dvdtrko,dwdtrko)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          up(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          wp(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
          !
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end kernels 
    !@cuf istat=cudaDeviceSynchronize()
#else
    !$OMP END PARALLEL DO
#endif
    !
    ! add the incremental pressure
    !
    call pres_src(nx,ny,nz,dxi,dyi,dzi,rho0i,rho,p,pp,dudtrk,dvdtrk,dwdtrk)
    !
    ! add the source terms
    !
    call grav_src(nx,ny,nz,rho,rho_av, &
#if defined(_BOUSSINESQ)
                  nh_t,tmp, & 
#endif
                  dudtrk,dvdtrk,dwdtrk)
#if defined(_USE_VOF)
    call surft_src(nx,ny,nz,dxi,dyi,dzi,kappa,psi,rho,dudtrk,dvdtrk,dwdtrk)
#endif
#if defined(_TURB_FORCING)
    call forc_src(nx,ny,nz,1._rp/dxi,1._rp/dyi,1._rp/dzi,dudtrk,dvdtrk,dwdtrk)
#endif
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          up(i,j,k) = up(i,j,k) + factor12*dudtrk(i,j,k)
          vp(i,j,k) = vp(i,j,k) + factor12*dvdtrk(i,j,k)
          wp(i,j,k) = wp(i,j,k) + factor12*dwdtrk(i,j,k)
          !
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end kernels 
    !!!!!!!@cuf istat=cudaDeviceSynchronize()
#else
    !$OMP END PARALLEL DO
#endif
    !
    ! compute the bulk-velocity forcing
    !
    if(present(f)) then
      f(1:3) = 0._rp
    endif
    if(is_forced(1)) then
      call chkmean(nx,ny,nz,dims,nh_d,nh_up,1._rp/(dzfi*lz),up,meanvel)
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            up(i,j,k) = up(i,j,k) + (bulk_velx-meanvel)
          enddo
        enddo
      enddo
      !$acc end kernels
      if(present(f)) then
        f(1) = bulk_velx-meanvel
      endif
    endif
    if(is_forced(2)) then
      call chkmean(nx,ny,nz,dims,nh_d,nh_up,1._rp/(dzfi*lz),vp,meanvel)
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            vp(i,j,k) = vp(i,j,k) + (bulk_vely-meanvel)
          enddo
        enddo
      enddo
      !$acc end kernels
      if(present(f)) then
        f(2) = bulk_vely-meanvel
      endif
    endif
    if(is_forced(3)) then
      call chkmean(nx,ny,nz,dims,nh_d,nh_up,1._rp/(dzci*lz),wp,meanvel)
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            wp(i,j,k) = wp(i,j,k) + (bulk_velz-meanvel)
          enddo
        enddo
      enddo
      !$acc end kernels
      if(present(f)) then
        f(3) = bulk_velz-meanvel
      endif
    endif
    !
    return
  end subroutine rk
  !
end module mod_rk
