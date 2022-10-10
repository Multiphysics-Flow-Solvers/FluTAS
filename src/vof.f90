!
! SPDX-License-Identifier: MIT
!
module mod_vof
  !
  ! based on:
  ! Ii, Satoshi, et al. "An interface capturing method with a continuous function: 
  ! The THINC method with multi-dimensional reconstruction." 
  ! Journal of Computational Physics 231.5 (2012): 2328-2358.
  !
  use mod_types, only: rp
  !@cuf use cudafor
  !
  implicit none
  !
  real(rp), parameter :: b_th  = 2.0_rp, qu_th  = 1.0_rp, limit = 10._rp**(-8)
  real(rp), parameter :: b_thi = 1.0_rp/b_th
  real(rp), parameter :: small = 10._rp**(-16)
  !
  real(rp), parameter :: rmm = 0.5_rp*(1.0_rp-1.0_rp/sqrt(3.0_rp)), &
                         rpp = 0.5_rp*(1.0_rp+1.0_rp/sqrt(3.0_rp))
  !
  ! Those are constants 
  real(rp), parameter, dimension(3,3) :: xmm = reshape((/0.0_rp,rmm,rmm, &
                                                        rmm,0.0_rp,rmm, &
                                                        rmm,rmm,0.0_rp/),shape(xmm))
  real(rp), parameter, dimension(3,3) :: xmp = reshape((/0.0_rp,rmm,rpp, &
                                                        rmm,0.0_rp,rpp, &
                                                        rmm,rpp,0.0_rp/),shape(xmp))
  real(rp), parameter, dimension(3,3) :: xpm = reshape((/0.0_rp,rpp,rmm, &
                                                        rpp,0.0_rp,rmm, &
                                                        rpp,rmm,0.0_rp/),shape(xpm))
  real(rp), parameter, dimension(3,3) :: xpp = reshape((/0.0_rp,rpp,rpp, &
                                                        rpp,0.0_rp,rpp, &
                                                        rpp,rpp,0.0_rp/),shape(xpp))
  !@cuf attributes(managed) :: xmm, xmp, xpm, xpp
  !
  real(rp), allocatable, dimension(:,:,:) :: dvof1,dvof2,flux
  !@cuf attributes(managed) :: dvof1,dvof2,flux
  !
  private
  public  :: advvof,update_vof,update_property,initvof
  !
  contains
  !
  subroutine advvof(n,dli,dt,halo,nh_d,nh_u,dzc,dzf,ug,vg,wg,vof,nor,cur,kappa,d_thinc)
    !
    use mod_param, only: cbcvof,bcvof
    use mod_bound, only: boundp
    !
    ! VoF advection
    !
    implicit none
    !
    integer , intent(in   ), dimension(3)                       :: n
    real(rp), intent(in   ), dimension(3)                       :: dli
    real(rp), intent(in   )                                     :: dt
    integer , intent(in   ), dimension(3)                       :: halo
    integer , intent(in   )                                     :: nh_d,nh_u
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: ug,vg,wg ! interface velocity
    real(rp), intent(inout), dimension(0:,0:,0:)                :: vof
    real(rp), intent(inout), dimension(0:,0:,0:,1:)             :: nor,cur
    real(rp), intent(inout), dimension(0:,0:,0:)                :: kappa    ! curvature
    real(rp), intent(inout), dimension(0:,0:,0:)                :: d_thinc
    !
    real(rp), dimension(3) :: dl
    real(rp) :: dli1, dli2, dli3
    integer :: i,j,k,im,jm,km
    integer :: n1, n2, n3
    !@cuf attributes(managed) :: nor, cur, kappa, d_thinc, vof, ug, vg, wg, dzc, dzf
    !
    ! we allocate dvof1, dvof2 and flux only once
    !
    if(.not.allocated(dvof1)) allocate(dvof1(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    if(.not.allocated(dvof2)) allocate(dvof2(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    if(.not.allocated(flux )) allocate(flux (0:n(1)+0,0:n(2)+0,0:n(3)+0))
    !
    dl(:) = dli(:)**(-1)
    dli1 = dli(1)
    dli2 = dli(2)
    dli3 = dli(3)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)
    !
    ! TODO: prefetcing
    !
    ! flux in x
    !
#if !defined(_TWOD)
    call cmpt_vof_flux(n(1), n(2), n(3),dli(1),dt,nh_u,vof,nor,cur,d_thinc,1,ug,flux) 
#endif
    !
    !$acc kernels
    do k=1,n3
      do j=1,n2
        do i=1,n1
          im = i-1
          !
#if defined(_TWOD)
          dvof1(i,j,k) = vof(i,j,k)
#else
          dvof1(i,j,k) = (vof(i,j,k)-(flux(i,j,k)-flux(im,j,k))*dli1)/(1.0_rp-dt*dli1*(ug(i,j,k)-ug(im,j,k)))
#endif
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    ! update vof
    !
#if !defined(_TWOD)
    call clip_vof(n,dvof1)
    call boundp(cbcvof,n,bcvof,nh_d,+1,halo,dl,dzc,dzf,dvof1)
    call update_vof(n,dli,nh_d,dzc,dzf,+1,halo,dvof1,nor,cur,kappa,d_thinc)
#endif
    !
    ! flux in y
    !
    call cmpt_vof_flux(n(1), n(2), n(3),dli(2),dt,nh_u,dvof1,nor,cur,d_thinc,2,vg,flux) 
    !
    !$acc kernels
    do k=1,n3
      do j=1,n2
        do i=1,n1
          jm = j-1
          !
          dvof2(i,j,k) = (dvof1(i,j,k)-(flux(i,j,k)-flux(i,jm,k))*dli2)/(1.0_rp-dt*dli2*(vg(i,j,k)-vg(i,jm,k)))
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    ! update vof
    !
    call clip_vof(n,dvof2)
    call boundp(cbcvof,n,bcvof,nh_d,+1,halo,dl,dzc,dzf,dvof2)
    call update_vof(n,dli,nh_d,dzc,dzf,+1,halo,dvof2,nor,cur,kappa,d_thinc)
    !
    ! flux in z
    !
    call cmpt_vof_flux(n(1), n(2), n(3),dli(3),dt,nh_u,dvof2,nor,cur,d_thinc,3,wg,flux) 
    !
    !$acc kernels
    do k=1,n3
      do j=1,n2
        do i=1,n1
          km = k-1
          !
          vof(i,j,k) = (dvof2(i,j,k)-(flux(i,j,k)-flux(i,j,km))*dli3)/(1.0_rp-dt*dli3*(wg(i,j,k)-wg(i,j,km)))
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    ! update vof
    !
    call clip_vof(n,vof)
    !call boundp(cbcvof,n,bcvof,nh_d,+1,halo,dl,dzc,dzf,vof) ! it can be skipped
    !call update_vof(n,dli,nh_d,dzc,dzf,+1,halo,vof,nor,cur,kappa,d_thinc) ! it can be skipped
    !
    ! divergence correction step
    !
    !$acc kernels
    do k=1,n3
      do j=1,n2
        do i=1,n1
          km = k-1
          jm = j-1
          im = i-1
          !
          vof(i,j,k) = vof(i,j,k) - dt*( dvof1(i,j,k)*dli1*(ug(i,j,k)-ug(im,j,k)) + &
                                         dvof2(i,j,k)*dli2*(vg(i,j,k)-vg(i,jm,k)) + &
                                           vof(i,j,k)*dli3*(wg(i,j,k)-wg(i,j,km)) )
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    ! update vof
    !
    call clip_vof(n,vof)
    call boundp(cbcvof,n,bcvof,nh_d,+1,halo,dl,dzc,dzf,vof)
    call update_vof(n,dli,nh_d,dzc,dzf,+1,halo,vof,nor,cur,kappa,d_thinc)
    !
    return
  end subroutine advvof
  !
  subroutine update_vof(n,dli,nh_d,dzc,dzf,nh_p,halo,vof,nor,cur,kappa,d_thinc)
    !
    use mod_param, only: cbcvof,bcvof
    use mod_bound, only: boundp
    !
    implicit none
    !
    integer , intent(in ), dimension(3)           :: n
    real(rp), intent(in ), dimension(3)           :: dli
    integer , intent(in )                         :: nh_d
    real(rp), intent(in ), dimension(1-nh_d:)     :: dzc,dzf
    integer , intent(in )                         :: nh_p
    integer , intent(in ), dimension(3)           :: halo
    real(rp), intent(in ), dimension(0:,0:,0:   ) :: vof
    real(rp), intent(out), dimension(0:,0:,0:,1:) :: nor
    real(rp), intent(out), dimension(0:,0:,0:,1:) :: cur
    real(rp), intent(out), dimension(0:,0:,0:   ) :: kappa
    real(rp), intent(out), dimension(0:,0:,0:   ) :: d_thinc
    !
    real(rp), dimension(3) :: dl
    integer :: p,i,j,k
    integer :: n1, n2, n3
    !@cuf attributes(managed) :: nor, cur, kappa, d_thinc, vof, dzc, dzf
    !
    dl(:) = dli(:)**(-1)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)
    !
    !call cmpt_nor_curv_1o(n(1),n(2),n(3),dli,vof,nor,cur,kappa)
    call cmpt_nor_curv_2o(n(1),n(2),n(3),dli,vof,nor,cur,kappa)
    !
    do p=1,3
      call boundp(cbcvof,n,bcvof,nh_d,nh_p,halo,dl,dzc,dzf,nor(:,:,:,p))
    enddo
    do p=1,6
      call boundp(cbcvof,n,bcvof,nh_d,nh_p,halo,dl,dzc,dzf,cur(:,:,:,p))
    enddo
    call boundp(cbcvof,n,bcvof,nh_d,nh_p,halo,dl,dzc,dzf,kappa)
    !
    call cmpt_d_thinc(n1, n2, n3,nor,cur,vof,d_thinc)
    call boundp(cbcvof,n,bcvof,nh_d,nh_p,halo,dl,dzc,dzf,d_thinc)
    !
    return
  end subroutine update_vof
  !
  subroutine update_property(n,prop12,vof,prop)
    !
    implicit none
    !
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(2)        :: prop12
    real(rp), intent(in ), dimension(0:,0:,0:) :: vof
    real(rp), intent(out), dimension(0:,0:,0:) :: prop
    !
    integer :: i,j,k
    integer :: n1, n2, n3
    real(rp) :: prop12_1, prop12_2
    ! 
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)
    !
    prop12_1 = prop12(1)
    prop12_2 = prop12(2)
    !
    !$acc parallel loop collapse(3)
    do k=1,n3
      do j=1,n2
        do i=1,n1
          prop(i,j,k) = vof(i,j,k)*prop12_1+(1.0_rp-vof(i,j,k))*prop12_2
        enddo
      enddo
    enddo
    !$acc end parallel loop 
    !
    return
  end subroutine update_property
  !
  subroutine clip_vof(n,vof)
    !
    implicit none
    !
    integer , intent(in   ), dimension(3)        :: n
    real(rp), intent(inout), dimension(0:,0:,0:) :: vof
    !
    integer :: i,j,k
    integer :: n1, n2, n3
    !
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)
    !
    !$acc parallel loop collapse(3)
    do k=0,n3+1
      do j=0,n2+1
        do i=0,n1+1
          vof(i,j,k) = min(max(0.0_rp,vof(i,j,k)),1.0_rp) 
        enddo
      enddo
    enddo
    !$acc end parallel loop  
    !
    return
  end subroutine clip_vof
  !
  subroutine cmpt_nor_curv_1o(n1,n2,n3,dli,vof,nor,cur,kappa)
    !
    ! cmpt_nor_curv, 1st order (GPU accelerated)
    !
    implicit none
    !
    integer , intent(in )                         :: n1,n2,n3
    real(rp), intent(in ), dimension(3)           :: dli
    real(rp), intent(in ), dimension(0:,0:,0:   ) :: vof
    real(rp), intent(out), dimension(0:,0:,0:,1:) :: nor,cur
    real(rp), intent(out), dimension(0:,0:,0:   ) :: kappa
    !
    real(rp) :: dndx, dndy, dndz
    real(rp) :: nor1, nor2, nor3
    real(rp) :: norm
    real(rp) :: dli1, dli2, dli3
    integer  :: i,j,k,p
    !@cuf attributes(managed) :: vof, nor, cur, kappa
    !
    dli1 = dli(1)
    dli2 = dli(2)
    dli3 = dli(3)
    !
    !$acc kernels
    do k=1,n3
      do j=1,n2
        do i=1,n1
          !
          !$acc cache(vof(i-1:i+1,j-1:j+1,k-1:k+1))
          nor1 = (4*(vof(i + 1, j, k) - vof(i - 1, j, k)) &
                 + 2*(vof(i + 1, j + 1, k) + vof(i + 1, j, k + 1) + &
                      vof(i + 1, j - 1, k) + vof(i + 1, j, k - 1)) &
                 - 2*(vof(i - 1, j + 1, k) + vof(i - 1, j, k + 1) + &
                      vof(i - 1, j - 1, k) + vof(i - 1, j, k - 1)) &
                 + 1*(vof(i + 1, j + 1, k + 1) + vof(i + 1, j - 1, k + 1) + &
                      vof(i + 1, j + 1, k - 1) + vof(i + 1, j - 1, k - 1)) &
                 - 1*(vof(i - 1, j + 1, k + 1) + vof(i - 1, j - 1, k + 1) + &
                      vof(i - 1, j + 1, k - 1) + vof(i - 1, j - 1, k - 1)))*0.25_rp*0.125_rp*dli1
          nor2 = (4*(vof(i, j + 1, k) - vof(i, j - 1, k)) &
                 + 2*(vof(i + 1, j + 1, k) + vof(i, j + 1, k + 1) + &
                      vof(i - 1, j + 1, k) + vof(i, j + 1, k - 1)) &
                 - 2*(vof(i + 1, j - 1, k) + vof(i, j - 1, k + 1) + &
                      vof(i - 1, j - 1, k) + vof(i, j - 1, k - 1)) &
                 + 1*(vof(i + 1, j + 1, k + 1) + vof(i - 1, j + 1, k - 1) + &
                      vof(i - 1, j + 1, k + 1) + vof(i + 1, j + 1, k - 1)) &
                 - 1*(vof(i + 1, j - 1, k + 1) + vof(i - 1, j - 1, k - 1) + &
                      vof(i - 1, j - 1, k + 1) + vof(i + 1, j - 1, k - 1)))*0.25_rp*0.125_rp*dli2
          nor3 = (4*(vof(i, j, k + 1) - vof(i, j, k - 1)) &
                 + 2*(vof(i + 1, j, k + 1) + vof(i, j + 1, k + 1) + &
                      vof(i - 1, j, k + 1) + vof(i, j - 1, k + 1)) &
                 - 2*(vof(i + 1, j, k - 1) + vof(i, j + 1, k - 1) + &
                      vof(i - 1, j, k - 1) + vof(i, j - 1, k - 1)) &
                 + 1*(vof(i + 1, j + 1, k + 1) + vof(i - 1, j - 1, k + 1) + &
                      vof(i - 1, j + 1, k + 1) + vof(i + 1, j - 1, k + 1)) &
                 - 1*(vof(i + 1, j + 1, k - 1) + vof(i - 1, j - 1, k - 1) + &
                      vof(i - 1, j + 1, k - 1) + vof(i + 1, j - 1, k - 1)))*0.25_rp*0.125_rp*dli3
          !
          norm = sqrt((nor1**2) + (nor2**2) + (nor3**2) + small)
          !
          nor(i,j,k,1) = nor1/norm!*dx
          nor(i,j,k,2) = nor2/norm!*dy
          nor(i,j,k,3) = nor3/norm!*dz
          !if (abs(nor(i,j,k,1)-norchk(i,j,k,1)).gt.0.5d-14) print*, 'ERROR 1=',nor(i,j,k,1)/norchk(i,j,k,1)
          !if (abs(nor(i,j,k,2)-norchk(i,j,k,2)).gt.0.5d-14) print*, 'ERROR 2=',nor(i,j,k,2)/norchk(i,j,k,2)
          !if (abs(nor(i,j,k,3)-norchk(i,j,k,3)).gt.0.5d-14) print*, 'ERROR 3=',nor(i,j,k,3)/norchk(i,j,k,3)
        enddo
      enddo
    enddo
    !
    do k=1,n3
      do j=1,n2
        do i=1,n1
          !
          !$acc cache(nor(i-1:i+1, j      , k      , 1))
          !$acc cache(nor(i      , j-1:j+1, k      , 2))
          !$acc cache(nor(i      , j      , k-1:k+1, 3))
          dndx = (nor(i+1,j,k,1)-nor(i-1,j,k,1))*dli1*0.5_rp
          dndy = (nor(i,j+1,k,2)-nor(i,j-1,k,2))*dli2*0.5_rp
          dndz = (nor(i,j,k+1,3)-nor(i,j,k-1,3))*dli3*0.5_rp
          !
          cur(i,j,k,1) = 0.0_rp
          cur(i,j,k,2) = 0.0_rp
          cur(i,j,k,3) = 0.0_rp
          cur(i,j,k,4) = 0.0_rp
          cur(i,j,k,5) = 0.0_rp
          cur(i,j,k,6) = 0.0_rp
          kappa(i,j,k) = -(dndx + dndy + dndz) ! curvature
          !
         enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine cmpt_nor_curv_1o
  !
#if defined(_FAST_KERNELS_2) 
  subroutine cmpt_nor_curv_2o(n1,n2,n3,dli,vof,nor,cur,kappa)
    !
    ! cmpt_nor_curv, 2nd order (GPU accelerated)
    !
    implicit none
    !
    integer , intent(in )                         :: n1,n2,n3
    real(rp), intent(in ), dimension(3)           :: dli
    real(rp), intent(in ), dimension(0:,0:,0:   ) :: vof
    real(rp), intent(out), dimension(0:,0:,0:,1:) :: nor,cur
    real(rp), intent(out), dimension(0:,0:,0:   ) :: kappa
    !
    real(rp), dimension(8) :: nx,ny,nz,mx,my,mz
    real(rp), dimension(3) :: dl
    real(rp) :: norm
    integer  :: i,j,k,p
    !@cuf attributes(managed) :: vof, nor, cur, kappa
    !
    dl(:) = dli(:)**(-1)
    !
    !$acc parallel loop collapse(3) private(nx,ny,nz,mx,my,mz)    
    do k=1,n3
      do j=1,n2
        do i=1,n1
          !
#if defined(_TWOD) 
          mx(1:8) = 0.0_rp
#else
          !i+1/2 j+1/2 k+1/2
          mx(1) = ((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1)))*dli(1)*0.25_rp
          !i+1/2 j-1/2 k+1/2
          mx(2) = ((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j-1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1)))*dli(1)*0.25_rp
          !i+1/2 j+1/2 k-1/2
          mx(3) = ((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j+1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1)))*dli(1)*0.25_rp
          !i+1/2 j-1/2 k-1/2
          mx(4) = ((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j-1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1)))*dli(1)*0.25_rp
          !i-1/2 j+1/2 k+1/2
          mx(5) = ((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j+1,k+1)))*dli(1)*0.25_rp
          !i-1/2 j-1/2 k+1/2
          mx(6) = ((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j-1,k+1)))*dli(1)*0.25_rp
          !i-1/2 j+1/2 k-1/2
          mx(7) = ((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j+1,k-1)))*dli(1)*0.25_rp
          !i-1/2 j-1/2 k-1/2
          mx(8) = ((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j-1,k-1)))*dli(1)*0.25_rp
#endif
          ! 
          !i+1/2 j+1/2 k+1/2
          my(1) = ((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)))*dli(2)*0.25_rp
          !i+1/2 j-1/2 k+1/2
          my(2) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)) - &
                   (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1)))*dli(2)*0.25_rp
          !i+1/2 j+1/2 k-1/2
          my(3) = ((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)))*dli(2)*0.25_rp
          !i+1/2 j-1/2 k-1/2
          my(4) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)) - &
                   (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli(2)*0.25_rp
          !i-1/2 j+1/2 k+1/2
          my(5) = ((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)))*dli(2)*0.25_rp
          !i-1/2 j-1/2 k+1/2
          my(6) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)) - &
                   (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1)))*dli(2)*0.25_rp
          !i-1/2 j+1/2 k-1/2
          my(7) = ((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)))*dli(2)*0.25_rp
          !i-1/2 j-1/2 k-1/2
          my(8) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)) - &
                   (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli(2)*0.25_rp
          !
          !i+1/2 j+1/2 k+1/2
          mz(1) = ((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  )))*dli(3)*0.25_rp
          !i+1/2 j-1/2 k+1/2
          mz(2) = ((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  )))*dli(3)*0.25_rp
          !i+1/2 j+1/2 k-1/2
          mz(3) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1)))*dli(3)*0.25_rp
          !i+1/2 j-1/2 k-1/2
          mz(4) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli(3)*0.25_rp
          !i-1/2 j+1/2 k+1/2
          mz(5) = ((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  )))*dli(3)*0.25_rp
          !i-1/2 j-1/2 k+1/2
          mz(6) = ((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  )))*dli(3)*0.25_rp
          !i-1/2 j+1/2 k-1/2
          mz(7) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1)))*dli(3)*0.25_rp
          !i-1/2 j-1/2 k-1/2
          mz(8) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli(3)*0.25_rp
          !
          !$acc loop seq
          do p=1,8
            norm  = sqrt(mx(p)**2+my(p)**2+mz(p)**2+small)
            nx(p) = mx(p)/norm
            ny(p) = my(p)/norm
            nz(p) = mz(p)/norm
          enddo
          !$acc end loop
          !
          ! compute the normal vector
          !
          nor(i,j,k,1) = 0._rp
          nor(i,j,k,2) = 0._rp
          nor(i,j,k,3) = 0._rp
          !$acc loop seq
          do p=1,8
            nor(i,j,k,1) = nor(i,j,k,1) + 0.125_rp*mx(p)
            nor(i,j,k,2) = nor(i,j,k,2) + 0.125_rp*my(p)
            nor(i,j,k,3) = nor(i,j,k,3) + 0.125_rp*mz(p)
          enddo
          !$acc end loop
          !
          norm         = sqrt(nor(i,j,k,1)**2+nor(i,j,k,2)**2+nor(i,j,k,3)**2+small)
          nor(i,j,k,1) = nor(i,j,k,1)/norm
          nor(i,j,k,2) = nor(i,j,k,2)/norm
          nor(i,j,k,3) = nor(i,j,k,3)/norm
          ! 
          ! compute the curvature tensor
          !
          cur(i,j,k,1) = ((nx(1)+nx(2)+nx(3)+nx(4))-(nx(5)+nx(6)+nx(7)+nx(8)))*dl(1)*0.25_rp
          cur(i,j,k,2) = ((ny(1)+ny(3)+ny(5)+ny(7))-(ny(2)+ny(4)+ny(6)+ny(8)))*dl(2)*0.25_rp
          cur(i,j,k,3) = ((nz(1)+nz(2)+nz(5)+nz(6))-(nz(3)+nz(4)+nz(7)+nz(8)))*dl(3)*0.25_rp
          cur(i,j,k,4) = ((ny(1)+ny(2)+ny(3)+ny(4))-(ny(5)+ny(6)+ny(7)+ny(8)))*dl(2)*0.25_rp*0.5_rp+&
                         ((nx(1)+nx(3)+nx(5)+nx(7))-(nx(2)+nx(4)+nx(6)+nx(8)))*dl(1)*0.25_rp*0.5_rp
          cur(i,j,k,5) = ((nz(1)+nz(2)+nz(3)+nz(4))-(nz(5)+nz(6)+nz(7)+nz(8)))*dl(3)*0.25_rp*0.5_rp+&
                         ((nx(1)+nx(2)+nx(5)+nx(6))-(nx(3)+nx(4)+nx(7)+nx(8)))*dl(1)*0.25_rp*0.5_rp
          cur(i,j,k,6) = ((nz(1)+nz(3)+nz(5)+nz(7))-(nz(2)+nz(4)+nz(6)+nz(8)))*dl(3)*0.25_rp*0.5_rp+&
                         ((ny(1)+ny(2)+ny(5)+ny(6))-(ny(3)+ny(4)+ny(7)+ny(8)))*dl(2)*0.25_rp*0.5_rp
          !
          kappa(i,j,k) = -(cur(i,j,k,1)*dli(1)**2+cur(i,j,k,2)*dli(2)**2+cur(i,j,k,3)*dli(3)**2) ! curvature
          !$acc loop seq
          do p=1,6
            cur(i,j,k,p) = cur(i,j,k,p)*qu_th
          enddo
          !$acc end loop
          !
        enddo
      enddo
    enddo
    !$acc end parallel loop
    !
    return
  end subroutine cmpt_nor_curv_2o
  !
#else
  !
  subroutine cmpt_nor_curv_2o(n1,n2,n3,dli,vof,nor,cur,kappa)
    !
    ! cmpt_nor_curv, 2nd order (GPU accelerated)
    !
    implicit none
    !
    integer , intent(in )                         :: n1,n2,n3
    real(rp), intent(in ), dimension(3)           :: dli
    real(rp), intent(in ), dimension(0:,0:,0:   ) :: vof
    real(rp), intent(out), dimension(0:,0:,0:,1:) :: nor,cur
    real(rp), intent(out), dimension(0:,0:,0:   ) :: kappa
    !
    real(rp) :: mx1, mx2, mx3, mx4, mx5, mx6, mx7, mx8
    real(rp) :: my1, my2, my3, my4, my5, my6, my7, my8
    real(rp) :: mz1, mz2, mz3, mz4, mz5, mz6, mz7, mz8
    real(rp) :: nx1, nx2, nx3, nx4, nx5, nx6, nx7, nx8
    real(rp) :: ny1, ny2, ny3, ny4, ny5, ny6, ny7, ny8
    real(rp) :: nz1, nz2, nz3, nz4, nz5, nz6, nz7, nz8
    real(rp) :: dli1, dli2, dli3
    real(rp) :: dl1, dl2, dl3
    real(rp) :: norm
    integer  :: i,j,k,p
    !@cuf attributes(managed) :: vof, nor, cur, kappa
    !
    dli1=dli(1)
    dli2=dli(2)
    dli3=dli(3)
    dl1=dli1**(-1)
    dl2=dli2**(-1)
    dl3=dli3**(-1)
    !
    !$acc parallel loop collapse(3)
    do k=1,n3
      do j=1,n2
        do i=1,n1
          !
#if defined(_TWOD) 
          mx1 = 0.0_rp
          mx2 = 0.0_rp
          mx3 = 0.0_rp
          mx4 = 0.0_rp
          mx5 = 0.0_rp
          mx6 = 0.0_rp
          mx7 = 0.0_rp
          mx8 = 0.0_rp
#else
          !i+1/2 j+1/2 k+1/2
          mx1 = ((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1)))*dli1*0.25_rp
          !i+1/2 j-1/2 k+1/2
          mx2 = ((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j-1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1)))*dli1*0.25_rp
          !i+1/2 j+1/2 k-1/2
          mx3 = ((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j+1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1)))*dli1*0.25_rp
          !i+1/2 j-1/2 k-1/2
          mx4 = ((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j-1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1)))*dli1*0.25_rp
          !i-1/2 j+1/2 k+1/2
          mx5 = ((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j+1,k+1)))*dli1*0.25_rp
          !i-1/2 j-1/2 k+1/2
          mx6 = ((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j-1,k+1)))*dli1*0.25_rp
          !i-1/2 j+1/2 k-1/2
          mx7 = ((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j+1,k-1)))*dli1*0.25_rp
          !i-1/2 j-1/2 k-1/2
          mx8 = ((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j-1,k-1)))*dli1*0.25_rp
#endif
          ! 
          !i+1/2 j+1/2 k+1/2
          my1 = ((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)))*dli2*0.25_rp
          !i+1/2 j-1/2 k+1/2
          my2 = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)) - &
                   (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1)))*dli2*0.25_rp
          !i+1/2 j+1/2 k-1/2
          my3 = ((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)))*dli2*0.25_rp
          !i+1/2 j-1/2 k-1/2
          my4 = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)) - &
                   (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli2*0.25_rp
          !i-1/2 j+1/2 k+1/2
          my5 = ((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)))*dli2*0.25_rp
          !i-1/2 j-1/2 k+1/2
          my6 = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)) - &
                   (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1)))*dli2*0.25_rp
          !i-1/2 j+1/2 k-1/2
          my7 = ((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)))*dli2*0.25_rp
          !i-1/2 j-1/2 k-1/2
          my8 = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)) - &
                   (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli2*0.25_rp
          !
          !i+1/2 j+1/2 k+1/2
          mz1 = ((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  )))*dli3*0.25_rp
          !i+1/2 j-1/2 k+1/2
          mz2 = ((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  )))*dli3*0.25_rp
          !i+1/2 j+1/2 k-1/2
          mz3 = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1)))*dli3*0.25_rp
          !i+1/2 j-1/2 k-1/2
          mz4 = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli3*0.25_rp
          !i-1/2 j+1/2 k+1/2
          mz5 = ((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  )))*dli3*0.25_rp
          !i-1/2 j-1/2 k+1/2
          mz6 = ((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  )))*dli3*0.25_rp
          !i-1/2 j+1/2 k-1/2
          mz7 = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1)))*dli3*0.25_rp
          !i-1/2 j-1/2 k-1/2
          mz8 = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli3*0.25_rp
          !
          nx1 = mx1 / ( sqrt(mx1**2 + my1**2+mz1**2+small) )
          ny1 = my1 / ( sqrt(mx1**2 + my1**2+mz1**2+small) )
          nz1 = mz1 / ( sqrt(mx1**2 + my1**2+mz1**2+small) )
          nx2 = mx2 / ( sqrt(mx2**2 + my2**2+mz2**2+small) )
          ny2 = my2 / ( sqrt(mx2**2 + my2**2+mz2**2+small) )
          nz2 = mz2 / ( sqrt(mx2**2 + my2**2+mz2**2+small) )
          nx3 = mx3 / ( sqrt(mx3**2 + my3**2+mz3**2+small) )
          ny3 = my3 / ( sqrt(mx3**2 + my3**2+mz3**2+small) )
          nz3 = mz3 / ( sqrt(mx3**2 + my3**2+mz3**2+small) )
          nx4 = mx4 / ( sqrt(mx4**2 + my4**2+mz4**2+small) )
          ny4 = my4 / ( sqrt(mx4**2 + my4**2+mz4**2+small) )
          nz4 = mz4 / ( sqrt(mx4**2 + my4**2+mz4**2+small) )
          nx5 = mx5 / ( sqrt(mx5**2 + my5**2+mz5**2+small) )
          ny5 = my5 / ( sqrt(mx5**2 + my5**2+mz5**2+small) )
          nz5 = mz5 / ( sqrt(mx5**2 + my5**2+mz5**2+small) )
          nx6 = mx6 / ( sqrt(mx6**2 + my6**2+mz6**2+small) )
          ny6 = my6 / ( sqrt(mx6**2 + my6**2+mz6**2+small) )
          nz6 = mz6 / ( sqrt(mx6**2 + my6**2+mz6**2+small) )
          nx7 = mx7 / ( sqrt(mx7**2 + my7**2+mz7**2+small) )
          ny7 = my7 / ( sqrt(mx7**2 + my7**2+mz7**2+small) )
          nz7 = mz7 / ( sqrt(mx7**2 + my7**2+mz7**2+small) )
          nx8 = mx8 / ( sqrt(mx8**2 + my8**2+mz8**2+small) )
          ny8 = my8 / ( sqrt(mx8**2 + my8**2+mz8**2+small) )
          nz8 = mz8 / ( sqrt(mx8**2 + my8**2+mz8**2+small) )
          !
          ! compute the normal vector
          !
          !nor(i,j,k,1) = nor(i,j,k,1) + 0.125_rp*(mx1 + mx2 + mx3 + mx4 + mx5 + mx6 + mx7 + mx8)
          !nor(i,j,k,2) = nor(i,j,k,2) + 0.125_rp*(my1 + my2 + my3 + my4 + my5 + my6 + my7 + my8)
          !nor(i,j,k,3) = nor(i,j,k,3) + 0.125_rp*(mz1 + mz2 + mz3 + mz4 + mz5 + mz6 + mz7 + mz8)
          nor(i,j,k,1) = 0.125_rp*(mx1 + mx2 + mx3 + mx4 + mx5 + mx6 + mx7 + mx8)
          nor(i,j,k,2) = 0.125_rp*(my1 + my2 + my3 + my4 + my5 + my6 + my7 + my8)
          nor(i,j,k,3) = 0.125_rp*(mz1 + mz2 + mz3 + mz4 + mz5 + mz6 + mz7 + mz8)
          !
          norm         = sqrt(nor(i,j,k,1)**2+nor(i,j,k,2)**2+nor(i,j,k,3)**2+small)
          nor(i,j,k,1) = nor(i,j,k,1)/norm
          nor(i,j,k,2) = nor(i,j,k,2)/norm
          nor(i,j,k,3) = nor(i,j,k,3)/norm
          ! 
          ! compute the curvature tensor
          !
          cur(i,j,k,1) = ((nx1+nx2+nx3+nx4)-(nx5+nx6+nx7+nx8))*dl1*0.25_rp
          cur(i,j,k,2) = ((ny1+ny3+ny5+ny7)-(ny2+ny4+ny6+ny8))*dl2*0.25_rp
          cur(i,j,k,3) = ((nz1+nz2+nz5+nz6)-(nz3+nz4+nz7+nz8))*dl3*0.25_rp
          cur(i,j,k,4) = ((ny1+ny2+ny3+ny4)-(ny5+ny6+ny7+ny8))*dl2*0.25_rp*0.5_rp+&
                         ((nx1+nx3+nx5+nx7)-(nx2+nx4+nx6+nx8))*dl1*0.25_rp*0.5_rp
          cur(i,j,k,5) = ((nz1+nz2+nz3+nz4)-(nz5+nz6+nz7+nz8))*dl3*0.25_rp*0.5_rp+&
                         ((nx1+nx2+nx5+nx6)-(nx3+nx4+nx7+nx8))*dl1*0.25_rp*0.5_rp
          cur(i,j,k,6) = ((nz1+nz3+nz5+nz7)-(nz2+nz4+nz6+nz8))*dl3*0.25_rp*0.5_rp+&
                         ((ny1+ny2+ny5+ny6)-(ny3+ny4+ny7+ny8))*dl2*0.25_rp*0.5_rp
          !
          kappa(i,j,k) = -(cur(i,j,k,1)*dli1**2+cur(i,j,k,2)*dli2**2+cur(i,j,k,3)*dli3**2) ! curvature
          !
          !$acc loop seq
          do p=1,6
            cur(i,j,k,p) = cur(i,j,k,p)*qu_th
          enddo
          !$acc end loop
          !
        enddo
      enddo
    enddo
    !$acc end parallel loop
    !
    return
  end subroutine cmpt_nor_curv_2o
#endif
  !
#if defined(_FAST_KERNELS_3)
  subroutine cmpt_d_thinc(n1,n2,n3,nor,cur,vof,d_thinc)
    !
    implicit none
    !
    integer , intent(in )                         :: n1,n2,n3
    real(rp), intent(in ), dimension(0:,0:,0:,1:) :: nor,cur
    real(rp), intent(in ), dimension(0:,0:,0:)    :: vof
    real(rp), intent(out), dimension(0:,0:,0:)    :: d_thinc
    !
    real(rp), dimension(3,3) :: av,bv,ev
!    real(rp), dimension(3) :: nor_v,cv,dv
!    real(rp), dimension(6) :: cur_v
    real(rp) :: nor_v1, nor_v2, nor_v3
    real(rp) :: cv1, cv2, cv3, cv_curr
    real(rp) :: dv1, dv2, dv3, dv_curr
    real(rp) :: cur_v1, cur_v2, cur_v3, cur_v4, cur_v5, cur_v6 
#if defined(_TWOD)
    real(rp) :: fm,fp,a2,b2,c2
#else
    real(rp) :: fpp,fmm,fpm,fmp,a4,b4,c4,d4,e4
#endif
    real(rp) :: aa, qq, surf
    real(rp) :: dtemp
    real(rp) :: vof_ijk
    integer  :: i,j,k,p,ind2
    !
    real(rp) :: b3,b2,b1,b0,a4i
    real(rp) :: c2,c1,c0
    real(rp) :: z1,check,a,b
    real(rp) :: aa,bb,cc,dd,dd_s
    !
    !@cuf attributes(managed) :: vof, d_thinc, nor, cur
    !
    !$acc parallel loop collapse(3) private(av,bv,ev)
    do k=1,n3
      do j=1,n2
        do i=1,n1
          !
          vof_ijk = vof(i,j,k)
          !
          if((vof_ijk .le.limit).or.(vof_ijk.ge.1.0_rp-limit)) then
            !
            d_thinc(i,j,k) = -1000.0_rp
            !
          else
            !
!            nor_v(:) = nor(i,j,k,:)
!            cur_v(:) = cur(i,j,k,:)
            nor_v1 = nor(i,j,k,1) 
            nor_v2 = nor(i,j,k,2) 
            nor_v3 = nor(i,j,k,3) 
            cur_v1 = cur(i,j,k,1)
            cur_v2 = cur(i,j,k,2)
            cur_v3 = cur(i,j,k,3)
            cur_v4 = cur(i,j,k,4)
            cur_v5 = cur(i,j,k,5)
            cur_v6 = cur(i,j,k,6)
            !
!            cv(:) = 1.0_rp
            !
            ind2 = 0
            if(abs(nor_v1).eq.max(abs(nor_v1),abs(nor_v2),abs(nor_v3))) ind2 = 1
            if(abs(nor_v2).eq.max(abs(nor_v1),abs(nor_v2),abs(nor_v3))) ind2 = 2
            if(abs(nor_v3).eq.max(abs(nor_v1),abs(nor_v2),abs(nor_v3))) ind2 = 3
            !
!            cv(ind2) = 0.0_rp            
            select case(ind2)
              case(1)
                cv1 = 0.0_rp
                cv2 = 1.0_rp
                cv3 = 1.0_rp
              case(2) 
                cv1 = 1.0_rp
                cv2 = 0.0_rp
                cv3 = 1.0_rp
              case(3)
                cv1 = 1.0_rp
                cv2 = 1.0_rp
                cv3 = 0.0_rp
            end select
            !
            ! calculation of the coefficients for the surface function
            !
            ev(:,1) = (/ 0.0_rp, cur_v2, cur_v3/)*0.5_rp
            ev(:,2) = (/ cur_v1, 0.0_rp, cur_v3/)*0.5_rp
            ev(:,3) = (/ cur_v1, cur_v2, 0.0_rp/)*0.5_rp
            av(:,1) = (/ 0.0_rp, 0.0_rp, cur_v6/)
            av(:,2) = (/ 0.0_rp, cur_v5, 0.0_rp/)
            av(:,3) = (/ cur_v4, 0.0_rp, 0.0_rp/)
            bv(:,1) = (/                        0.0_rp, nor_v2-0.5_rp*(cur_v2+cur_v6), nor_v3-0.5_rp*(cur_v3+cur_v6) /)
            bv(:,2) = (/ nor_v1-0.5_rp*(cur_v1+cur_v5),                        0.0_rp, nor_v3-0.5_rp*(cur_v3+cur_v5) /)
            bv(:,3) = (/ nor_v1-0.5_rp*(cur_v1+cur_v4), nor_v2-0.5_rp*(cur_v2+cur_v4),                        0.0_rp /)
            dv1   = ( nor_v1-0.5_rp*cv1*(cur_v1+cv2*cur_v4+cv3*cur_v5) )
            dv2   = ( nor_v2-0.5_rp*cv2*(cur_v2+cv1*cur_v4+cv3*cur_v6) )
            dv3   = ( nor_v3-0.5_rp*cv3*(cur_v3+cv1*cur_v5+cv2*cur_v6) )
            !
            ! build the polynomial equation and find its roots
            !
            aa = 0.0_rp
            qq = 0.0_rp
            !
#if defined(_TWOD)
            !
            ! --> 2D: quadratic equation
            !
            fm = 0.0_rp
            fp = 0.0_rp
            !$acc loop seq
            do p=2,3
              select case(p)
                case(2) 
                  cv_curr = cv2
                  dv_curr = dv2
                case(3)
                  cv_curr = cv3
                  dv_curr = dv3
              end select
              !
              aa = aa + (1.0_rp-cv_curr)*exp(2.0_rp*b_th*dv_curr)
              qq = qq + (1.0_rp-cv_curr)*exp(2.0_rp*b_th*dv_curr*(2.0_rp*vof_ijk-1.0_rp))
              !
              surf = sum(ev(:,p)*xmm(:,p)**2) + sum(bv(:,p)*xmm(:,p)) + &
                     av(1,p)*xmm(1,p)*xmm(2,p)+av(2,p)*xmm(1,p)*xmm(3,p)+av(3,p)*xmm(2,p)*xmm(3,p)
              fm   = fm + (1.0_rp-cv_curr)*exp(2.0_rp*b_th*surf)
              surf = sum(ev(:,p)*xpp(:,p)**2) + sum(bv(:,p)*xpp(:,p)) + &
                     av(1,p)*xpp(1,p)*xpp(2,p)+av(2,p)*xpp(1,p)*xpp(3,p)+av(3,p)*xpp(2,p)*xpp(3,p)
              fp   = fp + (1.0_rp-cv_curr)*exp(2.0_rp*b_th*surf)
              !
            enddo
            !$acc end loop 
            !
            a2 = aa*fm*fp*(aa-qq)
            b2 = aa*(fm+fp)*(1.0_rp-qq)
            c2 = 1.0_rp-aa*qq
            call solve_quad(c2,b2,a2,dtemp)
            !
#else
            !
            ! --> 3D: quartic equation
            !
            fmm = 0.0_rp
            fmp = 0.0_rp
            fpm = 0.0_rp
            fpp = 0.0_rp
            !$acc loop seq
            do p=1,3
              !
              select case(p)
                case(1) 
                  cv_curr = cv1
                  dv_curr = dv1
                case(2) 
                  cv_curr = cv2
                  dv_curr = dv2
                case(3)
                  cv_curr = cv3
                  dv_curr = dv3
              end select
              !
              aa = aa + (1.0_rp-cv_curr)*exp(2.0_rp*b_th*dv_curr)
              qq = qq + (1.0_rp-cv_curr)*exp(4.0_rp*b_th*dv_curr*(2.0_rp*vof_ijk-1.0_rp))
              !
              surf = sum(ev(:,p)*xmm(:,p)**2) + sum(bv(:,p)*xmm(:,p)) + &
                     av(1,p)*xmm(1,p)*xmm(2,p)+av(2,p)*xmm(1,p)*xmm(3,p)+av(3,p)*xmm(2,p)*xmm(3,p)
              fmm  = fmm + (1.0_rp-cv_curr)*exp(2.0_rp*b_th*surf)
              surf = sum(ev(:,p)*xmp(:,p)**2) + sum(bv(:,p)*xmp(:,p)) + &
                     av(1,p)*xmp(1,p)*xmp(2,p)+av(2,p)*xmp(1,p)*xmp(3,p)+av(3,p)*xmp(2,p)*xmp(3,p)
              fmp  = fmp + (1.0_rp-cv_curr)*exp(2.0_rp*b_th*surf)
              surf = sum(ev(:,p)*xpm(:,p)**2) + sum(bv(:,p)*xpm(:,p)) + &
                     av(1,p)*xpm(1,p)*xpm(2,p)+av(2,p)*xpm(1,p)*xpm(3,p)+av(3,p)*xpm(2,p)*xpm(3,p)
              fpm  = fpm + (1.0_rp-cv_curr)*exp(2.0_rp*b_th*surf)
              surf = sum(ev(:,p)*xpp(:,p)**2) + sum(bv(:,p)*xpp(:,p)) + &
                     av(1,p)*xpp(1,p)*xpp(2,p)+av(2,p)*xpp(1,p)*xpp(3,p)+av(3,p)*xpp(2,p)*xpp(3,p)
              fpp  = fpp + (1.0_rp-cv_curr)*exp(2.0_rp*b_th*surf)
              !
            enddo
            !$acc end loop 
            !
            a4 = aa**2*(fmm*fmp*fpm*fpp)*(aa**2-qq)
            b4 = aa**2*(fmm*fpm*fpp+fmp*fpm*fpp+fmm*fmp*fpm+fmm*fmp*fpp)*(aa-qq)
            c4 = aa**2*(fpm*fpp+fmm*fpm+fmp*fpm+fmm*fpp+fmm*fmp+fmp*fpp)*(1.0_rp-qq)
            d4 = aa*(fpm+fpp+fmm+fmp)*(1.0_rp-aa*qq)
            e4 = 1.0_rp-aa**2*qq
            !
            call solve_quar_paper(e4,d4,c4,b4,a4,dtemp)
#if 0 
            ! *** DISABLED CODE ***
            ! [NOTE] Enable this piece of code generates a crash at run-time for 512^3 test-case, not sure why ...
             associate( a0 => e4 ,a1 => d4 ,a2 => c4 ,a3 => b4, x => dtemp )
              ! 
              ! First change the coefficients in the form of Eq. (B.2)
              !
              a4i = 1._rp/a4
              b3  = a3*a4i
              b2  = a2*a4i
              b1  = a1*a4i
              b0  = a0*a4i
              !
              ! Calculate the coefficients of Eq. (B.5)
              !  note: in c1 it is not b2 (as in the paper), but b3!
              !
              c2  = -b2
              c1  = b1*b3-4._rp*b0
              c0  = b0*(4._rp*b2-b3**2)-b1**2
              !
              ! Calculate z1, Eq. (B.7)
              !
              a   = -c2**2/9._rp+c1/3._rp
              b   = 2._rp*c2**3/27._rp-c1*c2/3._rp+c0
              check = b**2+4._rp*a**3
              !
              if(check.ge.0._rp) then
                z1 = sign(1._rp,(-b+sqrt(check)))*abs(0.5_rp*(-b+sqrt(check)))**(1._rp/3._rp) + &
                    sign(1._rp,(-b-sqrt(check)))*abs(0.5_rp*(-b-sqrt(check)))**(1._rp/3._rp) - &
                    c2/3._rp
              else
                z1 = 2._rp*sqrt(-a)*cos(atan(sqrt(-check)/(-b))/3._rp)-c2/3._rp
              endif
              !
              ! Find new coefficients, Eq. (B.9)
              !
              aa   = 0.5_rp*b3
              bb   = 0.5_rp*z1
              dd_s = bb**2-b0 ! as said in the paper, this should be always positive but still we put a check
              if(dd_s.le.0._rp) then
                dd = limit    ! to avoid singularity in case of imaginary solution (which we filter out)
              else
                dd = sqrt(dd_s)
              endif
              cc = (-0.5_rp*b1+aa*bb)/dd
              !
              ! Finally, calculate solution from Eq. (B.11)
              ! Be aware of a small error in the paper, i.e. 
              ! the (+) sign inside the square root should be (-)
              !
              x = 0.5_rp*(-(aa-cc) + sqrt((aa-cc)**2-4._rp*(bb-dd)))
              !
            end associate  
#endif    
            !
#endif
            !
            d_thinc(i,j,k) = 0.5_rp*b_thi*log(dtemp)
            !
          endif
          !
        enddo
      enddo
    enddo
    !$acc end parallel loop 
    !
    return
  end subroutine cmpt_d_thinc
  !
#else
  !
  subroutine cmpt_d_thinc(n1,n2,n3,nor,cur,vof,d_thinc)
    !
    implicit none
    !
    integer , intent(in )                         :: n1,n2,n3
    real(rp), intent(in ), dimension(0:,0:,0:,1:) :: nor,cur
    real(rp), intent(in ), dimension(0:,0:,0:)    :: vof
    real(rp), intent(out), dimension(0:,0:,0:)    :: d_thinc
    !
    real(rp), dimension(3,3) :: av,bv,ev
    real(rp), dimension(3) :: nor_v,cv,dv
    real(rp), dimension(6) :: cur_v
#if defined(_TWOD)
    real(rp) :: fm,fp,a2,b2,c2
#else
    real(rp) :: fpp,fmm,fpm,fmp,a4,b4,c4,d4,e4
#endif
    real(rp) :: aa,qq,surf
    real(rp) :: dtemp
    integer  :: i,j,k,p,ind2
    !@cuf attributes(managed) :: vof, d_thinc, nor, cur
    !
    !$acc parallel loop collapse(3) private(nor_v,cv,dv,av,bv,ev,cur_v)
    do k=1,n3
      do j=1,n2
        do i=1,n1
          !
          if((vof(i,j,k).le.limit).or.(vof(i,j,k).ge.1.0_rp-limit)) then
            !
            d_thinc(i,j,k) = -1000.0_rp
            !
          else
            !
            nor_v(:) = nor(i,j,k,:)
            cur_v(:) = cur(i,j,k,:)
            cv(:) = 1.0_rp
            !
            ind2 = 0
            if(abs(nor(i,j,k,1)).eq.max(abs(nor(i,j,k,1)),abs(nor(i,j,k,2)),abs(nor(i,j,k,3)))) ind2 = 1
            if(abs(nor(i,j,k,2)).eq.max(abs(nor(i,j,k,1)),abs(nor(i,j,k,2)),abs(nor(i,j,k,3)))) ind2 = 2
            if(abs(nor(i,j,k,3)).eq.max(abs(nor(i,j,k,1)),abs(nor(i,j,k,2)),abs(nor(i,j,k,3)))) ind2 = 3
            !
            cv(ind2) = 0.0_rp            
            !
            ! calculation of the coefficients for the surface function
            !
            ev(:,1) = (/0.0_rp,cur_v(2),cur_v(3)/)*0.5_rp
            ev(:,2) = (/cur_v(1),0.0_rp,cur_v(3)/)*0.5_rp
            ev(:,3) = (/cur_v(1),cur_v(2),0.0_rp/)*0.5_rp
            av(:,1) = (/0.0_rp,0.0_rp,cur_v(6)/)
            av(:,2) = (/0.0_rp,cur_v(5),0.0_rp/)
            av(:,3) = (/cur_v(4),0.0_rp,0.0_rp/)
            bv(:,1) = (/0.0_rp,nor_v(2)-0.5_rp*(cur_v(2)+cur_v(6)),nor_v(3)-0.5_rp*(cur_v(3)+cur_v(6))/)
            bv(:,2) = (/nor_v(1)-0.5_rp*(cur_v(1)+cur_v(5)),0.0_rp,nor_v(3)-0.5_rp*(cur_v(3)+cur_v(5))/)
            bv(:,3) = (/nor_v(1)-0.5_rp*(cur_v(1)+cur_v(4)),nor_v(2)-0.5_rp*(cur_v(2)+cur_v(4)),0.0_rp/)
            dv(1)   = (nor_v(1)-0.5_rp*cv(1)*(cur_v(1)+cv(2)*cur_v(4)+cv(3)*cur_v(5)))
            dv(2)   = (nor_v(2)-0.5_rp*cv(2)*(cur_v(2)+cv(1)*cur_v(4)+cv(3)*cur_v(6)))
            dv(3)   = (nor_v(3)-0.5_rp*cv(3)*(cur_v(3)+cv(1)*cur_v(5)+cv(2)*cur_v(6)))
            !
            ! build the polynomial equation and find its roots
            !
            aa = 0.0_rp
            qq = 0.0_rp
            !
#if defined(_TWOD)
            !
            ! --> 2D: quadratic equation
            !
            fm = 0.0_rp
            fp = 0.0_rp
            !$acc loop seq
            do p=2,3
              !
              aa = aa + (1.0_rp-cv(p))*exp(2.0_rp*b_th*dv(p))
              qq = qq + (1.0_rp-cv(p))*exp(2.0_rp*b_th*dv(p)*(2.0_rp*vof(i,j,k)-1.0_rp))
              !
              surf = sum(ev(:,p)*xmm(:,p)**2) + sum(bv(:,p)*xmm(:,p)) + &
                     av(1,p)*xmm(1,p)*xmm(2,p)+av(2,p)*xmm(1,p)*xmm(3,p)+av(3,p)*xmm(2,p)*xmm(3,p)
              fm   = fm + (1.0_rp-cv(p))*exp(2.0_rp*b_th*surf)
              surf = sum(ev(:,p)*xpp(:,p)**2) + sum(bv(:,p)*xpp(:,p)) + &
                     av(1,p)*xpp(1,p)*xpp(2,p)+av(2,p)*xpp(1,p)*xpp(3,p)+av(3,p)*xpp(2,p)*xpp(3,p)
              fp   = fp + (1.0_rp-cv(p))*exp(2.0_rp*b_th*surf)
              !
            enddo
            !$acc end loop 
            !
            a2 = aa*fm*fp*(aa-qq)
            b2 = aa*(fm+fp)*(1.0_rp-qq)
            c2 = 1.0_rp-aa*qq
            call solve_quad(c2,b2,a2,dtemp)
            !
#else
            !
            ! --> 3D: quartic equation
            !
            fmm = 0.0_rp
            fmp = 0.0_rp
            fpm = 0.0_rp
            fpp = 0.0_rp
            do p=1,3
              !
              aa = aa + (1.0_rp-cv(p))*exp(2.0_rp*b_th*dv(p))
              qq = qq + (1.0_rp-cv(p))*exp(4.0_rp*b_th*dv(p)*(2.0_rp*vof(i,j,k)-1.0_rp))
              !
              surf = sum(ev(:,p)*xmm(:,p)**2) + sum(bv(:,p)*xmm(:,p)) + &
                     av(1,p)*xmm(1,p)*xmm(2,p)+av(2,p)*xmm(1,p)*xmm(3,p)+av(3,p)*xmm(2,p)*xmm(3,p)
              fmm  = fmm + (1.0_rp-cv(p))*exp(2.0_rp*b_th*surf)
              surf = sum(ev(:,p)*xmp(:,p)**2) + sum(bv(:,p)*xmp(:,p)) + &
                     av(1,p)*xmp(1,p)*xmp(2,p)+av(2,p)*xmp(1,p)*xmp(3,p)+av(3,p)*xmp(2,p)*xmp(3,p)
              fmp  = fmp + (1.0_rp-cv(p))*exp(2.0_rp*b_th*surf)
              surf = sum(ev(:,p)*xpm(:,p)**2) + sum(bv(:,p)*xpm(:,p)) + &
                     av(1,p)*xpm(1,p)*xpm(2,p)+av(2,p)*xpm(1,p)*xpm(3,p)+av(3,p)*xpm(2,p)*xpm(3,p)
              fpm  = fpm + (1.0_rp-cv(p))*exp(2.0_rp*b_th*surf)
              surf = sum(ev(:,p)*xpp(:,p)**2) + sum(bv(:,p)*xpp(:,p)) + &
                     av(1,p)*xpp(1,p)*xpp(2,p)+av(2,p)*xpp(1,p)*xpp(3,p)+av(3,p)*xpp(2,p)*xpp(3,p)
              fpp  = fpp + (1.0_rp-cv(p))*exp(2.0_rp*b_th*surf)
              !
            enddo
            !
            a4 = aa**2*(fmm*fmp*fpm*fpp)*(aa**2-qq)
            b4 = aa**2*(fmm*fpm*fpp+fmp*fpm*fpp+fmm*fmp*fpm+fmm*fmp*fpp)*(aa-qq)
            c4 = aa**2*(fpm*fpp+fmm*fpm+fmp*fpm+fmm*fpp+fmm*fmp+fmp*fpp)*(1.0_rp-qq)
            d4 = aa*(fpm+fpp+fmm+fmp)*(1.0_rp-aa*qq)
            e4 = 1.0_rp-aa**2*qq
            call solve_quar_paper(e4,d4,c4,b4,a4,dtemp)
            !
#endif
            !
            d_thinc(i,j,k) = 0.5_rp*b_thi*log(dtemp)
            !
          endif
          !
        enddo
      enddo
    enddo
    !$acc end parallel loop 
    !
    return
  end subroutine cmpt_d_thinc
#endif
  !
#if defined(_FAST_KERNELS_4)
  subroutine cmpt_vof_flux(n1, n2, n3,dli,dt,nh_u,vof,nor,cur,d_thinc,dir,vel,flux)
    !
    implicit none
    !
    integer , intent(in )                                     :: n1,n2,n3
    real(rp), intent(in )                                     :: dli
    real(rp), intent(in )                                     :: dt
    integer , intent(in )                                     :: nh_u
    real(rp), intent(in ), dimension(0:,0:,0:   )             :: vof
    real(rp), intent(in ), dimension(0:,0:,0:,1:)             :: nor,cur
    real(rp), intent(in ), dimension(0:,0:,0:   )             :: d_thinc
    integer , intent(in )                                     :: dir
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: vel
    real(rp), intent(out), dimension(     0:,     0:,     0:) :: flux
    !
!    real(rp), dimension(3)   :: nor_v,cv
    real(rp) :: nor_v1, nor_v2, nor_v3
    real(rp) :: cv1, cv2, cv3 
!    real(rp), dimension(6)   :: cur_v
    real(rp)  :: cur_v1, cur_v2, cur_v3, cur_v4, cur_v5, cur_v6 
    real(rp), dimension(6,3) :: f2_l
    real(rp), dimension(8,3) :: xv
!    real(rp), dimension(6)   :: crd
    real(rp) :: crd1, crd2, crd3, crd4, crd5, crd6
!    integer , dimension(3)   :: f1_l
    integer  :: f1_l1, f1_l2, f1_l3
    real(rp) :: xa,xb,ya,yb,za,zb,cf,vel_sign,dl
    real(rp) :: a,b,rm1,rp1,rm2,rp2
    real(rp) :: cxxa,cyya,czza,cxya,cyza,cxza,a100,a010,a001,sum_int_func
    real(rp) :: surf_a,surf_b
    integer  :: ind2,ii,jj,kk,i,j,k,p,q
    !
    real(rp), parameter :: gm = -1.0_rp/sqrt(3.0_rp), &
                           gp = +1.0_rp/sqrt(3.0_rp)
    !@cuf attributes(managed) :: vof, nor, cur, d_thinc, vel, flux
    !
    dl = 1.0_rp/dli
    !
    select case(dir)
      case(1)
!       crd(:) = (/0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp,1.0_rp/) ! along x
        crd1 = 0.0_rp
        crd2 = 0.0_rp
        crd3 = 0.0_rp
        crd4 = 1.0_rp
        crd5 = 0.0_rp
        crd6 = 1.0_rp
      case(2) 
!       crd(:) = (/0.0_rp,1.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp/) ! along y
        crd1 = 0.0_rp
        crd2 = 1.0_rp
        crd3 = 0.0_rp
        crd4 = 0.0_rp
        crd5 = 0.0_rp
        crd6 = 1.0_rp
      case(3)
!       crd(:) = (/0.0_rp,1.0_rp,0.0_rp,1.0_rp,0.0_rp,0.0_rp/) ! along z
        crd1 = 0.0_rp
        crd2 = 1.0_rp
        crd3 = 0.0_rp
        crd4 = 1.0_rp
        crd5 = 0.0_rp
        crd6 = 0.0_rp
    end select
    !
    !$acc parallel loop collapse(3) private(f2_l, xv) create(f2_l, xv)
    do k=0,n3
      do j=0,n2
        do i=0,n1
          !
          f2_l(:,:) = 0.0_rp
          !
          ! decide the upwind path
          !
          vel_sign = sign(1.0_rp,vel(i,j,k))
          !
          ! compute the indexes, ii,jj,kk
          !  note: they have to computed for each (i,j,k) since
          !        they depend on the local velocity and direction
          !
          ! f1_l(dir) = nint(0.5_rp*(-vel_sign+1.0_rp)) ! i.e., f1_l(dir) = 0 (vel>=0) or 1 (vel<0)
          select case(dir)
            case(1)
              f1_l1   = nint(0.5_rp*(-vel_sign+1.0_rp)) ! i.e., f1_l(dir) = 0 (vel>=0) or 1 (vel<0)
              f1_l2   = 0
              f1_l3   = 0
            case(2) 
              f1_l1   = 0
              f1_l2   = nint(0.5_rp*(-vel_sign+1.0_rp)) ! i.e., f1_l(dir) = 0 (vel>=0) or 1 (vel<0)
              f1_l3   = 0
            case(3)
              f1_l1   = 0
              f1_l2   = 0
              f1_l3   = nint(0.5_rp*(-vel_sign+1.0_rp)) ! i.e., f1_l(dir) = 0 (vel>=0) or 1 (vel<0)
          end select
          !
          ii = i + f1_l1
          jj = j + f1_l2
          kk = k + f1_l3
          !
          ! compute the integration extrema, xa,xb,ya,yb,za,zb
          !  note: they have to computed for each (i,j,k) since
          !        they depend on the local velocity and direction
          !
          cf = 0.5_rp*(vel_sign+1.0_rp)
          !f2_l(2*dir-1+f1_l(dir),dir) = cf-dt*dli*vel(i,j,k) ! we touch only the integration direction
          !f2_l(2*dir+0-f1_l(dir),dir) = cf
          select case(dir)
            case(1)
              f2_l(2*dir-1+f1_l1,dir) = cf-dt*dli*vel(i,j,k) ! we touch only the integration direction
              f2_l(2*dir+0-f1_l1,dir) = cf
            case(2) 
              f2_l(2*dir-1+f1_l2,dir) = cf-dt*dli*vel(i,j,k) ! we touch only the integration direction
              f2_l(2*dir+0-f1_l2,dir) = cf
            case(3)
              f2_l(2*dir-1+f1_l3,dir) = cf-dt*dli*vel(i,j,k) ! we touch only the integration direction
              f2_l(2*dir+0-f1_l3,dir) = cf
          end select
          !
          xa = crd1 + f2_l(1,dir)
          xb = crd2 + f2_l(2,dir)
          ya = crd3 + f2_l(3,dir)
          yb = crd4 + f2_l(4,dir)
          za = crd5 + f2_l(5,dir)
          zb = crd6 + f2_l(6,dir)
          !
          ! decide dir2
          !
          if((vof(ii,jj,kk).le.limit).or.(vof(ii,jj,kk).ge.1.0_rp-limit)) then
            !
            flux(i,j,k) = vof(ii,jj,kk)*(xb-xa)*(yb-ya)*(zb-za)
            !
          else
            !
!            nor_v(:) = nor(ii,jj,kk,:) 
            nor_v1 = nor(ii,jj,kk,1) 
            nor_v2 = nor(ii,jj,kk,2) 
            nor_v3 = nor(ii,jj,kk,3) 
!            cur_v(:) = cur(ii,jj,kk,:)
            cur_v1 = cur(ii,jj,kk,1)
            cur_v2 = cur(ii,jj,kk,2)
            cur_v3 = cur(ii,jj,kk,3)
            cur_v4 = cur(ii,jj,kk,4)
            cur_v5 = cur(ii,jj,kk,5)
            cur_v6 = cur(ii,jj,kk,6)
            !
            ! ind2 = maxloc(abs(nor_v(:)),1)
            if(    abs(nor(ii,jj,kk,1)).eq.max(abs(nor(ii,jj,kk,1)),abs(nor(ii,jj,kk,2)),abs(nor(ii,jj,kk,3)))) then
              ind2 = 1
            elseif(abs(nor(ii,jj,kk,2)).eq.max(abs(nor(ii,jj,kk,1)),abs(nor(ii,jj,kk,2)),abs(nor(ii,jj,kk,3)))) then
              ind2 = 2
            elseif(abs(nor(ii,jj,kk,3)).eq.max(abs(nor(ii,jj,kk,1)),abs(nor(ii,jj,kk,2)),abs(nor(ii,jj,kk,3)))) then
              ind2 = 3
            end if
            !
            if (ind2.eq.1) then
              !
              a    = xa
              b    = xb
              rm1  = 0.5_rp*((yb-ya)*gm+(ya+yb))
              rp1  = 0.5_rp*((yb-ya)*gp+(ya+yb))
              rm2  = 0.5_rp*((zb-za)*gm+(za+zb))
              rp2  = 0.5_rp*((zb-za)*gp+(za+zb))
              flux(i,j,k) = (yb-ya)*(zb-za)
              !
              xv(1, 1) = b; xv(1, 2) = rm1; xv(1, 3) = rm2
              xv(2, 1) = a; xv(2, 2) = rm1; xv(2, 3) = rm2
              xv(3, 1) = b; xv(3, 2) = rm1; xv(3, 3) = rp2
              xv(4, 1) = a; xv(4, 2) = rm1; xv(4, 3) = rp2
              xv(5, 1) = b; xv(5, 2) = rp1; xv(5, 3) = rm2
              xv(6, 1) = a; xv(6, 2) = rp1; xv(6, 3) = rm2
              xv(7, 1) = b; xv(7, 2) = rp1; xv(7, 3) = rp2
              xv(8, 1) = a; xv(8, 2) = rp1; xv(8, 3) = rp2
              !
              !xv(1,:) = (/b,rm1,rm2/) 
              !xv(2,:) = (/a,rm1,rm2/)
              !xv(3,:) = (/b,rm1,rp2/)
              !xv(4,:) = (/a,rm1,rp2/)
              !xv(5,:) = (/b,rp1,rm2/)
              !xv(6,:) = (/a,rp1,rm2/)
              !xv(7,:) = (/b,rp1,rp2/)
              !xv(8,:) = (/a,rp1,rp2/)
              !
            elseif (ind2.eq.2) then
              !
              a    = ya
              b    = yb
              rm1  = 0.5_rp*((xb-xa)*gm+(xa+xb))
              rp1  = 0.5_rp*((xb-xa)*gp+(xa+xb))
              rm2  = 0.5_rp*((zb-za)*gm+(za+zb))
              rp2  = 0.5_rp*((zb-za)*gp+(za+zb))
              flux(i,j,k) = (xb-xa)*(zb-za)
              !
              xv(1, 1) = rm1; xv(1, 2) = b; xv(1, 3) = rm2
              xv(2, 1) = rm1; xv(2, 2) = a; xv(2, 3) = rm2
              xv(3, 1) = rm1; xv(3, 2) = b; xv(3, 3) = rp2
              xv(4, 1) = rm1; xv(4, 2) = a; xv(4, 3) = rp2
              xv(5, 1) = rp1; xv(5, 2) = b; xv(5, 3) = rm2
              xv(6, 1) = rp1; xv(6, 2) = a; xv(6, 3) = rm2
              xv(7, 1) = rp1; xv(7, 2) = b; xv(7, 3) = rp2
              xv(8, 1) = rp1; xv(8, 2) = a; xv(8, 3) = rp2
              !
              !xv(1,:) = (/rm1,b,rm2/) 
              !xv(2,:) = (/rm1,a,rm2/)
              !xv(3,:) = (/rm1,b,rp2/)
              !xv(4,:) = (/rm1,a,rp2/)
              !xv(5,:) = (/rp1,b,rm2/)
              !xv(6,:) = (/rp1,a,rm2/)
              !xv(7,:) = (/rp1,b,rp2/)
              !xv(8,:) = (/rp1,a,rp2/)
              !
            elseif (ind2.eq.3) then
              !
              a    = za
              b    = zb
              rm1  = 0.5_rp*((xb-xa)*gm+(xa+xb))
              rp1  = 0.5_rp*((xb-xa)*gp+(xa+xb))
              rm2  = 0.5_rp*((yb-ya)*gm+(ya+yb))
              rp2  = 0.5_rp*((yb-ya)*gp+(ya+yb))
              flux(i,j,k) = (xb-xa)*(yb-ya)
              !
              xv(1, 1) = rm1; xv(1, 2) = rm2; xv(1, 3) = b
              xv(2, 1) = rm1; xv(2, 2) = rm2; xv(2, 3) = a
              xv(3, 1) = rm1; xv(3, 2) = rp2; xv(3, 3) = b
              xv(4, 1) = rm1; xv(4, 2) = rp2; xv(4, 3) = a
              xv(5, 1) = rp1; xv(5, 2) = rm2; xv(5, 3) = b
              xv(6, 1) = rp1; xv(6, 2) = rm2; xv(6, 3) = a
              xv(7, 1) = rp1; xv(7, 2) = rp2; xv(7, 3) = b
              xv(8, 1) = rp1; xv(8, 2) = rp2; xv(8, 3) = a
              !
              !xv(1,:) = (/rm1,rm2,b/) 
              !xv(2,:) = (/rm1,rm2,a/)
              !xv(3,:) = (/rm1,rp2,b/)
              !xv(4,:) = (/rm1,rp2,a/)
              !xv(5,:) = (/rp1,rm2,b/)
              !xv(6,:) = (/rp1,rm2,a/)
              !xv(7,:) = (/rp1,rp2,b/)
              !xv(8,:) = (/rp1,rp2,a/)
              !
            endif
            !
!            cv(:)    = 1.0_rp
!            cv(ind2) = 0.0_rp
            select case(ind2)
              case(1)
                cv1 = 0.0_rp
                cv2 = 1.0_rp
                cv3 = 1.0_rp
              case(2) 
                cv1 = 1.0_rp
                cv2 = 0.0_rp
                cv3 = 1.0_rp
              case(3)
                cv1 = 1.0_rp
                cv2 = 1.0_rp
                cv3 = 0.0_rp
            end select
            !
            cxxa = cv1*0.5_rp*cv1*cur_v1
            cyya = cv2*0.5_rp*cv2*cur_v2
            czza = cv3*0.5_rp*cv3*cur_v3
            cxya = cv1*cv2*cur_v4
            cxza = cv1*cv3*cur_v5
            cyza = cv2*cv3*cur_v6
            a100 = (nor_v1-0.5_rp*cv1*(cur_v1+cv2*cur_v4+cv3*cur_v5))
            a010 = (nor_v2-0.5_rp*cv2*(cur_v2+cv1*cur_v4+cv3*cur_v6))
            a001 = (nor_v3-0.5_rp*cv3*(cur_v3+cv1*cur_v5+cv2*cur_v6))
            !
            sum_int_func = 0.0_rp
            !$acc loop seq
            do p=1,7,2
              !
              q = p+1
              surf_a = cxxa*xv(p,1)*xv(p,1) + cyya*xv(p,2)*xv(p,2) + czza*xv(p,3)*xv(p,3) + &
                       cxya*xv(p,1)*xv(p,2) + cyza*xv(p,2)*xv(p,3) + cxza*xv(p,3)*xv(p,1) + &
                       a100*xv(p,1)+a010*xv(p,2)+a001*xv(p,3)
              surf_b = cxxa*xv(q,1)*xv(q,1) + cyya*xv(q,2)*xv(q,2) + czza*xv(q,3)*xv(q,3) + &
                       cxya*xv(q,1)*xv(q,2) + cyza*xv(q,2)*xv(q,3) + cxza*xv(q,3)*xv(q,1) + &
                       a100*xv(q,1)+a010*xv(q,2)+a001*xv(q,3)
              !

              select case(ind2)
                case(1)
                   sum_int_func = sum_int_func + &
                             (b-a+b_thi/(nor_v1+limit)*log( &
                             cosh(b_th*(surf_a+d_thinc(ii,jj,kk))) / &
                             cosh(b_th*(surf_b+d_thinc(ii,jj,kk)))))
                case(2) 
                  sum_int_func = sum_int_func + &
                             (b-a+b_thi/(nor_v2+limit)*log( &
                             cosh(b_th*(surf_a+d_thinc(ii,jj,kk))) / &
                             cosh(b_th*(surf_b+d_thinc(ii,jj,kk)))))
                case(3)
                  sum_int_func = sum_int_func + &
                             (b-a+b_thi/(nor_v3+limit)*log( &
                             cosh(b_th*(surf_a+d_thinc(ii,jj,kk))) / &
                             cosh(b_th*(surf_b+d_thinc(ii,jj,kk)))))
              end select
              !
            enddo
            !$acc end loop
            !
            flux(i,j,k) = 0.125_rp*flux(i,j,k)*sum_int_func
            !
          endif
          !
          flux(i,j,k) = flux(i,j,k)*vel_sign*dl
          !
        enddo
      enddo
    enddo
    !$acc end parallel loop
    !
    return
  end subroutine cmpt_vof_flux
  !
#else
  !
  subroutine cmpt_vof_flux(n1, n2, n3,dli,dt,nh_u,vof,nor,cur,d_thinc,dir,vel,flux)
    !
    implicit none
    !
    integer , intent(in )                                     :: n1,n2,n3
    real(rp), intent(in )                                     :: dli
    real(rp), intent(in )                                     :: dt
    integer , intent(in )                                     :: nh_u
    real(rp), intent(in ), dimension(0:,0:,0:   )             :: vof
    real(rp), intent(in ), dimension(0:,0:,0:,1:)             :: nor,cur
    real(rp), intent(in ), dimension(0:,0:,0:   )             :: d_thinc
    integer , intent(in )                                     :: dir
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: vel
    real(rp), intent(out), dimension(     0:,     0:,     0:) :: flux
    !
    real(rp), dimension(3)   :: nor_v,cv
    real(rp), dimension(6)   :: cur_v
    real(rp), dimension(6,3) :: f2_l
    real(rp), dimension(8,3) :: xv
    real(rp), dimension(6)   :: crd
    integer , dimension(3)   :: f1_l
    real(rp) :: xa,xb,ya,yb,za,zb,cf,vel_sign,dl
    real(rp) :: a,b,rm1,rp1,rm2,rp2
    real(rp) :: cxxa,cyya,czza,cxya,cyza,cxza,a100,a010,a001,sum_int_func
    real(rp) :: surf_a,surf_b
    integer  :: ind2,ii,jj,kk,i,j,k,p,q
    !
    real(rp), parameter :: gm = -1.0_rp/sqrt(3.0_rp), &
                           gp = +1.0_rp/sqrt(3.0_rp)
    !@cuf attributes(managed) :: vof, nor, cur, d_thinc, vel, flux
    !
    dl = 1.0_rp/dli
    !
#if 0
    ! *DISABLED CODE*
    select case(dir)
    case(1)
     crd(:) = (/0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp,1.0_rp/) ! along x
    case(2) 
     crd(:) = (/0.0_rp,1.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp/) ! along y
    case(3)
     crd(:) = (/0.0_rp,1.0_rp,0.0_rp,1.0_rp,0.0_rp,0.0_rp/) ! along z
    end select
#endif
    !
    ! initialize two auxiliary arrays
    !
#if 0
    ! *DISABLED CODE*
    f1_l(:)   = 0
    f2_l(:,:) = 0.0_rp
#endif
    !
    !$acc parallel loop collapse(3) private(nor_v, cv,cur_v, f2_l, xv, crd, f1_l)
    do k=0,n3
      do j=0,n2
        do i=0,n1
          !
          select case(dir)
          case(1)
           crd(:) = (/0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp,1.0_rp/) ! along x
          case(2) 
           crd(:) = (/0.0_rp,1.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp/) ! along y
          case(3)
           crd(:) = (/0.0_rp,1.0_rp,0.0_rp,1.0_rp,0.0_rp,0.0_rp/) ! along z
          end select
          !
          f1_l(:)   = 0
          f2_l(:,:) = 0.0_rp
          !
          ! decide the upwind path
          !
          vel_sign = sign(1.0_rp,vel(i,j,k))
          !
          ! compute the indexes, ii,jj,kk
          !  note: they have to computed for each (i,j,k) since
          !        they depend on the local velocity and direction
          !
          f1_l(dir) = nint(0.5_rp*(-vel_sign+1.0_rp)) ! i.e., f1_l(dir) = 0 (vel>=0) or 1 (vel<0)
          ii = i+f1_l(1)
          jj = j+f1_l(2)
          kk = k+f1_l(3)
          !
          ! compute the integration extrema, xa,xb,ya,yb,za,zb
          !  note: they have to computed for each (i,j,k) since
          !        they depend on the local velocity and direction
          !
          cf = 0.5_rp*(vel_sign+1.0_rp) ! i.e., cf = 1.0 (vel>=0) or 0.0 (vel<0)
          f2_l(2*dir-1+f1_l(dir),dir) = cf-dt*dli*vel(i,j,k) ! we touch only the integration direction
          f2_l(2*dir+0-f1_l(dir),dir) = cf
          !
          xa = crd(1)+f2_l(1,dir)
          xb = crd(2)+f2_l(2,dir)
          ya = crd(3)+f2_l(3,dir)
          yb = crd(4)+f2_l(4,dir)
          za = crd(5)+f2_l(5,dir)
          zb = crd(6)+f2_l(6,dir)
          !
          ! decide dir2
          !
          if((vof(ii,jj,kk).le.limit).or.(vof(ii,jj,kk).ge.1.0_rp-limit)) then
            !
            flux(i,j,k) = vof(ii,jj,kk)*(xb-xa)*(yb-ya)*(zb-za)
            !
          else
            !
            nor_v(:) = nor(ii,jj,kk,:) 
            cur_v(:) = cur(ii,jj,kk,:)
            !
            ! ind2 = maxloc(abs(nor_v(:)),1)
            if(    abs(nor(ii,jj,kk,1)).eq.max(abs(nor(ii,jj,kk,1)),abs(nor(ii,jj,kk,2)),abs(nor(ii,jj,kk,3)))) then
              ind2 = 1
            elseif(abs(nor(ii,jj,kk,2)).eq.max(abs(nor(ii,jj,kk,1)),abs(nor(ii,jj,kk,2)),abs(nor(ii,jj,kk,3)))) then
              ind2 = 2
            elseif(abs(nor(ii,jj,kk,3)).eq.max(abs(nor(ii,jj,kk,1)),abs(nor(ii,jj,kk,2)),abs(nor(ii,jj,kk,3)))) then
              ind2 = 3
            end if
            !
            if(    ind2.eq.1) then
              !
              a    = xa
              b    = xb
              rm1  = 0.5_rp*((yb-ya)*gm+(ya+yb))
              rp1  = 0.5_rp*((yb-ya)*gp+(ya+yb))
              rm2  = 0.5_rp*((zb-za)*gm+(za+zb))
              rp2  = 0.5_rp*((zb-za)*gp+(za+zb))
              flux(i,j,k) = (yb-ya)*(zb-za)
              !
              xv(1, 1) = b; xv(1, 2) = rm1; xv(1, 3) = rm2
              xv(2, 1) = a; xv(2, 2) = rm1; xv(2, 3) = rm2
              xv(3, 1) = b; xv(3, 2) = rm1; xv(3, 3) = rp2
              xv(4, 1) = a; xv(4, 2) = rm1; xv(4, 3) = rp2
              xv(5, 1) = b; xv(5, 2) = rp1; xv(5, 3) = rm2
              xv(6, 1) = a; xv(6, 2) = rp1; xv(6, 3) = rm2
              xv(7, 1) = b; xv(7, 2) = rp1; xv(7, 3) = rp2
              xv(8, 1) = a; xv(8, 2) = rp1; xv(8, 3) = rp2
              !
              !xv(1,:) = (/b,rm1,rm2/) 
              !xv(2,:) = (/a,rm1,rm2/)
              !xv(3,:) = (/b,rm1,rp2/)
              !xv(4,:) = (/a,rm1,rp2/)
              !xv(5,:) = (/b,rp1,rm2/)
              !xv(6,:) = (/a,rp1,rm2/)
              !xv(7,:) = (/b,rp1,rp2/)
              !xv(8,:) = (/a,rp1,rp2/)
              !
            elseif(ind2.eq.2) then
              !
              a    = ya
              b    = yb
              rm1  = 0.5_rp*((xb-xa)*gm+(xa+xb))
              rp1  = 0.5_rp*((xb-xa)*gp+(xa+xb))
              rm2  = 0.5_rp*((zb-za)*gm+(za+zb))
              rp2  = 0.5_rp*((zb-za)*gp+(za+zb))
              flux(i,j,k) = (xb-xa)*(zb-za)
              !
              xv(1, 1) = rm1; xv(1, 2) = b; xv(1, 3) = rm2
              xv(2, 1) = rm1; xv(2, 2) = a; xv(2, 3) = rm2
              xv(3, 1) = rm1; xv(3, 2) = b; xv(3, 3) = rp2
              xv(4, 1) = rm1; xv(4, 2) = a; xv(4, 3) = rp2
              xv(5, 1) = rp1; xv(5, 2) = b; xv(5, 3) = rm2
              xv(6, 1) = rp1; xv(6, 2) = a; xv(6, 3) = rm2
              xv(7, 1) = rp1; xv(7, 2) = b; xv(7, 3) = rp2
              xv(8, 1) = rp1; xv(8, 2) = a; xv(8, 3) = rp2
              !
              !xv(1,:) = (/rm1,b,rm2/) 
              !xv(2,:) = (/rm1,a,rm2/)
              !xv(3,:) = (/rm1,b,rp2/)
              !xv(4,:) = (/rm1,a,rp2/)
              !xv(5,:) = (/rp1,b,rm2/)
              !xv(6,:) = (/rp1,a,rm2/)
              !xv(7,:) = (/rp1,b,rp2/)
              !xv(8,:) = (/rp1,a,rp2/)
              !
            elseif(ind2.eq.3) then
              !
              a    = za
              b    = zb
              rm1  = 0.5_rp*((xb-xa)*gm+(xa+xb))
              rp1  = 0.5_rp*((xb-xa)*gp+(xa+xb))
              rm2  = 0.5_rp*((yb-ya)*gm+(ya+yb))
              rp2  = 0.5_rp*((yb-ya)*gp+(ya+yb))
              flux(i,j,k) = (xb-xa)*(yb-ya)
              !
              xv(1, 1) = rm1; xv(1, 2) = rm2; xv(1, 3) = b
              xv(2, 1) = rm1; xv(2, 2) = rm2; xv(2, 3) = a
              xv(3, 1) = rm1; xv(3, 2) = rp2; xv(3, 3) = b
              xv(4, 1) = rm1; xv(4, 2) = rp2; xv(4, 3) = a
              xv(5, 1) = rp1; xv(5, 2) = rm2; xv(5, 3) = b
              xv(6, 1) = rp1; xv(6, 2) = rm2; xv(6, 3) = a
              xv(7, 1) = rp1; xv(7, 2) = rp2; xv(7, 3) = b
              xv(8, 1) = rp1; xv(8, 2) = rp2; xv(8, 3) = a
              !
              !xv(1,:) = (/rm1,rm2,b/) 
              !xv(2,:) = (/rm1,rm2,a/)
              !xv(3,:) = (/rm1,rp2,b/)
              !xv(4,:) = (/rm1,rp2,a/)
              !xv(5,:) = (/rp1,rm2,b/)
              !xv(6,:) = (/rp1,rm2,a/)
              !xv(7,:) = (/rp1,rp2,b/)
              !xv(8,:) = (/rp1,rp2,a/)
              !
            endif
            !
            cv(:)    = 1.0_rp
            cv(ind2) = 0.0_rp
            !
            cxxa = cv(1)*0.5_rp*cv(1)*cur_v(1)
            cyya = cv(2)*0.5_rp*cv(2)*cur_v(2)
            czza = cv(3)*0.5_rp*cv(3)*cur_v(3)
            cxya = cv(1)*cv(2)*cur_v(4)
            cxza = cv(1)*cv(3)*cur_v(5)
            cyza = cv(2)*cv(3)*cur_v(6)
            a100 = (nor_v(1)-0.5_rp*cv(1)*(cur_v(1)+cv(2)*cur_v(4)+cv(3)*cur_v(5)))
            a010 = (nor_v(2)-0.5_rp*cv(2)*(cur_v(2)+cv(1)*cur_v(4)+cv(3)*cur_v(6)))
            a001 = (nor_v(3)-0.5_rp*cv(3)*(cur_v(3)+cv(1)*cur_v(5)+cv(2)*cur_v(6)))
            !
            sum_int_func = 0.0_rp
            !$acc loop seq
            do p=1,7,2
              !
              q = p+1
              surf_a = cxxa*xv(p,1)*xv(p,1) + cyya*xv(p,2)*xv(p,2) + czza*xv(p,3)*xv(p,3) + &
                       cxya*xv(p,1)*xv(p,2) + cyza*xv(p,2)*xv(p,3) + cxza*xv(p,3)*xv(p,1) + &
                       a100*xv(p,1)+a010*xv(p,2)+a001*xv(p,3)
              surf_b = cxxa*xv(q,1)*xv(q,1) + cyya*xv(q,2)*xv(q,2) + czza*xv(q,3)*xv(q,3) + &
                       cxya*xv(q,1)*xv(q,2) + cyza*xv(q,2)*xv(q,3) + cxza*xv(q,3)*xv(q,1) + &
                       a100*xv(q,1)+a010*xv(q,2)+a001*xv(q,3)
              !
              sum_int_func = sum_int_func + &
                             (b-a+b_thi/(nor_v(ind2)+limit)*log( &
                             cosh(b_th*(surf_a+d_thinc(ii,jj,kk))) / &
                             cosh(b_th*(surf_b+d_thinc(ii,jj,kk)))))
              !
            enddo
            !$acc end loop
            !
            flux(i,j,k) = 0.125_rp*flux(i,j,k)*sum_int_func
            !
          endif
          !
          flux(i,j,k) = flux(i,j,k)*vel_sign*dl
          !
        enddo
      enddo
    enddo
    !$acc end parallel loop
    !
    return
  end subroutine cmpt_vof_flux
#endif
  !
  subroutine solve_quar_paper(a0,a1,a2,a3,a4,x)
    !
    ! subroutine to solve the quartic equation
    !
    implicit none
    !$acc routine(solve_quar_paper) seq
    !
    real(rp), intent(in ) :: a0,a1,a2,a3,a4
    real(rp), intent(out) :: x
    !
    real(rp) :: b3,b2,b1,b0,a4i
    real(rp) :: c2,c1,c0
    real(rp) :: z1,check,a,b
    real(rp) :: aa,bb,cc,dd,dd_s
    ! 
    ! First change the coefficients in the form of Eq. (B.2)
    !
    a4i = 1._rp/a4
    b3  = a3*a4i
    b2  = a2*a4i
    b1  = a1*a4i
    b0  = a0*a4i
    !
    ! Calculate the coefficients of Eq. (B.5)
    !  note: in c1 it is not b2 (as in the paper), but b3!
    !
    c2  = -b2
    c1  = b1*b3-4._rp*b0
    c0  = b0*(4._rp*b2-b3**2)-b1**2
    !
    ! Calculate z1, Eq. (B.7)
    !
    a   = -c2**2/9._rp+c1/3._rp
    b   = 2._rp*c2**3/27._rp-c1*c2/3._rp+c0
    check = b**2+4._rp*a**3
    !
    if(check.ge.0._rp) then
      z1 = sign(1._rp,(-b+sqrt(check)))*abs(0.5_rp*(-b+sqrt(check)))**(1._rp/3._rp) + &
           sign(1._rp,(-b-sqrt(check)))*abs(0.5_rp*(-b-sqrt(check)))**(1._rp/3._rp) - &
           c2/3._rp
    else
      z1 = 2._rp*sqrt(-a)*cos(atan(sqrt(-check)/(-b))/3._rp)-c2/3._rp
    endif
    !
    ! Find new coefficients, Eq. (B.9)
    !
    aa   = 0.5_rp*b3
    bb   = 0.5_rp*z1
    dd_s = bb**2-b0 ! as said in the paper, this should be always positive but still we put a check
    if(dd_s.le.0._rp) then
      dd = limit    ! to avoid singularity in case of imaginary solution (which we filter out)
    else
      dd = sqrt(dd_s)
    endif
    cc = (-0.5_rp*b1+aa*bb)/dd
    !
    ! Finally, calculate solution from Eq. (B.11)
    ! Be aware of a small error in the paper, i.e. 
    ! the (+) sign inside the square root should be (-)
    !
    x = 0.5_rp*(-(aa-cc) + sqrt((aa-cc)**2-4._rp*(bb-dd)))
    !
    return
  end subroutine solve_quar_paper
  !
  subroutine solve_quad(a0,a1,a2,x)
    !
    ! subroutine to solve the quadratic equation
    !
    implicit none
    !$acc routine(solve_quad) seq
    !
    real(rp), intent(in ) :: a0,a1,a2
    real(rp), intent(out) :: x
    !
    real(rp)              :: x1,x2
    ! 
    x1 = 0.5_rp*(-a1+sqrt(a1**2-4.0_rp*a2*a0))/a2
    x2 = 0.5_rp*(-a1-sqrt(a1**2-4.0_rp*a2*a0))/a2
    x  = max(x1,x2)
    !
    return
  end subroutine solve_quad
  !
  ! code devoted to the VoF initialization
  !
#if defined(_INIT_MONTECARLO)
  !--------------------------------------C interface--------------------------------------------
  interface
    subroutine init_MTHINC(xc, yc, zc, rs, xl, xu, dx, beta, vof) BIND(C, name="init_MTHINC_")
      use, intrinsic :: ISO_C_BINDING, ONLY: C_INT
      implicit none
      ! integer (C_INT) ::  dimM(0:2), Ng, maxit, it, ierr, norm,drank, dnProcNode
      real(rp) :: xc, yc, zc, rs,  dx, beta, vof
      real(rp) :: xl(3), xu(3)
      ! real(rp) :: maxError, beta, tres2
    end subroutine init_MTHINC
  end interface
#endif
  !
  subroutine initvof(n,dli,vof)
    !
    use mod_param     , only: inivof,lx,ly,lz,cbcvof,xc,yc,zc,r,nbub
    use mod_common_mpi, only: myid,ierr,ijk_start
    use mod_sanity    , only: flutas_error
    !
    ! computes initial conditions for the VoF field
    !
    implicit none
    !
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(out), dimension(0:,0:,0:) :: vof
    !
    integer  :: i,j,k,q,ii,jj,kk, i_b
    real(rp) :: x,y,z,xl,yl,zl,xx,yy,zz
    real(rp) :: sdist,sdistmin
    real(rp) :: zfilm_top,zfilm_bot,sdist1,sdist2,dfilm
    real(rp) :: xfilm_bot,yfilm_bot,csi1,csi2
    integer  :: nbox
    real(rp) :: eps,epsbox
    real(rp), dimension(3) :: dl,dlbox
    integer , dimension(3) :: iperiod
    real(rp) :: grid_vol_ratio
    integer , dimension(2):: nx_b, ny_b, nz_b
    integer :: n1, n2, n3
    real(rp):: sign_v
    !@cuf attributes(managed) :: vof
    !
    nbox = 100
    dl(:) = dli(:)**(-1)
    dlbox(:) = dl(:)/(1.0_rp*nbox)
    !
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)
    !
#if defined(_TWOD)
    eps    = sqrt(              dl(2)**2    + dl(3)**2   )*0.5_rp
    epsbox = sqrt(              dlbox(2)**2 + dlbox(3)**2)*0.5_rp
#else
    eps    = sqrt(dl(1)**2    + dl(2)**2    + dl(3)**2   )*0.5_rp
    epsbox = sqrt(dlbox(1)**2 + dlbox(2)**2 + dlbox(3)**2)*0.5_rp
#endif
    grid_vol_ratio = product(dlbox(:))/product(dl(:))
    iperiod(:) = 0
    !
    do q=1,3
      if(cbcvof(0,q)//cbcvof(1,q).eq.'PP') iperiod(q) = 1
    enddo
    !
    !$acc kernels
    vof(:,:,:) = 0.0_rp
    !$acc end kernels
    !
    select case(trim(inivof))
    case('bub')
      ! if (nbub.ne.1) then
      !   print*, "Error in the number of bubble nbub. Use VOFI (bub_vofi) for nbub>1!"
      !   call exit
      ! endif
      !
      do i_b=1,nbub
        call bub_lim(nx_b,xc(i_b),r(i_b),dl(1),n(1),ijk_start(1))
        call bub_lim(ny_b,yc(i_b),r(i_b),dl(2),n(2),ijk_start(2))
        call bub_lim(nz_b,zc(i_b),r(i_b),dl(3),n(3),ijk_start(3))
        !$acc kernels
        do k=nz_b(1),nz_b(2)
          z = (k+ijk_start(3)-0.5_rp)*dl(3)
          do j=ny_b(1),ny_b(2)
            y = (j+ijk_start(2)-0.5_rp)*dl(2)
            do i=nx_b(1),nx_b(2)
              if (vof(i,j,k).eq.0.) then
                x = (i+ijk_start(1)-0.5_rp)*dl(1)
#if defined(_TWOD)
                sdistmin = max(   ly,lz)*2.0_rp
#else
                sdistmin = max(lx,ly,lz)*2.0_rp
#endif
                do kk = -1,1
                  do jj = -1,1
                    do ii = -1,1
#if defined(_TWOD)
                      sdist = sqrt(  &
                                    (y+jj*iperiod(2)*ly-yc(i_b))**2 + &
                                    (z+kk*iperiod(3)*lz-zc(i_b))**2 ) - r(i_b)
#else    
                      sdist = sqrt( (x+ii*iperiod(1)*lx-xc(i_b))**2 + &
                                    (y+jj*iperiod(2)*ly-yc(i_b))**2 + &
                                    (z+kk*iperiod(3)*lz-zc(i_b))**2 ) - r(i_b)
#endif
                        if(abs(sdist).lt.sdistmin) sdistmin = sdist
                    enddo
                  enddo
                enddo
                sdist = sdistmin
                if(     sdist.lt.-eps ) then
                  vof(i,j,k) = 1.0_rp
                elseif( sdist.gt. eps ) then
                  vof(i,j,k) = 0.0_rp
                else
#if !defined(_INIT_MONTECARLO)
                  zl = z-dl(3)/2.0_rp
                  yl = y-dl(2)/2.0_rp
                  xl = x-dl(1)/2.0_rp
                  do kk=1,nbox
                    zz = zl + (kk-0.5_rp)*dlbox(3)
                    do jj=1,nbox
                      yy = yl + (jj-0.5_rp)*dlbox(2)
                      do ii=1,nbox
                        xx = xl + (ii-0.5_rp)*dlbox(1)
#if defined(_TWOD)
                        sdist = sqrt(                  (yy-yc(i_b))**2 + (zz-zc(i_b))**2) - r(i_b)
#else
                        sdist = sqrt((xx-xc(i_b))**2 + (yy-yc(i_b))**2 + (zz-zc(i_b))**2) - r(i_b)
#endif
                        if(sdist.lt.-epsbox) vof(i,j,k) = vof(i,j,k) + grid_vol_ratio
                      enddo
                    enddo
                  enddo
#else
                  call init_MTHINC(xc(i_b), yc(i_b), zc(i_b), r(i_b), &
                    (/x-dl(1)*0.5_rp, y-dl(2)*0.5_rp, z-dl(3)*0.5_rp/), &
                    (/x+dl(1)*0.5_rp, y+dl(2)*0.5_rp, z+dl(3)*0.5_rp/), & 
                    dl(1), beta_thinc, vof(i,j,k))
                  vof(i,j,k) = vof(i,j,k)*dli(1)*dli(2)*dli(3)
#endif
                endif
              endif
            enddo
          enddo
        enddo
        !$acc end kernels
        !
      enddo
      !
    case('uni')
      !$acc kernels
      do k=1,n3
        do j=1,n2
          do i=1,n1
            vof(i,j,k) = 1.0_rp
           enddo
         enddo
       enddo
       !$acc end kernels
       !
       !
    case('zer')
      !$acc kernels
      do k=1,n3
        do j=1,n2
          do i=1,n1
            vof(i,j,k) = 0.0_rp
          enddo
        enddo
      enddo
      !$acc end kernels
      !
      !
    case('flm')
      eps    = dl(3)
      epsbox = dlbox(3)
      grid_vol_ratio = product(dlbox(:))/product(dl(:))
      dfilm = 1._rp/12._rp*.5_rp
      zfilm_top = lz/2._rp+dfilm/2._rp
      zfilm_bot = lz/2._rp-dfilm/2._rp
      do k=1,n3
        z = (k-0.5_rp)*dl(3)
        do j=1,n2
          do i=1,n1
            sdist1 =  (z - zfilm_top)
            sdist2 = -(z - zfilm_bot)
            if(     all((/sdist1,sdist2/) .lt.-eps) ) then
              vof(i,j,k) = 1.0_rp
            elseif( all((/sdist1,sdist2/) .gt. eps) ) then
              vof(i,j,k) = 0.0_rp
            else
              zl = z-dl(3)/2.0_rp
              do kk=1,nbox
                zz = zl + kk*dlbox(3)
                do jj=1,nbox
                  do ii=1,nbox
                    sdist1 =  (zz - zfilm_top)
                    sdist2 = -(zz - zfilm_bot)
                    if( all((/sdist1,sdist2/).lt.-epsbox) ) vof(i,j,k) = vof(i,j,k) + grid_vol_ratio
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
      !
    case('twx')
      eps    = dl(1)
      epsbox = dlbox(1)
      grid_vol_ratio = product(dlbox(:))/product(dl(:))
      xfilm_bot = lx/2.0_rp
      !$acc kernels
      do k=1,n3
        do j=1,n2
          do i=1,n1
            ii = ijk_start(1) + i
            x = (ii-0.5_rp)*dl(1)
            sdist1 =  (x - xfilm_bot)
            if(     sdist1 .gt. +eps)  then
              vof(i,j,k) = 1.0_rp
            elseif( sdist1 .lt. -eps)  then
              vof(i,j,k) = 0.0_rp
            else
              xl = x-dl(1)/2.0_rp
              do ii=1,nbox
                xx = xl + ii*dlbox(1)
                do jj=1,nbox
                  do kk=1,nbox
                    sdist1 =  (xx - xfilm_bot)
                    if( sdist1.gt.+epsbox)  vof(i,j,k) = vof(i,j,k) + grid_vol_ratio
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
      !$acc end kernels
      !
    case('twy')
      eps    = dl(2)
      epsbox = dlbox(2)
      grid_vol_ratio = product(dlbox(:))/product(dl(:))
      yfilm_bot = ly/2.0_rp
      !$acc kernels
      do j=1,n2
        jj = ijk_start(2) + j
        y = (jj-0.5_rp)*dl(2)
        do k=1,n3
          do i=1,n1
            sdist1 =  (y - yfilm_bot)
            if(     sdist1 .gt. +eps)  then
              vof(i,j,k) = 1.0_rp
            elseif( sdist1 .lt. -eps)  then
              vof(i,j,k) = 0.0_rp
            else
              yl = y-dl(2)/2.0_rp
              do jj=1,nbox
                yy = yl + jj*dlbox(2)
                do kk=1,nbox
                  do ii=1,nbox
                    sdist1 =  (yy - yfilm_bot)
                    if( sdist1.gt.+epsbox)  vof(i,j,k) = vof(i,j,k) + grid_vol_ratio
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
      !$acc end kernels
      !
    case('twz')
      eps    = dl(3)
      epsbox = dlbox(3)
      grid_vol_ratio = product(dlbox(:))/product(dl(:))
      zfilm_bot = lz/2.0_rp
      !$acc kernels
      do k=1,n3
        kk = ijk_start(3) + k
        z = (kk-0.5_rp)*dl(3)
        do j=1,n2
          do i=1,n1
            sdist1 =  (z - zfilm_bot)
            if(     sdist1 .gt. +eps)  then
              vof(i,j,k) = 1.0_rp
            elseif( sdist1 .lt. -eps)  then
              vof(i,j,k) = 0.0_rp
            else
              zl = z-dl(3)/2.0_rp
              do kk=1,nbox
                zz = zl + kk*dlbox(3)
                do jj=1,nbox
                  do ii=1,nbox
                    sdist1 =  (zz - zfilm_bot)
                    if( sdist1.gt.+epsbox)  vof(i,j,k) = vof(i,j,k) + grid_vol_ratio
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
      !$acc end kernels
      !
    case('tax')
      !
      ! note: sign_v = 1 (vof=1 top, vof=0 bottom), sign_v = -1 (vof=0 top, vof=1 bottom)
      !
      sign_v    = +1._rp
      xfilm_bot = lx/2.0_rp
      !$acc kernels
      do k=1,n3
        do j=1,n2
          do i=1,n1
            ii = ijk_start(1) + i
            x  = (ii-0.5_rp)*dl(1)
            !
            if(     x.gt.(xfilm_bot+0.05_rp*lx))  then
              vof(i,j,k) = (1._rp+1._rp*sign_v)/2._rp
            elseif( x.lt.(xfilm_bot-0.05_rp*lx))  then
              vof(i,j,k) = (1._rp-1._rp*sign_v)/2._rp
            else
              csi1 = (b_th/dl(1))*(sign_v*(x-dl(1)*0.5_rp-xfilm_bot))
              csi2 = (b_th/dl(1))*(sign_v*(x+dl(1)*0.5_rp-xfilm_bot))
              vof(i,j,k) = 1._rp*sign_v*((0.5_rp*csi2+0.5_rp*log(cosh(csi2))) - &
                                         (0.5_rp*csi1+0.5_rp*log(cosh(csi1))))
              vof(i,j,k) = vof(i,j,k)*dl(1)/(dl(1)*b_th)
            endif
            !
            vof(i,j,k) = min(max(0.0_rp,vof(i,j,k)),1.0_rp) ! also in Matlab, noticed overshoot of 1.0+1e-12
            !
          enddo
        enddo
      enddo
      !$acc end kernels
      !
    case('tay')
      !
      ! note: sign_v = 1 (vof=1 top, vof=0 bottom), sign_v = -1 (vof=0 top, vof=1 bottom)
      !
      sign_v    = +1._rp
      yfilm_bot = ly/2.0_rp
      !$acc kernels
      do k=1,n3
        do j=1,n2
          do i=1,n1
            jj = ijk_start(2) + j
            y  = (jj-0.5_rp)*dl(2)
            !
            if(    y.gt.(yfilm_bot+0.05_rp*ly))  then
              vof(i,j,k) = (1._rp+1._rp*sign_v)/2._rp
            elseif(y.lt.(yfilm_bot-0.05_rp*ly))  then
              vof(i,j,k) = (1._rp-1._rp*sign_v)/2._rp
            else
              csi1 = (b_th/dl(2))*(sign_v*(y-dl(2)*0.5_rp-yfilm_bot))
              csi2 = (b_th/dl(2))*(sign_v*(y+dl(2)*0.5_rp-yfilm_bot))
              vof(i,j,k) = 1._rp*sign_v*((0.5_rp*csi2+0.5_rp*log(cosh(csi2))) - &
                                         (0.5_rp*csi1+0.5_rp*log(cosh(csi1))))
              vof(i,j,k) = vof(i,j,k)*dl(2)/(dl(2)*b_th)
            endif
            !
            vof(i,j,k) = min(max(0.0_rp,vof(i,j,k)),1.0_rp) ! also in Matlab, noticed overshoot of 1.0+1e-12
            !
          enddo
        enddo
      enddo
      !$acc end kernels
      !
    case('taz')
      !
      ! note: sign_v = 1 (vof=1 top, vof=0 bottom), sign_v = -1 (vof=0 top, vof=1 bottom)
      !
      sign_v    = +1._rp
      zfilm_bot = lz/2.0_rp
      !$acc kernels
      do k=1,n3
        do j=1,n2
          do i=1,n1
            kk = ijk_start(3) + k
            z  = (kk-0.5)*dl(3)
            !
            if(    z.gt.(zfilm_bot+0.05_rp*lz))  then
              vof(i,j,k) = (1._rp+1._rp*sign_v)/2._rp
            elseif(z.lt.(zfilm_bot-0.05_rp*lz))  then
              vof(i,j,k) = (1._rp-1._rp*sign_v)/2._rp
            else
              csi1 = (b_th/dl(3))*(sign_v*(z-dl(3)*0.5-zfilm_bot))*1._rp
              csi2 = (b_th/dl(3))*(sign_v*(z+dl(3)*0.5-zfilm_bot))*1._rp
              vof(i,j,k) = 1._rp*sign_v*((0.5_rp*csi2+0.5_rp*log(cosh(csi2))) - &
                                         (0.5_rp*csi1+0.5_rp*log(cosh(csi1))))
              vof(i,j,k) = vof(i,j,k)*dl(3)/(dl(3)*b_th)
            endif
            !
            vof(i,j,k) = min(max(0._rp,vof(i,j,k)),1._rp) ! also in Matlab, noticed overshoot of 1.0+1e-12
            !
          enddo
        enddo
      enddo
      !$acc end kernels
      !
    case default  
      call flutas_error('Error: invalid name of the initial VoF field. Simulation aborted. Check vof.in')
    end select
    return
  end subroutine initvof
  !
  subroutine bub_lim(nx_b,xc,r,dx,n,ijk_start)
    !
    implicit none
    !
    integer , intent(in ) :: n,ijk_start
    real(rp), intent(in ) :: xc,r,dx
    integer , intent(out), dimension(2) :: nx_b
    !
    integer :: lb, ub
    !
    lb = floor((xc-r)/dx)-1  
    if    (lb.le.(ijk_start+1)) then
      nx_b(1) = 1
    elseif(lb.ge.(ijk_start+n)) then
      nx_b(1) = n
    else
      nx_b(1) = lb-ijk_start 
    endif
    !
    ub = ceiling((xc+r)/dx)+1
    if    (ub.le.(ijk_start+1)) then
     nx_b(2)  = 1
    elseif(ub.ge.(ijk_start+n)) then
      nx_b(2) = n
    else
      nx_b(2) = ub-ijk_start 
    endif
    !
    return
  end subroutine bub_lim
  !
end module mod_vof
