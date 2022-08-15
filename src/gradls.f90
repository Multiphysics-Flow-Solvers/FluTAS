!
! SPDX-License-Identifier: MIT
!
module mod_gradls
  !
  use mod_types, only: rp
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: weno5,weno5_old
  !
  contains
  !
  subroutine weno5_old(n,dli,qmin,is_f,phi,ux,uy,uz,dphidt)
    !
    implicit none
    !
    real(rp), parameter, dimension(3,3) :: c = 1._rp/6._rp*reshape((/ 2._rp,-7._rp,11._rp, &
                                                                     -1._rp, 5._rp, 2._rp, &
                                                                      2._rp, 5._rp,-1._rp/),shape(c))
    real(rp), parameter, dimension(3)   :: sigma = (/0.1_rp,0.6_rp,0.3_rp/) 
    real(rp), parameter :: eps = 10._rp**(-6)
    !
    integer , intent(in ), dimension(3)                    :: n
    real(rp), intent(in ), dimension(3)                    :: dli
    integer , intent(in )                                  :: qmin
    logical , intent(in )                                  :: is_f
    real(rp), intent(in ), dimension(   -2:,   -2:,   -2:) :: phi
    real(rp), intent(in ), dimension(-qmin:,-qmin:,-qmin:) :: ux,uy,uz
    real(rp), intent(out), dimension(    1:,    1:,    1:) :: dphidt
    !
    real(rp), dimension(-2:2) :: f
    real(rp), dimension(3) :: beta,we,dfdlh
    real(rp) :: dphidx,dphidy,dphidz
    real(rp) :: uxc,uyc,uzc
    integer  :: a,i,j,k,p,q
    !
    if(is_f) then
      q = 1
    else
      q = 0
    endif
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          uxc = 0.5_rp*(ux(i-q,j,k)+ux(i,j,k))
#if defined(_TWOD)
          dphidx = 0._rp
#else
          a = nint(sign(1._rp,uxc))
          f(-2) = a*(phi(i-2*a,j,k) - phi(i-3*a,j,k))*dli(1)
          f(-1) = a*(phi(i-1*a,j,k) - phi(i-2*a,j,k))*dli(1)
          f( 0) = a*(phi(i+0*a,j,k) - phi(i-1*a,j,k))*dli(1)
          f( 1) = a*(phi(i+1*a,j,k) - phi(i+0*a,j,k))*dli(1)
          f( 2) = a*(phi(i+2*a,j,k) - phi(i+1*a,j,k))*dli(1)
          beta(1) = (13._rp/12._rp)*(      f(-2) -2._rp*f(-1) +       f( 0))**2 + &
                    ( 1._rp/4._rp )*(      f(-2) -4._rp*f(-1) + 3._rp*f( 0))**2
          beta(2) = (13._rp/12._rp)*(      f(-1) -2._rp*f( 0) +       f( 1))**2 + &
                    ( 1._rp/4._rp )*(      f(-1)              -       f( 1))**2
          beta(3) = (13._rp/12._rp)*(      f( 0) -2._rp*f( 1) +       f( 2))**2 + &
                    ( 1._rp/4._rp )*(3._rp*f( 0) -4._rp*f( 1) +       f( 2))**2
          !we(:) = sigma(:)/(beta(:)+eps)**2
          we(:) = sigma(:)*(1._rp+(abs(beta(1)-beta(3))/(eps+beta(:))))
          we(:) = we(:)/sum(we(:))
          dfdlh(1) = sum(c(:,1)*f(-2:0))
          dfdlh(2) = sum(c(:,2)*f(-1:1))
          dfdlh(3) = sum(c(:,3)*f( 0:2))
          dphidx   = sum(we(:)*dfdlh(:))
#endif
          !
          uyc = 0.5_rp*(uy(i,j-q,k)+uy(i,j,k))
          a = nint(sign(1._rp,uyc))
          f(-2) = a*(phi(i,j-2*a,k) - phi(i,j-3*a,k))*dli(2)
          f(-1) = a*(phi(i,j-1*a,k) - phi(i,j-2*a,k))*dli(2)
          f( 0) = a*(phi(i,j+0*a,k) - phi(i,j-1*a,k))*dli(2)
          f( 1) = a*(phi(i,j+1*a,k) - phi(i,j+0*a,k))*dli(2)
          f( 2) = a*(phi(i,j+2*a,k) - phi(i,j+1*a,k))*dli(2)
          beta(1) = (13._rp/12._rp)*(      f(-2) - 2._rp*f(-1) +       f( 0))**2 + &
                    ( 1._rp/4._rp )*(      f(-2) - 4._rp*f(-1) + 3._rp*f( 0))**2
          beta(2) = (13._rp/12._rp)*(      f(-1) - 2._rp*f( 0) +       f( 1))**2 + &
                    ( 1._rp/4._rp )*(      f(-1)               -       f( 1))**2
          beta(3) = (13._rp/12._rp)*(      f( 0) - 2._rp*f( 1) +       f( 2))**2 + &
                    ( 1._rp/4._rp )*(3._rp*f( 0) - 4._rp*f( 1) +       f( 2))**2
          !we(:) = sigma(:)/(beta(:)+eps)**2
          we(:) = sigma(:)*(1._rp+(abs(beta(1)-beta(3))/(eps+beta(:))))
          we(:) = we(:)/sum(we(:))
          dfdlh(1) = sum(c(:,1)*f(-2:0))
          dfdlh(2) = sum(c(:,2)*f(-1:1))
          dfdlh(3) = sum(c(:,3)*f( 0:2))
          dphidy   = sum(we(:)*dfdlh(:))
          !
          uzc = 0.5_rp*(uz(i,j,k-q)+uz(i,j,k))
          a = nint(sign(1._rp,uzc))
          f(-2) = a*(phi(i,j,k-2*a) - phi(i,j,k-3*a))*dli(3)
          f(-1) = a*(phi(i,j,k-1*a) - phi(i,j,k-2*a))*dli(3)
          f( 0) = a*(phi(i,j,k+0*a) - phi(i,j,k-1*a))*dli(3)
          f( 1) = a*(phi(i,j,k+1*a) - phi(i,j,k+0*a))*dli(3)
          f( 2) = a*(phi(i,j,k+2*a) - phi(i,j,k+1*a))*dli(3)
          beta(1) = (13._rp/12._rp)*(      f(-2) - 2._rp*f(-1) +       f( 0))**2 + &
                    ( 1._rp/4._rp )*(      f(-2) - 4._rp*f(-1) + 3._rp*f( 0))**2
          beta(2) = (13._rp/12._rp)*(      f(-1) - 2._rp*f( 0) +       f( 1))**2 + &
                    ( 1._rp/4._rp )*(      f(-1)               -       f( 1))**2
          beta(3) = (13._rp/12._rp)*(      f( 0) - 2._rp*f( 1) +       f( 2))**2 + &
                    ( 1._rp/4._rp )*(3._rp*f( 0) - 4._rp*f( 1) +       f( 2))**2
          !we(:) = sigma(:)/(beta(:)+eps)**2
          we(:) = sigma(:)*(1._rp+(abs(beta(1)-beta(3))/(eps+beta(:))))
          we(:) = we(:)/sum(we(:))
          dfdlh(1) = sum(c(:,1)*f(-2:0))
          dfdlh(2) = sum(c(:,2)*f(-1:1))
          dfdlh(3) = sum(c(:,3)*f( 0:2))
          dphidz   = sum(we(:)*dfdlh(:))
          !
          dphidt(i,j,k) = - (uxc*dphidx + uyc*dphidy + uzc*dphidz)
          !
        enddo
      enddo
    enddo
    return
  end subroutine weno5_old
  !
  subroutine weno5(nx,ny,nz,dxi,dyi,dzi,nh_u,is_f,phi,ux,uy,uz,dphidt)
    !
    implicit none
    !
    real(rp), parameter :: c11    =  2._rp/6._rp, &
                           c21    = -7._rp/6._rp, & 
                           c31    = 11._rp/6._rp, & 
                           c12    = -1._rp/6._rp, & 
                           c22    =  5._rp/6._rp, & 
                           c32    =  2._rp/6._rp, & 
                           c13    =  2._rp/6._rp, & 
                           c23    =  5._rp/6._rp, & 
                           c33    = -1._rp/6._rp 
    real(rp), parameter :: sigma1 = 0.1_rp, &
                           sigma2 = 0.6_rp, & 
                           sigma3 = 0.3_rp 
    real(rp), parameter :: eps    = 10._rp**(-6)
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_u
    logical , intent(in )                                     :: is_f
    real(rp), intent(in ), dimension(    -2:,    -2:,    -2:) :: phi
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: ux,uy,uz
    real(rp), intent(out), dimension(      :,      :,      :) :: dphidt
    !
    real(rp) :: fm2,fm1,f0,fp1,fp2
    real(rp) :: beta1,beta2,beta3
    real(rp) :: we1,we2,we3,sum_we
    real(rp) :: dfdlh1,dfdlh2,dfdlh3
    real(rp) :: dphidx,dphidy,dphidz
    integer  :: a,i,j,k,p
    real(rp) :: uxc,uyc,uzc
    integer  :: q
    !@cuf attributes(managed) :: ux, uy, uz, phi, dphidt
    !
    if(is_f) then
      q = 1
    else
      q = 0
    endif
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(uxc,uyc,uzc,a) &
    !$OMP PRIVATE(fm2,fm1,f0,fp1,fp2) &
    !$OMP PRIVATE(beta1,beta2,beta3,we1,we2,we3,eps) &
    !$OMP PRIVATE(c11,c21,c31,c12,c22,c32,c13,c23,c33) &
    !$OMP PRIVATE(dfdlh1,dfdlh2,dfdlh3,dphidx,dphidy,dphidz) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,ux,uy,uz,phi,dphidt)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          uxc = 0.5_rp*(ux(i-q,j,k)+ux(i,j,k))
#if defined(_TWOD)
          dphidx = 0._rp
#else
          a   = nint(sign(1._rp,uxc))
          fm2 = a*(phi(i-2*a,j,k) - phi(i-3*a,j,k))*dxi
          fm1 = a*(phi(i-1*a,j,k) - phi(i-2*a,j,k))*dxi
          f0  = a*(phi(i+0*a,j,k) - phi(i-1*a,j,k))*dxi
          fp1 = a*(phi(i+1*a,j,k) - phi(i+0*a,j,k))*dxi
          fp2 = a*(phi(i+2*a,j,k) - phi(i+1*a,j,k))*dxi
          beta1 = (13._rp/12._rp)*(      fm2 - 2._rp*fm1 +       f0 )**2 + &
                  ( 1._rp/4._rp )*(      fm2 - 4._rp*fm1 + 3._rp*f0 )**2
          beta2 = (13._rp/12._rp)*(      fm1 - 2._rp*f0  +       fp1)**2 + &
                  ( 1._rp/4._rp )*(      fm1             -       fp1)**2
          beta3 = (13._rp/12._rp)*(      f0  - 2._rp*fp1 +       fp2)**2 + &
                  ( 1._rp/4._rp )*(3._rp*f0  - 4._rp*fp1 +       fp2)**2
          we1 = sigma1*(1._rp+(abs(beta1-beta3)/(eps+beta1)))
          we2 = sigma2*(1._rp+(abs(beta1-beta3)/(eps+beta2)))
          we3 = sigma3*(1._rp+(abs(beta1-beta3)/(eps+beta3)))
          sum_we = we1+we2+we3
          we1 = we1/sum_we
          we2 = we2/sum_we
          we3 = we3/sum_we
          dfdlh1 = c11*fm2+c21*fm1+c31*f0
          dfdlh2 = c12*fm1+c22*f0 +c32*fp1
          dfdlh3 = c13*f0 +c23*fp1+c33*fp2
          dphidx = we1*dfdlh1+we2*dfdlh2+we3*dfdlh3
#endif
          !
          uyc = 0.5_rp*(uy(i,j-q,k)+uy(i,j,k))
          a   = nint(sign(1._rp,uyc))
          fm2 = a*(phi(i,j-2*a,k) - phi(i,j-3*a,k))*dyi
          fm1 = a*(phi(i,j-1*a,k) - phi(i,j-2*a,k))*dyi
          f0  = a*(phi(i,j+0*a,k) - phi(i,j-1*a,k))*dyi
          fp1 = a*(phi(i,j+1*a,k) - phi(i,j+0*a,k))*dyi
          fp2 = a*(phi(i,j+2*a,k) - phi(i,j+1*a,k))*dyi
          beta1 = (13._rp/12._rp)*(      fm2 - 2._rp*fm1 +       f0 )**2 + &
                  ( 1._rp/4._rp )*(      fm2 - 4._rp*fm1 + 3._rp*f0 )**2
          beta2 = (13._rp/12._rp)*(      fm1 - 2._rp*f0  +       fp1)**2 + &
                  ( 1._rp/4._rp )*(      fm1             -       fp1)**2
          beta3 = (13._rp/12._rp)*(      f0  - 2._rp*fp1 +       fp2)**2 + &
                  ( 1._rp/4._rp )*(3._rp*f0  - 4._rp*fp1 +       fp2)**2
          we1 = sigma1*(1._rp+(abs(beta1-beta3)/(eps+beta1)))
          we2 = sigma2*(1._rp+(abs(beta1-beta3)/(eps+beta2)))
          we3 = sigma3*(1._rp+(abs(beta1-beta3)/(eps+beta3)))
          sum_we = we1+we2+we3
          we1 = we1/sum_we
          we2 = we2/sum_we
          we3 = we3/sum_we
          dfdlh1 = c11*fm2+c21*fm1+c31*f0
          dfdlh2 = c12*fm1+c22*f0 +c32*fp1
          dfdlh3 = c13*f0 +c23*fp1+c33*fp2
          dphidy = we1*dfdlh1+we2*dfdlh2+we3*dfdlh3
          !
          uzc = 0.5_rp*(uz(i,j,k-q)+uz(i,j,k))
          a   = nint(sign(1._rp,uzc))
          fm2 = a*(phi(i,j,k-2*a) - phi(i,j,k-3*a))*dzi
          fm1 = a*(phi(i,j,k-1*a) - phi(i,j,k-2*a))*dzi
          f0  = a*(phi(i,j,k+0*a) - phi(i,j,k-1*a))*dzi
          fp1 = a*(phi(i,j,k+1*a) - phi(i,j,k+0*a))*dzi
          fp2 = a*(phi(i,j,k+2*a) - phi(i,j,k+1*a))*dzi
          beta1 = (13._rp/12._rp)*(      fm2 - 2._rp*fm1 +       f0 )**2 + &
                  ( 1._rp/4._rp )*(      fm2 - 4._rp*fm1 + 3._rp*f0 )**2
          beta2 = (13._rp/12._rp)*(      fm1 - 2._rp*f0  +       fp1)**2 + &
                  ( 1._rp/4._rp )*(      fm1             -       fp1)**2
          beta3 = (13._rp/12._rp)*(      f0  - 2._rp*fp1 +       fp2)**2 + &
                  ( 1._rp/4._rp )*(3._rp*f0  - 4._rp*fp1 +       fp2)**2
          we1 = sigma1*(1._rp+(abs(beta1-beta3)/(eps+beta1)))
          we2 = sigma2*(1._rp+(abs(beta1-beta3)/(eps+beta2)))
          we3 = sigma3*(1._rp+(abs(beta1-beta3)/(eps+beta3)))
          sum_we = we1+we2+we3
          we1 = we1/sum_we
          we2 = we2/sum_we
          we3 = we3/sum_we
          dfdlh1 = c11*fm2+c21*fm1+c31*f0
          dfdlh2 = c12*fm1+c22*f0 +c32*fp1
          dfdlh3 = c13*f0 +c23*fp1+c33*fp2
          dphidz = we1*dfdlh1+we2*dfdlh2+we3*dfdlh3
          !
          dphidt(i,j,k) = - (uxc*dphidx + uyc*dphidy + uzc*dphidz)
          !
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end kernels 
#else
    !$OMP END PARALLEL DO
#endif
    !
    return
  end subroutine weno5
end module mod_gradls
