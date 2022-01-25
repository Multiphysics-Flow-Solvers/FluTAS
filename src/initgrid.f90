!
! SPDX-License-Identifier: MIT
!
module mod_initgrid
  !
  use mod_param, only: pi
  use mod_types
  !
  implicit none
  !
  private
  public  :: initgrid
  !
  contains
  !
  ! [NOTE] This could be done on GPU but is it worth it?
  subroutine initgrid(inivel,n,gr,lz,nh_d,dzc,dzf,zc,zf)
    !
    ! initializes the non-uniform grid in z
    !
    implicit none
    !
    character(len=3), intent(in )                           :: inivel
    integer         , intent(in )                           :: n
    real(rp)        , intent(in )                           :: gr,lz
    integer         , intent(in )                           :: nh_d
    real(rp)        , intent(out), dimension(1-nh_d:n+nh_d) :: dzc,dzf,zc,zf
    !
    real(rp) :: z0
    integer  :: k
    procedure (), pointer :: gridpoint => null()
    select case(inivel)
    case('zer','log','poi','cou')
      gridpoint => gridpoint_cluster_two_end
    case('hcl','hcp')
      gridpoint => gridpoint_cluster_one_end
    case default
      gridpoint => gridpoint_cluster_two_end
    end select
    !
    ! step 1) determine coordinates of cell faces zf
    !
    do k=1,n
      z0  = (k-0.0_rp)/(1.0_rp*n)
      call gridpoint(gr,z0,zf(k))
      zf(k) = zf(k)*lz
    enddo
    zf(0) = 0.0_rp
    !
    ! step 2) determine grid spacing between faces dzf
    !
    do k=1,n
      dzf(k) = zf(k)-zf(k-1)
    enddo
    dzf(0  ) = dzf(1)
    dzf(n+1) = dzf(n)
    !
    ! step 3) determine grid spacing between centers dzc
    !
    do k=0,n
      dzc(k) = 0.5_rp*(dzf(k)+dzf(k+1))
    enddo
    dzc(n+1) = dzc(n)
    !
    ! step 4) compute coordinates of cell centers zc and faces zf
    !
    zc(0) = -dzc(0)/2.0_rp
    zf(0) = 0.0_rp
    do k=1,n+1
      zc(k) = zc(k-1) + dzc(k-1)
      zf(k) = zf(k-1) + dzf(k)
    enddo
    !
    ! step 5) extension to 0,-1,... and n+1,n+2,... for dzf and dzc
    !
    do k=1-nh_d,0
      dzf(k) = dzf(-k+1)
      dzc(k) = dzc(-k  )
    enddo
    do k=n+1,n+nh_d
      dzf(k) = dzf(2*n-k-1)
      dzc(k) = dzc(2*n-k  )
    enddo
    !
    ! step 6) extension to 0,-1,... and n+1,n+2,... for zf and zc
    !
    do k=0,1-nh_d,-1
      zf(k) = zf(k+1)-dzf(k  )
      zc(k) = zf(k+1)-dzc(k+1)
    enddo
    do k=n+1,n+nh_d
      zf(k) = zf(k-1)+dzf(k  )
      zc(k) = zf(k-1)+dzc(k-1)
    enddo
    !
    return
  end subroutine initgrid
  !
  ! grid stretching functions 
  ! see e.g., Fluid Flow Phenomena -- A Numerical Toolkit, by P. Orlandi 
  !
  subroutine gridpoint_cluster_two_end(alpha,z0,z)
    !
    ! clustered at the two sides
    !
    implicit none
    !
    real(rp), intent(in ) :: alpha,z0
    real(rp), intent(out) :: z
    !
    if(alpha.ne.0._rp) then
      z = 0.5_rp*(1._rp+tanh((z0-0.5_rp)*alpha)/tanh(alpha/2.0_rp))
    else
      z = z0
    endif
    !
    return
  end subroutine gridpoint_cluster_two_end
  !
  subroutine gridpoint_cluster_one_end(alpha,z0,z)
    !
    ! clustered at the lower side
    !
    implicit none
    !
    real(rp), intent(in ) :: alpha,z0
    real(rp), intent(out) :: z
    !
    if(alpha.ne.0.0_rp) then
      z = 1.0_rp*(1.0_rp+tanh((z0-1.0_rp)*alpha)/tanh(alpha/1.0_rp))
    else
      z = z0
    endif
    !
    return
  end subroutine gridpoint_cluster_one_end
  !
  subroutine gridpoint_cluster_middle(alpha,z0,z)
    !
    ! clustered in the middle
    !
    implicit none
    !
    real(rp), intent(in ) :: alpha,z0
    real(rp), intent(out) :: z
    !
    if(alpha.ne.0.0_rp) then
      if(    z0.le.0.5_rp) then 
        z = 0.5_rp*(1.0_rp-1.0_rp+tanh(2.0_rp*alpha*(z0-0.0_rp))/tanh(alpha))
      elseif(z0.gt.0.5) then
        z = 0.5_rp*(1.0_rp+1.0_rp+tanh(2.0_rp*alpha*(z0-1.0_rp))/tanh(alpha))
      endif
    else
      z = z0
    endif
    !
    return
  end subroutine gridpoint_cluster_middle
  !
end module mod_initgrid
