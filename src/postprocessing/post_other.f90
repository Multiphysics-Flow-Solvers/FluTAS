module mod_post 
  !
  use mod_types
  use mod_common_mpi, only: myid, comm_Cart,ierr
  use mpi
  use mod_bound, only: boundp
!   use mod_param
!   use mod_common_mpi
!   use decomp_2d
!   use mod_phase_indicator
  implicit none
  private
  public  histograms, energyBalance_mf 
contains

subroutine histograms(u,v,w,n,bin_r, minmax,istep, doDeriv)
  use mod_param, only:dli
  implicit none
  real(rp), intent(in), dimension(0:,0:,0:)              :: u,v,w
  integer , intent(in), dimension(3)                     :: n
  integer , intent(in)                                   :: minmax
  integer , intent(in)                                   :: bin_r, istep
  logical , intent(in)                                   :: doDeriv
  real(rp),             dimension(1:n(1),1:n(2),1:n(3))  :: dudx,dvdy,dwdz
  integer                                                :: i,j,k,ip,jp,kp,im,jm,km

  if (doDeriv) then
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ip = i + 1
          jp = j + 1
          kp = k + 1
          im = i - 1
          jm = j - 1
          km = k - 1

          dudx(i,j,k) = ( u(i,j,k) - u(im,j,k) )*dli(1)
          dvdy(i,j,k) = ( v(i,j,k) - v(im,j,k) )*dli(2)
          dwdz(i,j,k) = ( w(i,j,k) - w(im,j,k) )*dli(3)
        enddo
      enddo
    enddo
  endif

call compHist3d(u(1:n(1),1:n(2),1:n(3)),n,bin_r, minmax,istep, 'uvel')
call compHist3d(v(1:n(1),1:n(2),1:n(3)),n,bin_r, minmax,istep, 'vvel')
call compHist3d(w(1:n(1),1:n(2),1:n(3)),n,bin_r, minmax,istep, 'wvel')
call compHist3d(dudx,n,bin_r, minmax,istep, 'dudx')
call compHist3d(dvdy,n,bin_r, minmax,istep, 'dvdy')
call compHist3d(dwdz,n,bin_r, minmax,istep, 'dwdz')
end subroutine histograms


subroutine compHist3d(var,n,bin_r, minmax,istep, varname)
  implicit none
  integer , intent(in), dimension(3)                     :: n
  real(rp), intent(in), dimension(1:,1:,1:)              :: var
  integer , intent(in)                                   :: bin_r, istep
  integer , intent(in)                                   :: minmax
  real(rp)                                               :: minV, maxV
  character(len=4), intent(in)                           :: varname
  integer                                                :: i,j,k, ind, error
  real(rp)                                               :: db, avg
  real(rp),  dimension(0:bin_r-1)                               :: binV
  integer(8) ,  dimension(0:bin_r-1)                               :: histVar, histVar_all
  character(len=7)                                       :: fldnum


  avg = sum(var)/product(n)
  select case(minmax)
  case(0)
    minV = minval(var)
    maxV = maxval(var)
    call mpi_allreduce(MPI_IN_PLACE,minV,1,mpi_real8,mpi_min,comm_cart,error)
    call mpi_allreduce(MPI_IN_PLACE,maxV,1,mpi_real8,mpi_max,comm_cart,error)
  end select
  db   = 0
  db   = (maxV-minV)/(1.0*(bin_r))
  if (db.eq.0) db=1

  do i=0,bin_r-1
    binV(i) = minV + i*db
  enddo
  histVar = 0
  histVar_all = 0
  do i=1,n(1)
    do j=1,n(2)
      do k=1,n(3)
        ind = int((var(i,j,k)-minV)/db)
        histVar(ind) = histVar(ind) +1
      enddo
    enddo
  enddo

  call mpi_allreduce(histVar(0),histVar_all(0),bin_r,mpi_integer8,mpi_sum,comm_cart,error)
  if(myid.eq.0) then
    write(fldnum,'(i7.7)') istep
    open(20, file='data/post/hist/'//varname//'-'//fldnum//'.out')
    write(20,'(A)') '# bin, hist'

    do i=0,bin_r-1
      write(20,'(1E15.7 , 1I8.7)') binV(i), histVar_all(i)
    enddo
      close(20)
  endif
  call mpi_barrier(comm_cart,error)
  return
end subroutine compHist3d




subroutine energyBalance_mf(u,v,w,p,psi,kappa,rho,mu,n,time,nh_d,nh_vof,halo_vof, dzc, dzf)
  use mod_param, only: gacc, dli,sigma, cbcvof, bcvof, is_outflow,dl, dims, rho1, rho2
  implicit none
  real(rp), intent(in),   dimension(0:,0:,0:)                     :: u,v,w, rho, mu,psi, kappa,p
  integer,  intent(in),   dimension(3)                            :: n
  real(rp), intent(in)                                            :: time
  integer,  intent(in)                                            :: nh_d,nh_vof
  integer,  intent(in),   dimension(2)                            :: halo_vof
  real(rp)         , intent(in   ), dimension(-2:)                :: dzc,dzf
  integer                                                         :: i,j,k, ip,im,jp,jm,kp,km,error
  real(rp),            dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1)      :: up,vp,wp
  real(rp)                                                        :: Uavg, Vavg, Wavg
  real(rp),            dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1,6)    :: S      
  real(rp)                                                        :: VGRux, VGRuy, VGRuz, &
                                                                     VGRvx, VGRvy, VGRvz, &
                                                                     VGRwx, VGRwy, VGRwz, &
                                                                     duTxx, duTxy, duTxz, &
                                                                     duTyx, duTyy, duTyz, &
                                                                     duTzx, duTzy, duTzz, &
                                                                     dupx , dupy,  dupz 
  ! real(rp),               dimension(6)                            :: S
  real(rp)                                                        :: S_S, E_E, &
                                                                     ens_l, eps_l, KEp_l, KE_l,  &
                                                                     ens_g, eps_g, KEp_g, KE_g,  &
                                                                     Tp_l,  Tp_g,  Tnu_l, Tnu_g, &
                                                                     Psi_nu, volF, urms, vrms, wrms


  E_E     = 0.0
  ens_l   = 0.0
  ens_g   = 0.0
  eps_l   = 0.0
  eps_g   = 0.0
  KEp_l   = 0.0
  KEp_g   = 0.0
  KE_l    = 0.0
  KE_g    = 0.0
  Tnu_l   = 0.0
  Tnu_g   = 0.0
  Tp_l    = 0.0
  Tp_g    = 0.0
  Psi_nu  = 0.0
  volF    = 0.0
  Uavg = sum(u)/product(n)
  Vavg = sum(v)/product(n)
  Wavg = sum(w)/product(n)
  call mpi_allreduce(MPI_IN_PLACE,Uavg,1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,Vavg,1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,Wavg,1,mpi_real8,mpi_sum,comm_cart,error)
  Uavg = Uavg/product(dims)
  Vavg = Vavg/product(dims)
  Wavg = Wavg/product(dims)

  up = u-Uavg
  vp = v-Uavg
  wp = w-Uavg
  urms = sum(up**2)/product(n)
  vrms = sum(vp**2)/product(n)
  wrms = sum(wp**2)/product(n)
  call mpi_allreduce(MPI_IN_PLACE,urms,1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,vrms,1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,wrms,1,mpi_real8,mpi_sum,comm_cart,error)
  urms = urms/product(dims)
  vrms = vrms/product(dims)
  wrms = wrms/product(dims)
  ! call bounduvw_b(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,up,vp,wp)

  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        ip = i + 1
        jp = j + 1
        kp = k + 1
        im = i - 1
        jm = j - 1
        km = k - 1

        VGRux = ( up(i,j,k) - up(im,j,k) )*dli(1)

        VGRuy = 0.25*( up(i,j,k) + up(i,jp,k) + up(im,j,k) + up(im,jp,k) )*dli(2) - &
                0.25*( up(i,j,k) + up(i,jm,k) + up(im,j,k) + up(im,jm,k) )*dli(2)  

        VGRuz = 0.25*( up(i,j,k) + up(i,j,kp) + up(im,j,k) + up(im,j,kp) )*dli(3) - &
                0.25*( up(i,j,k) + up(i,j,km) + up(im,j,k) + up(im,j,km) )*dli(3)

        VGRvx = 0.25*( vp(i,j,k) + vp(ip,j,k) + vp(i,jm,k) + vp(ip,jm,k) )*dli(1) - &
                0.25*( vp(i,j,k) + vp(im,j,k) + vp(i,jm,k) + vp(im,jm,k) )*dli(1)

        VGRvy = ( vp(i,j,k) - vp(i,jm,k) )*dli(2)

        VGRvz = 0.25*( vp(i,j,k) + vp(i,j,kp) + vp(i,jm,k) + vp(i,jm,kp) )*dli(3) - &
                0.25*( vp(i,j,k) + vp(i,j,km) + vp(i,jm,k) + vp(i,jm,km) )*dli(3)

        VGRwx = 0.25*( wp(i,j,k) + wp(ip,j,k) + wp(i,j,km) + wp(ip,j,km) )*dli(1) - &
                0.25*( wp(i,j,k) + wp(im,j,k) + wp(i,j,km) + wp(im,j,km) )*dli(1)

        VGRwy = 0.25*( wp(i,j,k) + wp(i,jp,k) + wp(i,j,km) + wp(i,jp,km) )*dli(2) - &
                0.25*( wp(i,j,k) + wp(i,jm,k) + wp(i,j,km) + wp(i,jm,km) )*dli(2)

        VGRwz = ( wp(i,j,k) - wp(i,j,km) )*dli(3)

        !-----S11 = dU/dx
        S(i,j,k,1) = VGRux
        !-----S12 = 0.5*(dU/dy+dV/dx) = S21
        S(i,j,k,2) = 0.5 * (VGRuy + VGRvx)
        !-----S13 = 0.5*(dU/dz+dW/dx) = S31
        S(i,j,k,3) = 0.5 * (VGRuz + VGRwx)
        !-----S22 = dV/dy
        S(i,j,k,4) = VGRvy
        !-----S23 = 0.5*(dV/dz+dW/dy) = S32
        S(i,j,k,5) = 0.5 * (VGRvz + VGRwy)
        !-----S33 = dW/dz
        S(i,j,k,6) = VGRwz

        E_E = (VGRwy-VGRvz)**2. + (VGRuz-VGRwx)**2. + (VGRvx-VGRuy)**2. 

        ens_l = ens_l + psi(i,j,k)*E_E 
        ens_g = ens_g + (1-psi(i,j,k))*E_E
      enddo
    enddo
  enddo

  ! WARNING: valid only for periodic BC
  !do i=1,6
    !call boundp(cbcvof,n,bcvof,dl,dzc,dzf,S(:,:,:,i))
    !call boundp(cbcvof,n,bcvof,nh_d,nh_vof,halo_vof,dl,dzc,dzf,S(:,:,:,i))
  !enddo

  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
          ip = i + 1
          jp = j + 1
          kp = k + 1
          im = i - 1
          jm = j - 1
          km = k - 1

          duTxx =  0.5* (up(i ,j,k)*(S(i,j,k,1)+S(ip,j,k,1))*(mu(i,j,k)+mu(ip,j,k)) - &
                         up(im,j,k)*(S(i,j,k,1)+S(im,j,k,1))*(mu(i,j,k)+mu(im,j,k)))*dli(1)

          duTxy = 0.125*(up(i ,j,k) + up(i ,jp,k) + up(im,j,k) + up(im,jp,k))*       &
                        (S(i,jp,k,2)+ S(i,j,k,2)) * (mu(i,jp,k)+mu(i,j,k))*dli(2) -    &
                  0.125*(up(i ,j,k) + up(i ,jm,k) + up(im,j,k)+ up(im,jm,k))*           &
                        (S(i,jm,k,2)+ S(i,j,k,2)) * (mu(i ,j,k)+mu(i ,jm,k))*dli(2)  

          duTxz = 0.125*(up(i ,j,k) + up(i ,j,kp) + up(im,j,k) + up(im,j,kp))*        &
                        (S(i,j,kp,3)+ S(i,j,k,3)) * (mu(i,j,kp)+mu(i ,j,k))*dli(3) -  &
                  0.125*(up(i ,j,k) + up(i ,j,km) + up(im,j,k) + up(im,j,km))*         &
                        (S(i,j,k,3)+ S(i,j,km,3)) * (mu(i ,j,k)+mu(i ,j,km))*dli(3)

          duTyx = 0.125*( vp(i,j ,k) + vp(ip,j ,k) + vp(i,jm,k)  + vp(ip,jm,k))*         &
                        (S(ip,j,k,2) + S(i,j,k,2)) * (mu(ip,j ,k)+mu(i,j ,k))*dli(1) - &
                  0.125*( vp(i,j ,k) + vp(im,j ,k) + vp(i,jm,k)  + vp(im,jm,k))*          &
                        ((S(i,j,k,2)+  S(im,j,k,2)) * (mu(i,j ,k) + mu(im,j ,k)))*dli(1)

          duTyy =   0.5*( vp(i,j ,k)*(S(i,j,k,4)+S(i,jp,k,4))*(mu(i,j,k)+mu(i,jp,k)) -   &
                        vp(i,jm,k)*(S(i,j,k,4)+S(i,jm,k,4))*(mu(i,j,k)+mu(i,jm,k)))*dli(2)

          duTyz = 0.125*( vp(i,j ,k) + vp(i,j,kp)  + vp(i,jm,k) + vp(i,jm,kp) )*          &
                        (S(i,j ,k,5) + S(i,j,kp,5)) * (mu(i,j ,k)+mu(i,j ,kp))*dli(3)-   &
                  0.125*( vp(i,j ,k) + vp(i,j ,km) + vp(i,jm,k) + vp(i,jm,km) )*           &
                        (S(i,j ,k,5) + S(i,j ,km,5))*(mu(i,j ,k) + mu(i,j ,km))*dli(3)

          duTzx = 0.125*( wp(i,j,k ) + wp(ip,j,k )  + wp(i,j,km) + wp(ip,j,km))*           &
                        (S(ip,j,k,3) + S(i,j,k ,3)) * (mu(ip,j,k) + mu(i,j,k ))*dli(1) -   &
                  0.125*( wp(i,j,k ) + wp(im,j,k ) + wp(i,j,km)  + wp(im,j,km))*            & 
                        (S(i,j,k ,3) + S(im,j,k,3))*(mu(i,j,k )  + mu(im,j,k ))*dli(1)

          duTzy = 0.125*( wp(i,j,k )  + wp(i,jp,k ) + wp(i,j,km) + wp(i,jp,km))*          &
                        (S(i,j,k,5)  + S(i,jp,k,5)) * (mu(i,j,k)+mu(i,jp,k ))* dli(2) -    &
                  0.125*( wp(i,j,k ) + wp(i,jm,k ) + wp(i,j,km) + wp(i,jm,km))*           &    
                        (S(i,j,k ,5) + S(i,jm,k ,5)) * (mu(i,j,k )+mu(i,jm,k ))*dli(2)

          duTzz =   0.5*( wp(i,j,k )*(S(i,j,kp,6)+S(i,j ,k,6))*(mu(i,j,k)+mu(i,j,kp))     - &
                        wp(i,j,km)*(S(i,j,k ,6)+S(i,j,km,6))*(mu(i,j,k)+mu(i,j,km)))*dli(3)

          dupx = 0.5*( up(i ,j,k)*(p(ip,j,k)+p(i,j,k)) - &
                       up(im,j,k)*(p(im,j,k)+p(i,j,k)) )*dli(1)
          dupy = 0.5*( vp(i,j ,k)*(p(i,jp,k)+p(i,j,k)) - &
                       vp(i,jm,k)*(p(i,jm,k)+p(i,j,k)) )*dli(2)
          dupz = 0.5*( wp(i,j,k )*(p(i,j,kp)+p(i,j,k)) -  &
                       wp(i,j,km)*(p(i,j,km)+p(i,j,k)) )*dli(3)

          



          S_S = S(i,j,k,1)*S(i,j,k,1) + & 
             2.*S(i,j,k,2)*S(i,j,k,2) + &
             2.*S(i,j,k,3)*S(i,j,k,3) + &
                S(i,j,k,4)*S(i,j,k,4) + &
             2.*S(i,j,k,5)*S(i,j,k,5) + &
                S(i,j,k,6)*S(i,j,k,6)


       
          eps_l = eps_l + psi(i,j,k)    *2. * mu(i,j,k) * S_S 
          eps_g = eps_g + (1-psi(i,j,k))*2. * mu(i,j,k) * S_S 



          KEp_g = KEp_g + rho2*&
          (1.0-psi(i,j,k))*0.5*(0.25*(up(i,j,k)+up(im,j,k))**2. + &
                                0.25*(vp(i,j,k)+vp(i,jm,k))**2. + &
                                0.25*(wp(i,j,k)+wp(i,j,km))**2.)
          KEp_l = KEp_l + rho1*&
                psi(i,j,k)*0.5*(0.25*(up(i,j,k)+up(im,j,k))**2. + &
                                0.25*(vp(i,j,k)+vp(i,jm,k))**2. + &
                                0.25*(wp(i,j,k)+wp(i,j,km))**2.)

          KE_g = KE_g + rho2*&
          (1.0-psi(i,j,k))*0.5*(0.25*(u(i,j,k)+u(im,j,k))**2. + &
                                0.25*(v(i,j,k)+v(i,jm,k))**2. + &
                                0.25*(w(i,j,k)+w(i,j,km))**2.)
          KE_l = KE_l + rho1*&
               psi(i,j,k)*0.5*(0.25*(u(i,j,k)+u(im,j,k))**2. + &
                               0.25*(v(i,j,k)+v(i,jm,k))**2. + &
                               0.25*(w(i,j,k)+w(i,j,km))**2.)

          Tnu_l = Tnu_l+    psi(i,j,k) *(duTxx+duTxy+duTxz+duTyx+duTyy+duTyz+duTzx+duTzy+duTzz)
          Tnu_g = Tnu_g+ (1-psi(i,j,k))*(duTxx+duTxy+duTxz+duTyx+duTyy+duTyz+duTzx+duTzy+duTzz)

          Tp_l = Tp_l+    psi(i,j,k)* (dupx+dupy+dupz)
          Tp_g = Tp_g+ (1-psi(i,j,k))*(dupx+dupy+dupz)

          Psi_nu = Psi_nu + sigma*kappa(i,j,k)*( &
                                  dli(1)*(0.5d0*(psi(ip,j,k)-psi(im,j,k)))*0.5d0*(up(i,j,k)+up(im,j,k)) + &
                                  dli(2)*(0.5d0*(psi(i,jp,k)-psi(i,jm,k)))*0.5d0*(vp(i,j,k)+vp(i,jm,k)) + &
                                  dli(3)*(0.5d0*(psi(i,j,kp)-psi(i,j,km)))*0.5d0*(wp(i,j,k)+wp(i,j,km))) 
          volF = volF + psi(i,j,k)

      enddo
    enddo
  enddo

  call mpi_allreduce(MPI_IN_PLACE,E_E,   1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,ens_l, 1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,ens_g, 1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,eps_l, 1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,eps_g, 1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,KEp_l, 1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,KEp_g, 1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,KE_l,  1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,KE_g,  1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,Tnu_l, 1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,Tnu_g, 1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,Tp_l,  1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,Tp_g,  1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,Psi_nu,1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,volF  ,1,mpi_real8,mpi_sum,comm_cart,error)

  E_E     = E_E   /product(n)/product(dims) 
  ens_l   = ens_l /product(n)/product(dims) 
  ens_g   = ens_g /product(n)/product(dims) 
  eps_l   = eps_l /product(n)/product(dims) 
  eps_g   = eps_g /product(n)/product(dims) 
  KEp_l   = KEp_l /product(n)/product(dims) 
  KEp_g   = KEp_g /product(n)/product(dims) 
  KE_l    = KE_l  /product(n)/product(dims) 
  KE_g    = KE_g  /product(n)/product(dims) 
  Tnu_l   = Tnu_l /product(n)/product(dims) 
  Tnu_g   = Tnu_g /product(n)/product(dims) 
  Tp_l    = Tp_l  /product(n)/product(dims) 
  Tp_g    = Tp_g  /product(n)/product(dims) 
  Psi_nu  = Psi_nu/product(n)/product(dims) 
  volF    = volF  /product(n)/product(dims) 
 
    if(myid .eq. 0 ) then
    open(92,file='data/post/balance/balance.out',position='append')
    write(92,'(16E15.7)') time,E_E, ens_l, ens_g, eps_l, eps_g, KE_l, KE_g, &
                                  KEp_l, KEp_g, Tnu_l, Tnu_g, Tp_l, Tp_g, Psi_nu, volF
    close(92)
  endif

  if(myid .eq. 0 ) then
     open(92,file='data/post/balance/rms.out',position='append')
    write(92,'(3E15.7)') urms,vrms,wrms
    close(92)
  endif
  return
end subroutine energyBalance_mf
!
end module mod_post
