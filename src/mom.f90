!
! SPDX-License-Identifier: MIT
!
module mod_mom
  !
  use mod_types, only: rp
  !@cuf use cudafor
  !
  implicit none
  !
  private
#if defined(_USE_VOF)
  public  :: momad_xyz_tw_cen,momad_xyz_tw_fll
#else
  public  :: momad_xyz_sp
#endif
  !
  contains
  !
#if defined(_USE_VOF)
  subroutine momad_xyz_tw_cen(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,mu,rho, &
                              dudt,dvdt,dwdt)  
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in ), dimension(1-nh_d:                ) :: dzci,dzfi
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in ), dimension(     0:,     0:,     0:) :: mu,rho
    real(rp), intent(out), dimension(     1:,     1:,     1:) :: dudt,dvdt,dwdt
    !
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real(rp) :: muxm,muxp,muym,muyp,muzm,muzp
    real(rp) :: rhox,rhoy,rhoz
    integer  :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm, &
                uvip,uvim,vvjp,vvjm,wvkp,wvkm, &
                uwip,uwim,vwjp,vwjm,wwkp,wwkm
    !
    !@cuf attributes(managed) :: u, v, w, mu, rho, dzci, dzfi, dudt, dvdt, dwdt
    !    
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i + 1
          im = i - 1
          jp = j + 1
          jm = j - 1
          kp = k + 1
          km = k - 1
          !
#if defined(_TWOD)
          dudt(i,j,k) = 0._rp
#else
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          ! 1. r.h.s. of x-momentum
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          uuip  = 0.25_rp*( u(ip,j,k)+u(i,j,k) )*( u(ip,j ,k )+u(i,j ,k ) )
          uuim  = 0.25_rp*( u(im,j,k)+u(i,j,k) )*( u(im,j ,k )+u(i,j ,k ) )
          uvjp  = 0.25_rp*( u(i,jp,k)+u(i,j,k) )*( v(ip,j ,k )+v(i,j ,k ) )
          uvjm  = 0.25_rp*( u(i,jm,k)+u(i,j,k) )*( v(ip,jm,k )+v(i,jm,k ) )
          uwkp  = 0.25_rp*( u(i,j,kp)+u(i,j,k) )*( w(ip,j ,k )+w(i,j ,k ) )
          uwkm  = 0.25_rp*( u(i,j,km)+u(i,j,k) )*( w(ip,j ,km)+w(i,j ,km) )
          !
          dudxp = (u(ip,j ,k)-u(i ,j ,k))*dxi
          dudxm = (u(i ,j ,k)-u(im,j ,k))*dxi
          dvdxp = (v(ip,j ,k)-v(i ,j ,k))*dxi
          dvdxm = (v(ip,jm,k)-v(i ,jm,k))*dxi
          dudyp = (u(i ,jp,k)-u(i ,j ,k))*dyi
          dudym = (u(i ,j ,k)-u(i ,jm,k))*dyi
          dudzp = (u(i ,j,kp)-u(i ,j,k ))*dzci(k)
          dudzm = (u(i ,j,k )-u(i ,j,km))*dzci(km)
          dwdxp = (w(ip,j,k )-w(i ,j,k ))*dxi
          dwdxm = (w(ip,j,km)-w(i ,j,km))*dxi
          !
          muxp = mu(ip,j,k)
          muxm = mu(i ,j,k)
          muyp = 0.25_rp*(mu(i,j,k)+mu(i,jp,k)+mu(ip,jp,k)+mu(ip,j,k))
          muym = 0.25_rp*(mu(i,j,k)+mu(i,jm,k)+mu(ip,jm,k)+mu(ip,j,k))
          muzp = 0.25_rp*(mu(i,j,k)+mu(i,j,kp)+mu(ip,j,kp)+mu(ip,j,k))
          muzm = 0.25_rp*(mu(i,j,k)+mu(i,j,km)+mu(ip,j,km)+mu(ip,j,k))
          !
          rhox = 0.5_rp*(rho(  ip,j,k)+rho(  i,j,k))
          !
          dudt(i,j,k) = &
                        dxi*(     -uuip + uuim ) + &
                        dyi*(     -uvjp + uvjm ) + &
                        dzfi(k)*( -uwkp + uwkm ) + &
                        ( &
                        dxi*(    (dudxp+dudxp)*muxp-(dudxm+dudxm)*muxm) + &
                        dyi*(    (dudyp+dvdxp)*muyp-(dudym+dvdxm)*muym) + &
                        dzfi(k)*((dudzp+dwdxp)*muzp-(dudzm+dwdxm)*muzm) ) /rhox
#endif 
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          ! 2. r.h.s. of y-momentum
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          uvip  = 0.25_rp*( u(i ,j,k)+u(i ,jp,k) )*( v(i,j,k )+v(ip,j ,k) )
          uvim  = 0.25_rp*( u(im,j,k)+u(im,jp,k) )*( v(i,j,k )+v(im,j ,k) )
          vvjp  = 0.25_rp*( v(i,j,k )+v(i,jp,k)  )*( v(i,j,k )+v(i ,jp,k) )
          vvjm  = 0.25_rp*( v(i,j,k )+v(i,jm,k)  )*( v(i,j,k )+v(i ,jm,k) )
          wvkp  = 0.25_rp*( w(i,j,k )+w(i,jp,k)  )*( v(i,j,kp)+v(i ,j ,k) )
          wvkm  = 0.25_rp*( w(i,j,km)+w(i,jp,km) )*( v(i,j,km)+v(i ,j ,k) )
          !
          dvdxp = (v(ip,j ,k)-v(i ,j ,k))*dxi
          dvdxm = (v(i ,j ,k)-v(im,j ,k))*dxi
          dudyp = (u(i ,jp,k)-u(i ,j ,k))*dyi
          dudym = (u(im,jp,k)-u(im,j ,k))*dyi
          dvdyp = (v(i,jp ,k)-v(i ,j ,k))*dyi
          dvdym = (v(i,j  ,k)-v(i ,jm,k))*dyi
          dvdzp = (v(i,j,kp) -v(i,j,k  ))*dzci(k)
          dvdzm = (v(i,j,k ) -v(i,j,km ))*dzci(km)
          dwdyp = (w(i,jp,k )-w(i,j,k  ))*dyi
          dwdym = (w(i,jp,km)-w(i,j,km ))*dyi
          !
          muxp = 0.25_rp*(mu(i,j,k)+mu(ip,j,k)+mu(ip,jp,k)+mu(i,jp,k))
          muxm = 0.25_rp*(mu(i,j,k)+mu(im,j,k)+mu(im,jp,k)+mu(i,jp,k))
          muyp = mu(i,jp,k)
          muym = mu(i,j,k)
          muzp = 0.25_rp*(mu(i,j,k)+mu(i,jp,k)+mu(i,jp,kp)+mu(i,j,kp))
          muzm = 0.25_rp*(mu(i,j,k)+mu(i,jp,k)+mu(i,jp,km)+mu(i,j,km))
          !
          rhoy = 0.5_rp*(rho(  i,jp,k )+rho(i,j,k ))
          !
          dvdt(i,j,k) = &
                        dxi*(     -uvip + uvim ) + &
                        dyi*(     -vvjp + vvjm ) + & 
                        dzfi(k)*( -wvkp + wvkm ) + & 
                        ( &
                        dxi*(    (dvdxp+dudyp)*muxp-(dvdxm+dudym)*muxm) + &
                        dyi*(    (dvdyp+dvdyp)*muyp-(dvdym+dvdym)*muym) + &
                        dzfi(k)*((dvdzp+dwdyp)*muzp-(dvdzm+dwdym)*muzm) ) /rhoy 
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          ! 3. r.h.s. of z-momentum
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          uwip  = 0.25_rp*( w(i,j,k)+w(ip,j,k) )*( u(i ,j ,k)+u(i ,j ,kp) )
          uwim  = 0.25_rp*( w(i,j,k)+w(im,j,k) )*( u(im,j ,k)+u(im,j ,kp) )
          vwjp  = 0.25_rp*( w(i,j,k)+w(i,jp,k) )*( v(i ,j ,k)+v(i ,j ,kp) )
          vwjm  = 0.25_rp*( w(i,j,k)+w(i,jm,k) )*( v(i ,jm,k)+v(i ,jm,kp) )
          wwkp  = 0.25_rp*( w(i,j,k)+w(i,j,kp) )*( w(i ,j ,k)+w(i ,j ,kp) )
          wwkm  = 0.25_rp*( w(i,j,k)+w(i,j,km) )*( w(i ,j ,k)+w(i ,j ,km) )
          !
          dwdxp = (w(ip,j,k )-w(i ,j,k))*dxi
          dwdxm = (w(i ,j,k )-w(im,j,k))*dxi
          dudzp = (u(i ,j,kp)-u(i ,j,k))*dzci(k)
          dudzm = (u(im,j,kp)-u(im,j,k))*dzci(k)
          dwdyp = (w(i,jp,k )-w(i,j ,k))*dyi
          dwdym = (w(i,j ,k )-w(i,jm,k))*dyi
          dvdzp = (v(i,j ,kp)-v(i,j ,k))*dzci(k)
          dvdzm = (v(i,jm,kp)-v(i,jm,k))*dzci(k)
          dwdzp = (w(i,j,kp )-w(i,j,k ))*dzfi(kp)
          dwdzm = (w(i,j,k  )-w(i,j,km))*dzfi(k)
          !
          muxp = 0.25_rp*(mu(i,j,k)+mu(i,j,kp)+mu(ip,j ,kp)+mu(ip,j ,k) )
          muxm = 0.25_rp*(mu(i,j,k)+mu(i,j,kp)+mu(im,j ,kp)+mu(im,j ,k) )
          muyp = 0.25_rp*(mu(i,j,k)+mu(i,j,kp)+mu(i ,jp,kp)+mu(i ,jp,k) )
          muym = 0.25_rp*(mu(i,j,k)+mu(i,j,kp)+mu(i ,jm,kp)+mu(i ,jm,k) )
          muzp = mu(i,j,kp)
          muzm = mu(i,j,k )
          !
          rhoz   = 0.5_rp*(rho(  i,j,kp)+rho(  i,j,k))
          !
          dwdt(i,j,k) = &
                        dxi*(     -uwip + uwim ) + &
                        dyi*(     -vwjp + vwjm ) + &
                        dzci(k)*( -wwkp + wwkm ) + &
                        ( &
                        dxi*(    (dwdxp+dudzp)*muxp-(dwdxm+dudzm)*muxm ) + &
                        dyi*(    (dwdyp+dvdzp)*muyp-(dwdym+dvdzm)*muym ) + &
                        dzci(k)*((dwdzp+dwdzp)*muzp-(dwdzm+dwdzm)*muzm ) )/rhoz  

          !
        enddo
      enddo
    enddo
    !$acc end kernels 
    !$OMP END PARALLEL DO
    !
    return
  end subroutine momad_xyz_tw_cen
  !
  subroutine momad_xyz_tw_fll(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,mu,rho, &
                              dudt,dvdt,dwdt)  
    !
    use mod_funcs, only: interp_g
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in ), dimension(1-nh_d:                ) :: dzci,dzfi
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in ), dimension(     0:,     0:,     0:) :: mu,rho
    real(rp), intent(out), dimension(     1:,     1:,     1:) :: dudt,dvdt,dwdt
    !
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real(rp) :: muxm,muxp,muym,muyp,muzm,muzp
    real(rp) :: rhox,rhoy,rhoz
    integer  :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: vec1,vec2,vec3,vec4,vec5
    real(rp) :: uc,vc,wc
    real(rp) :: ug_u,vg_u,wg_u, &
                ug_v,vg_v,wg_v, &
                ug_w,vg_w,wg_w
    integer  :: q
    !
    !@cuf attributes(managed) :: u, v, w, mu, rho, dzci, dzfi, dudt, dvdt, dwdt
    !    
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i + 1
          im = i - 1
          jp = j + 1
          jm = j - 1
          kp = k + 1
          km = k - 1
          !
#if defined(_TWOD)
          dudt(i,j,k) = 0._rp
#else
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          ! 1. r.h.s. of x-momentum
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          vec1 = u(i-2,j,k)
          vec2 = u(i-1,j,k)
          vec3 = u(i  ,j,k)
          vec4 = u(i+1,j,k)
          vec5 = u(i+2,j,k)
          uc   = u(i,j,k)
          ug_u = (uc)*(vec1-8._rp*vec2+8._rp*vec4-vec5)*(1._rp/12._rp)*dxi + & 
                 3._rp*abs(uc)*(vec1-4._rp*vec2+6._rp*vec3-4._rp*vec4+vec5)*(1._rp/12._rp)*dxi
          !
          vec1 = u(i,j-2,k)
          vec2 = u(i,j-1,k)
          vec3 = u(i,j  ,k)
          vec4 = u(i,j+1,k)
          vec5 = u(i,j+2,k)
          vc   = 0.25_rp*( v(i,j,k)+v(ip,j,k)+v(i,jm,k)+v(ip,jm,k) )
          vg_u = (vc)*(vec1-8._rp*vec2+8._rp*vec4-vec5)*(1._rp/12._rp)*dyi + & 
                 3._rp*abs(vc)*(vec1-4._rp*vec2+6._rp*vec3-4._rp*vec4+vec5)*(1._rp/12._rp)*dyi
          !
          vec1 = u(i,j,k-2)
          vec2 = u(i,j,k-1)
          vec3 = u(i,j,k  )
          vec4 = u(i,j,k+1)
          vec5 = u(i,j,k+2)
          wc   = 0.25_rp*( w(i,j,k)+w(ip,j,k)+w(i,j,km)+w(ip,j,km) )
          wg_u = (wc)*(vec1-8._rp*vec2+8._rp*vec4-vec5)*(1._rp/12._rp)*dzi + & 
                 3._rp*abs(wc)*(vec1-4._rp*vec2+6._rp*vec3-4._rp*vec4+vec5)*(1._rp/12._rp)*dzi
          !
          dudxp = (u(ip,j ,k)-u(i ,j ,k))*dxi
          dudxm = (u(i ,j ,k)-u(im,j ,k))*dxi
          dvdxp = (v(ip,j ,k)-v(i ,j ,k))*dxi
          dvdxm = (v(ip,jm,k)-v(i ,jm,k))*dxi
          dudyp = (u(i ,jp,k)-u(i ,j ,k))*dyi
          dudym = (u(i ,j ,k)-u(i ,jm,k))*dyi
          dudzp = (u(i ,j,kp)-u(i ,j,k ))*dzci(k)
          dudzm = (u(i ,j,k )-u(i ,j,km))*dzci(km)
          dwdxp = (w(ip,j,k )-w(i ,j,k ))*dxi
          dwdxm = (w(ip,j,km)-w(i ,j,km))*dxi
          !
          muxp = mu(ip,j,k)
          muxm = mu(i ,j,k)
          muyp = 0.25_rp*(mu(i,j,k)+mu(i,jp,k)+mu(ip,jp,k)+mu(ip,j,k))
          muym = 0.25_rp*(mu(i,j,k)+mu(i,jm,k)+mu(ip,jm,k)+mu(ip,j,k))
          muzp = 0.25_rp*(mu(i,j,k)+mu(i,j,kp)+mu(ip,j,kp)+mu(ip,j,k))
          muzm = 0.25_rp*(mu(i,j,k)+mu(i,j,km)+mu(ip,j,km)+mu(ip,j,k))
          !
          rhox = 0.5_rp*(rho(  ip,j,k)+rho(  i,j,k))
          !
          dudt(i,j,k) = &
                        - (ug_u + vg_u + wg_u) + &
                        ( &
                        dxi*(    (dudxp+dudxp)*muxp-(dudxm+dudxm)*muxm) + &
                        dyi*(    (dudyp+dvdxp)*muyp-(dudym+dvdxm)*muym) + &
                        dzfi(k)*((dudzp+dwdxp)*muzp-(dudzm+dwdxm)*muzm) ) /rhox 
#endif
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          ! 2. r.h.s. of y-momentum
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          vec1 = v(i-2,j,k)
          vec2 = v(i-1,j,k)
          vec3 = v(i  ,j,k)
          vec4 = v(i+1,j,k)
          vec5 = v(i+2,j,k)
          uc   = 0.25_rp*( u(i,j,k)+u(i,jp,k)+u(im,j,k)+u(im,jp,k) )
          ug_v = (uc)*(vec1-8._rp*vec2+8._rp*vec4-vec5)*(1._rp/12._rp)*dxi + & 
                 3._rp*abs(uc)*(vec1-4._rp*vec2+6._rp*vec3-4._rp*vec4+vec5)*(1._rp/12._rp)*dxi
          !
          vec1 = v(i,j-2,k)
          vec2 = v(i,j-1,k)
          vec3 = v(i,j  ,k)
          vec4 = v(i,j+1,k)
          vec5 = v(i,j+2,k)
          vc   = v(i,j,k) 
          vg_v = (vc)*(vec1-8._rp*vec2+8._rp*vec4-vec5)*(1._rp/12._rp)*dyi + & 
                 3._rp*abs(vc)*(vec1-4._rp*vec2+6._rp*vec3-4._rp*vec4+vec5)*(1._rp/12._rp)*dyi
          !
          vec1 = v(i,j,k-2)
          vec2 = v(i,j,k-1)
          vec3 = v(i,j,k  )
          vec4 = v(i,j,k+1)
          vec5 = v(i,j,k+2)
          wc   = 0.25_rp*( w(i,j,k)+w(i,jp,k)+w(i,j,km)+w(i,jp,km) )
          wg_v = (wc)*(vec1-8._rp*vec2+8._rp*vec4-vec5)*(1._rp/12._rp)*dzi + & 
                 3._rp*abs(wc)*(vec1-4._rp*vec2+6._rp*vec3-4._rp*vec4+vec5)*(1._rp/12._rp)*dzi
          !
          dvdxp = (v(ip,j ,k)-v(i ,j ,k))*dxi
          dvdxm = (v(i ,j ,k)-v(im,j ,k))*dxi
          dudyp = (u(i ,jp,k)-u(i ,j ,k))*dyi
          dudym = (u(im,jp,k)-u(im,j ,k))*dyi
          dvdyp = (v(i,jp ,k)-v(i ,j ,k))*dyi
          dvdym = (v(i,j  ,k)-v(i ,jm,k))*dyi
          dvdzp = (v(i,j,kp) -v(i,j,k  ))*dzci(k)
          dvdzm = (v(i,j,k ) -v(i,j,km ))*dzci(km)
          dwdyp = (w(i,jp,k )-w(i,j,k  ))*dyi
          dwdym = (w(i,jp,km)-w(i,j,km ))*dyi
          !
          muxp = 0.25_rp*(mu(i,j,k)+mu(ip,j,k)+mu(ip,jp,k)+mu(i,jp,k))
          muxm = 0.25_rp*(mu(i,j,k)+mu(im,j,k)+mu(im,jp,k)+mu(i,jp,k))
          muyp = mu(i,jp,k)
          muym = mu(i,j,k)
          muzp = 0.25_rp*(mu(i,j,k)+mu(i,jp,k)+mu(i,jp,kp)+mu(i,j,kp))
          muzm = 0.25_rp*(mu(i,j,k)+mu(i,jp,k)+mu(i,jp,km)+mu(i,j,km))
          !
          rhoy = 0.5_rp*(rho(  i,jp,k )+rho(i,j,k ))
          !
          dvdt(i,j,k) = &
                        - ( ug_v + vg_v + wg_v) + &
                        ( &
                        dxi*(    (dvdxp+dudyp)*muxp-(dvdxm+dudym)*muxm) + &
                        dyi*(    (dvdyp+dvdyp)*muyp-(dvdym+dvdym)*muym) + &
                        dzfi(k)*((dvdzp+dwdyp)*muzp-(dvdzm+dwdym)*muzm) ) /rhoy 
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          ! 3. r.h.s. of z-momentum
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          vec1 = w(i-2,j,k)
          vec2 = w(i-1,j,k)
          vec3 = w(i  ,j,k)
          vec4 = w(i+1,j,k)
          vec5 = w(i+2,j,k)
          uc   = 0.25_rp*( u(i,j,k)+u(i,j,kp)+u(im,j,k)+u(im,j,kp) )
          ug_w = (uc)*(vec1-8._rp*vec2+8._rp*vec4-vec5)*(1._rp/12._rp)*dxi + & 
                 3._rp*abs(uc)*(vec1-4._rp*vec2+6._rp*vec3-4._rp*vec4+vec5)*(1._rp/12._rp)*dxi
          !
          vec1 = w(i,j-2,k)
          vec2 = w(i,j-1,k)
          vec3 = w(i,j  ,k)
          vec4 = w(i,j+1,k)
          vec5 = w(i,j+2,k)
          vc   = 0.25_rp*( v(i,j,k)+v(i,j,kp)+v(i,jm,k)+v(i,jm,kp) )
          vg_w = (vc)*(vec1-8._rp*vec2+8._rp*vec4-vec5)*(1._rp/12._rp)*dyi + & 
                 3._rp*abs(vc)*(vec1-4._rp*vec2+6._rp*vec3-4._rp*vec4+vec5)*(1._rp/12._rp)*dyi
          !
          vec1 = w(i,j,k-2)
          vec2 = w(i,j,k-1)
          vec3 = w(i,j,k  )
          vec4 = w(i,j,k+1)
          vec5 = w(i,j,k+2)
          wc   = w(i,j,k)
          wg_w = (wc)*(vec1-8._rp*vec2+8._rp*vec4-vec5)*(1._rp/12._rp)*dzi + & 
                 3._rp*abs(wc)*(vec1-4._rp*vec2+6._rp*vec3-4._rp*vec4+vec5)*(1._rp/12._rp)*dzi
          !
          dwdxp = (w(ip,j,k )-w(i ,j,k))*dxi
          dwdxm = (w(i ,j,k )-w(im,j,k))*dxi
          dudzp = (u(i ,j,kp)-u(i ,j,k))*dzci(k)
          dudzm = (u(im,j,kp)-u(im,j,k))*dzci(k)
          dwdyp = (w(i,jp,k )-w(i,j ,k))*dyi
          dwdym = (w(i,j ,k )-w(i,jm,k))*dyi
          dvdzp = (v(i,j ,kp)-v(i,j ,k))*dzci(k)
          dvdzm = (v(i,jm,kp)-v(i,jm,k))*dzci(k)
          dwdzp = (w(i,j,kp )-w(i,j,k ))*dzfi(kp)
          dwdzm = (w(i,j,k  )-w(i,j,km))*dzfi(k)
          !
          muxp = 0.25_rp*(mu(i,j,k)+mu(i,j,kp)+mu(ip,j ,kp)+mu(ip,j ,k) )
          muxm = 0.25_rp*(mu(i,j,k)+mu(i,j,kp)+mu(im,j ,kp)+mu(im,j ,k) )
          muyp = 0.25_rp*(mu(i,j,k)+mu(i,j,kp)+mu(i ,jp,kp)+mu(i ,jp,k) )
          muym = 0.25_rp*(mu(i,j,k)+mu(i,j,kp)+mu(i ,jm,kp)+mu(i ,jm,k) )
          muzp = mu(i,j,kp)
          muzm = mu(i,j,k )
          !
          rhoz   = 0.5_rp*(rho(  i,j,kp)+rho(  i,j,k))
          !
          dwdt(i,j,k) = &
                        - ( ug_w + vg_w + wg_w) + &
                        ( &
                        dxi*(    (dwdxp+dudzp)*muxp-(dwdxm+dudzm)*muxm ) + &
                        dyi*(    (dwdyp+dvdzp)*muyp-(dwdym+dvdzm)*muym ) + &
                        dzci(k)*((dwdzp+dwdzp)*muzp-(dwdzm+dwdzm)*muzm ) )/rhoz  

          !
        enddo
      enddo
    enddo
    !$acc end kernels 
    !$OMP END PARALLEL DO
    !
    return
  end subroutine momad_xyz_tw_fll
#else
  subroutine momad_xyz_sp(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dudt,dvdt,dwdt)  
    !
    use mod_param, only: rho_sp,mu_sp
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in ), dimension(1-nh_d:                ) :: dzci,dzfi
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(out), dimension(     1:,     1:,     1:) :: dudt,dvdt,dwdt
    !
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm, &
                uvip,uvim,vvjp,vvjm,wvkp,wvkm, &
                uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(rp) :: visc_sp
    integer  :: im,ip,jm,jp,km,kp,i,j,k
    !
    !@cuf attributes(managed) :: u, v, w, dzci, dzfi, dudt, dvdt, dwdt
    !    
    visc_sp = mu_sp/rho_sp
    !
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i + 1
          im = i - 1
          jp = j + 1
          jm = j - 1
          kp = k + 1
          km = k - 1
          !
#if defined(_TWOD)
          dudt(i,j,k) = 0._rp
#else
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          ! 1. r.h.s. of x-momentum
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          uuip  = 0.25_rp*( u(ip,j,k)+u(i,j,k) )*( u(ip,j ,k )+u(i,j ,k ) )
          uuim  = 0.25_rp*( u(im,j,k)+u(i,j,k) )*( u(im,j ,k )+u(i,j ,k ) )
          uvjp  = 0.25_rp*( u(i,jp,k)+u(i,j,k) )*( v(ip,j ,k )+v(i,j ,k ) )
          uvjm  = 0.25_rp*( u(i,jm,k)+u(i,j,k) )*( v(ip,jm,k )+v(i,jm,k ) )
          uwkp  = 0.25_rp*( u(i,j,kp)+u(i,j,k) )*( w(ip,j ,k )+w(i,j ,k ) )
          uwkm  = 0.25_rp*( u(i,j,km)+u(i,j,k) )*( w(ip,j ,km)+w(i,j ,km) )
          !
          dudxp = (u(ip,j ,k)-u(i ,j ,k))*dxi
          dudxm = (u(i ,j ,k)-u(im,j ,k))*dxi
          dudyp = (u(i ,jp,k)-u(i ,j ,k))*dyi
          dudym = (u(i ,j ,k)-u(i ,jm,k))*dyi
          dudzp = (u(i ,j,kp)-u(i ,j,k ))*dzci(k )
          dudzm = (u(i ,j,k )-u(i ,j,km))*dzci(km)
          !
          dudt(i,j,k) = &
                        dxi*(     -uuip + uuim ) + &
                        dyi*(     -uvjp + uvjm ) + &
                        dzfi(k)*( -uwkp + uwkm ) + &
                        ( &
                        dxi*(    dudxp-dudxm) + &
                        dyi*(    dudyp-dudym) + &
                        dzfi(k)*(dudzp-dudzm) )*visc_sp
#endif
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          ! 2. r.h.s. of y-momentum
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          uvip  = 0.25_rp*( u(i ,j,k)+u(i ,jp,k) )*( v(i,j,k )+v(ip,j ,k) )
          uvim  = 0.25_rp*( u(im,j,k)+u(im,jp,k) )*( v(i,j,k )+v(im,j ,k) )
          vvjp  = 0.25_rp*( v(i,j,k )+v(i,jp,k)  )*( v(i,j,k )+v(i ,jp,k) )
          vvjm  = 0.25_rp*( v(i,j,k )+v(i,jm,k)  )*( v(i,j,k )+v(i ,jm,k) )
          wvkp  = 0.25_rp*( w(i,j,k )+w(i,jp,k)  )*( v(i,j,kp)+v(i ,j ,k) )
          wvkm  = 0.25_rp*( w(i,j,km)+w(i,jp,km) )*( v(i,j,km)+v(i ,j ,k) )
          !
          dvdxp = (v(ip,j ,k)-v(i ,j ,k))*dxi
          dvdxm = (v(i ,j ,k)-v(im,j ,k))*dxi
          dvdyp = (v(i,jp ,k)-v(i ,j ,k))*dyi
          dvdym = (v(i,j  ,k)-v(i ,jm,k))*dyi
          dvdzp = (v(i,j,kp) -v(i,j,k  ))*dzci(k)
          dvdzm = (v(i,j,k ) -v(i,j,km ))*dzci(km)
          !
          dvdt(i,j,k) = &
                        dxi*(     -uvip + uvim ) + &
                        dyi*(     -vvjp + vvjm ) + & 
                        dzfi(k)*( -wvkp + wvkm ) + &
                        ( &
                        dxi*(    dvdxp-dvdxm) + &
                        dyi*(    dvdyp-dvdym) + &
                        dzfi(k)*(dvdzp-dvdzm) )*visc_sp
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          ! 3. r.h.s. of z-momentum
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          uwip  = 0.25_rp*( w(i,j,k)+w(ip,j,k) )*( u(i ,j ,k)+u(i ,j ,kp) )
          uwim  = 0.25_rp*( w(i,j,k)+w(im,j,k) )*( u(im,j ,k)+u(im,j ,kp) )
          vwjp  = 0.25_rp*( w(i,j,k)+w(i,jp,k) )*( v(i ,j ,k)+v(i ,j ,kp) )
          vwjm  = 0.25_rp*( w(i,j,k)+w(i,jm,k) )*( v(i ,jm,k)+v(i ,jm,kp) )
          wwkp  = 0.25_rp*( w(i,j,k)+w(i,j,kp) )*( w(i ,j ,k)+w(i ,j ,kp) )
          wwkm  = 0.25_rp*( w(i,j,k)+w(i,j,km) )*( w(i ,j ,k)+w(i ,j ,km) )
          !
          dwdxp = (w(ip,j,k )-w(i ,j,k))*dxi
          dwdxm = (w(i ,j,k )-w(im,j,k))*dxi
          dwdyp = (w(i,jp,k )-w(i,j ,k))*dyi
          dwdym = (w(i,j ,k )-w(i,jm,k))*dyi
          dwdzp = (w(i,j,kp )-w(i,j,k ))*dzfi(kp)
          dwdzm = (w(i,j,k  )-w(i,j,km))*dzfi(k)
          !
          dwdt(i,j,k) = &
                        dxi*(     -uwip + uwim ) + &
                        dyi*(     -vwjp + vwjm ) + &
                        dzci(k)*( -wwkp + wwkm ) + &
                        ( &
                        dxi*(    dwdxp-dwdxm) + &
                        dyi*(    dwdyp-dwdym) + &
                        dzci(k)*(dwdzp-dwdzm) )*visc_sp 

          !
        enddo
      enddo
    enddo
    !$acc end kernels 
    !$OMP END PARALLEL DO
    !
    return
  end subroutine momad_xyz_sp
#endif
  !
end module mod_mom
