module mod_spectra
! use mpi
use decomp_2d
use mod_param, only: lx,ly,lz,doSpectra, doFlux, spectra_deltai, &
                     flux_deltai, l_flux, nflux, dli, nbin, &
                     doSF, nl_sf, dl_sf,nmom_sf, sf_deltai, &
                     doRR, RR_deltai
use mod_common_mpi, only: myid, comm_Cart,ierr
use mod_fftw_spectra, only: init_fft,fftr2c,fftc2r
use mod_types
use mpi
! use mod_post, only: compHist3d
! use mod_load
! use mod_spectra
!
logical :: init_status = .false.
private
public post_globalT
contains

subroutine post_globalT(ng, dims, n, u, v, w,vof,ux,vx,wx,uy,vy,wy,uz,vz,wz,vofx,vofy,vofz, istep)
  implicit none 
  integer,   intent(in), dimension(3)             :: n
  integer,   intent(in), dimension(2)             :: dims
  integer,   intent(in)                           :: istep
  real(rp),  intent(in),    dimension(0: ,0:, 0:) :: u,v,w,vof
  real(rp),  intent(inout), dimension(0: ,0:, 0:) :: ux,vx,wx,uy,vy,wy,uz,vz,wz,vofx,vofy,vofz
  integer,   dimension(3)                         :: nf, ng
  integer(8)                                      :: plan_r2c, plan_c2r


  ng(:) = n(:)
  ng(1:2) = ng(1:2)*dims(1:2)

  !!!!! Postprocessing with pencil-stencil oriented along Z
  uz(:,:,:) = u(1:n(1),1:n(2),1:n(3))
  vz(:,:,:) = v(1:n(1),1:n(2),1:n(3))
  wz(:,:,:) = w(1:n(1),1:n(2),1:n(3))

  if ((mod(istep,spectra_deltai).eq.0).and.(doSpectra.eqv..true.))then
    call init_fft(ng(3),plan_r2c, plan_c2r) !for the case we study, itot=jtot=ktot
    call compute_spectrum(istep,3,uz,uz,ng(3),lz,ng(1)/dims(1),ng(2)/dims(2),'11',dims,plan_r2c, plan_c2r)
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
    call init_fft(ng(3),plan_r2c, plan_c2r) !for the case we study, itot=jtot=ktot
    call compute_spectrum(istep,3,vz,vz,ng(3),lz,ng(1)/dims(1),ng(2)/dims(2),'22',dims,plan_r2c, plan_c2r)
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
    call init_fft(ng(3),plan_r2c, plan_c2r) !for the case we study, itot=jtot=ktot
    call compute_spectrum(istep,3,wz,wz,ng(3),lz,ng(1)/dims(1),ng(2)/dims(2),'33',dims,plan_r2c, plan_c2r)
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
  endif

  if ((mod(istep,flux_deltai).eq.0).and.(doFlux.eqv..true.))then
    call init_fft(ng(3),plan_r2c, plan_c2r) !for the case we study, itot=jtot=ktot
    nf = n
    call compute_flux(istep,3,uz,vz,wz,nf,ng,(/ng(1)/dims(1),ng(2)/dims(2),n(3)/2+1/),lz,dims,plan_r2c, plan_c2r)
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
  endif

  if ((mod(istep,sf_deltai).eq.0).and.(doSF.eqv..true.)) then
    call structureFunction(istep, 3, uz,vz,wz, ng,(/ng(1)/dims(1),ng(2)/dims(2),n(3)/),dims)
  endif

  if ((mod(istep,RR_deltai).eq.0).and.(doRR.eqv..true.)) then
    vofz(:,:,:) = vof(1:n(1),1:n(2),1:n(3))
    call decoupledAutocorr(istep, 3, uz,vz,wz,vofz, ng,(/ng(1)/dims(1),ng(2)/dims(2),n(3)/),dims)
  endif

  !!!!! Postprocessing with pencil-stencil oriented along Y
  call transpose_z_to_y(uz,uy)
  call transpose_z_to_y(vz,vy)
  call transpose_z_to_y(wz,wy)
  if ((mod(istep,spectra_deltai).eq.0).and.(doSpectra.eqv..true.))then
    call init_fft(ng(2),plan_r2c, plan_c2r) !for the case we study, itot=jtot=ktot
    call compute_spectrum(istep,2,uy,uy,ng(2),ly,ng(1)/dims(1),ng(3)/dims(2),'11',dims,plan_r2c, plan_c2r)
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
    call init_fft(ng(2),plan_r2c, plan_c2r) !for the case we study, itot=jtot=ktot
    call compute_spectrum(istep,2,vy,vy,ng(2),ly,ng(1)/dims(1),ng(3)/dims(2),'22',dims,plan_r2c, plan_c2r)
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
    call init_fft(ng(2),plan_r2c, plan_c2r) !for the case we study, itot=jtot=ktot
    call compute_spectrum(istep,2,wy,wy,ng(2),ly,ng(1)/dims(1),ng(3)/dims(2),'33',dims,plan_r2c, plan_c2r)
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
  endif

  if ((mod(istep,flux_deltai).eq.0).and.(doFlux.eqv..true.))then
    call init_fft(ng(2),plan_r2c, plan_c2r) !for the case we study, itot=jtot=ktot
    nf = (/ng(1)/dims(1),ng(2),n(3)/dims(2)/)
    call compute_flux(istep,2,uy,vy,wy,nf,ng,(/ng(1)/dims(1),ng(2)/2+1,n(3)/dims(2)/),ly,dims,plan_r2c, plan_c2r)
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
  endif

  if ((mod(istep,sf_deltai).eq.0).and.(doSF.eqv..true.)) then
    call structureFunction(istep, 2, uy,vy,wy, ng,(/ng(1)/dims(1),ng(2),n(3)/dims(2)/),dims)
  endif

  if ((mod(istep,RR_deltai).eq.0).and.(doRR.eqv..true.)) then
    call transpose_z_to_y(vofz,vofy)
    call decoupledAutocorr(istep, 2, uy,vy,wy,vofy, ng,(/ng(1)/dims(1),ng(2),n(3)/dims(2)/),dims)
  endif

  !!!!! Postprocessing with pencil-stencil oriented along Y
  call transpose_y_to_x(uy,ux)
  call transpose_y_to_x(vy,vx)
  call transpose_y_to_x(wy,wx)
  
  if ((mod(istep,spectra_deltai).eq.0).and.(doSpectra.eqv..true.))then
    call init_fft(ng(1),plan_r2c, plan_c2r) !for the case we study, itot=jtot=ktot
    call compute_spectrum(istep,1,ux,ux,ng(1),lx,ng(2)/dims(1),ng(3)/dims(2),'11',dims,plan_r2c, plan_c2r)
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
    call init_fft(ng(1),plan_r2c, plan_c2r) !for the case we study, itot=jtot=ktot
    call compute_spectrum(istep,1,vx,vx,ng(1),lx,ng(2)/dims(1),ng(3)/dims(2),'22',dims,plan_r2c, plan_c2r)
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
    call init_fft(ng(1),plan_r2c, plan_c2r) !for the case we study, itot=jtot=ktot
    call compute_spectrum(istep,1,wx,wx,ng(1),lx,ng(2)/dims(1),ng(3)/dims(2),'33',dims,plan_r2c, plan_c2r)
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
  endif

  if ((mod(istep,flux_deltai).eq.0).and.(doFlux.eqv..true.))then
    call init_fft(ng(1),plan_r2c, plan_c2r) !for the case we study, itot=jtot=ktot
    nf = (/ng(1),ng(2)/dims(1),n(3)/dims(2)/)
    call compute_flux(istep,1,ux,vx,wx,nf,ng,(/ng(1)/2+1,ng(2)/dims(1),n(3)/dims(2)/),lx,dims,plan_r2c, plan_c2r)
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
  endif

  if ((mod(istep,sf_deltai).eq.0).and.(doSF.eqv..true.)) then
    call structureFunction(istep, 1, ux,vx,wx, ng,(/ng(1),ng(2)/dims(1),n(3)/dims(2)/),dims)
  endif  

  if ((mod(istep,RR_deltai).eq.0).and.(doRR.eqv..true.)) then
    call transpose_y_to_x(vofy,vofx)
    call decoupledAutocorr(istep, 1, ux,vx,wx,vofx, ng,(/ng(1),ng(2)/dims(1),n(3)/dims(2)/),dims)
  endif

end subroutine post_globalT



subroutine compute_spectrum(istep,idir,a,b,nhet,lhet,nhom1,nhom2,index,dims,plan_r2c, plan_c2r)
  ! use mod_common_post
  implicit none
  !i,j -> components of the velocity
  integer, intent(in) :: istep,idir,nhet,nhom1,nhom2
  integer(8), intent(in)  ::plan_r2c, plan_c2r
  integer, intent(in), dimension(2):: dims
  real(rp), dimension(1:,1:,1:), intent(in) :: a,b
  character(len=2), intent(in) :: index
  real(rp), intent(in) :: lhet
  real(rp), dimension(nhet) :: Rab
  complex(rp), dimension(nhet/2+1) :: acmplx,bcmplx,Eabcmplx
  real(rp), dimension(nhet/2+1) :: Eab
  real(rp), dimension(nhet/2+1) :: Eab_av,Eab_av_all
  real(rp), dimension(nhet) :: Rab_av,Rab_av_all
  character(len=1) :: idirchar
  character(len=7) :: fldnum
  integer :: i,j,k,ii, error
  real(rp) :: wavenum
  real(rp) :: secondderiv,lab,lambda,var,delta,deltak
  real(rp), parameter :: pi = acos(-1.)
  real(rp) :: wavenummax
  real(rp) :: Eab_av_sum
  !real,dimension(nhet) :: a_aux,b_aux
  !
  !
  !

  Rab_av(:) = 0.
  Eab_av(:)  = 0.
  select case(idir)
  case(1)
    do k = 1,nhom2
      do j = 1,nhom1
          call fftr2c(nhet,a(1:nhet,j,k),acmplx(:),plan_r2c)
          call fftr2c(nhet,b(1:nhet,j,k),bcmplx(:),plan_r2c)
          Eabcmplx(:) = (acmplx(:)/(1.*nhet))*conjg(bcmplx(:)/(1.*nhet))
          Eab(:) = real(Eabcmplx(:))
         ! call fftc2r(nhet,Eabcmplx,Rab,idir)
          Eab_av(:) = Eab_av(:) + Eab(:)
         ! Rab_av(:) = Rab_av(:) + Rab(:)
      enddo
    enddo
  case(2)
    do k = 1,nhom2
      do i = 1,nhom1
          call fftr2c(nhet,a(i,1:nhet,k),acmplx(:),plan_r2c)
          call fftr2c(nhet,b(i,1:nhet,k),bcmplx(:),plan_r2c)
          Eabcmplx(:) = (acmplx(:)/(1.*nhet))*conjg(bcmplx(:)/(1.*nhet))
          Eab(:) = real(Eabcmplx(:))
         ! call fftc2r(nhet,Eabcmplx,Rab,idir)
          Eab_av(:) = Eab_av(:) + Eab(:)
         ! Rab_av(:) = Rab_av(:) + Rab(:)
      enddo
    enddo
  case(3)
    do j = 1,nhom2
      do i = 1,nhom1
          call fftr2c(nhet,a(i,j,1:nhet),acmplx(:),plan_r2c)
          call fftr2c(nhet,b(i,j,1:nhet),bcmplx(:),plan_r2c)
          Eabcmplx(:) = (acmplx(:)/(1.*nhet))*conjg(bcmplx(:)/(1.*nhet))
          Eab(:) = real(Eabcmplx(:))
         ! call fftc2r(nhet,Eabcmplx,Rab,idir)
          Eab_av(:) = Eab_av(:) + Eab(:)
         ! Rab_av(:) = Rab_av(:) + Rab(:)
      enddo
    enddo
  end select
  Eab_av(:) = Eab_av(:)/(1.*nhom1*nhom2)
  ! Rab_av(:) = Rab_av(:)/(1.*nhom1*nhom2)

  call mpi_allreduce(Eab_av(1),Eab_av_all(1),nhet/2+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,error)
  call mpi_allreduce(Rab_av(1),Rab_av_all(1),nhet,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,error)
  Eab_av(:) = Eab_av_all(:)/product(dims)
  ! Rab_av(:) = Rab_av_all(:)/product(dims)

  deltak = 2*pi/lhet
  wavenummax = nhet/2*deltak
  ! Eab_av(1) = Eab_av(1)
  ! do ii=2,nhet/2+1
  !   Eab_av(ii) = 2*Eab_av(ii)
  ! enddo
  ! lab = Rab_av(1)
  ! delta = lhet/nhet
  ! var = Rab_av(1)
  ! secondderiv = (Rab_av(nhet)-2*Rab_av(1)*2+Rab_av(2))/(delta**2.*var) ! periodicity 
  ! lambda = sqrt(-1./(secondderiv))


  if(myid.eq.0) then
  write(fldnum,'(i7.7)') istep
  write(idirchar,'(i1.1)') idir
  open(20, file='data/post/spectra/E'//index//'dir'//idirchar//'fld'//fldnum//'.out')
  write(20,'(A,I1,A)') '# E'//index//'(k',idir,')'
  write(20,'(A)') '# k, E'
  wavenum = 0.
  ! Eab_av_sum = sum(Eab_av(:))

  !write(*,*) 'small = ', var - Eab_av_sum/(1.*Nhet)

  do ii=0,nhet/2
    write(20,'(2E15.7)') wavenum,Eab_av(ii+1)
    wavenum = wavenum+deltak
  enddo
  close(20)
  !
  ! write
  !
  !open(20, 'data/post/spectra/R'//index//'dir'//idirchar//'fld'//fldnum//'.out')
  !write(20,'(A,I1,A)') '# R'//index//'(r',idir,')'
  !!write(20,'(A,E15.7)') '# Integral Length scale: ', lab
  !write(20,'(A,E15.7)') '# Taylor Microscale: ', lambda
  !write(20,'(A)') '# r, R'
  !
  !ii=0
  !write(20,'(2E15.7)') 0.*delta,Rab_av(ii+1)/var
  !do ii=1,nhet/2
  !  write(20,'(2E15.7)') 1.*ii*delta,0.5*(Rab_av(ii+1)+Rab_av(nhet-ii+1))/var
  !enddo
  !close(20)
  endif
  call mpi_barrier(comm_cart,error)
  !
  return
end subroutine compute_spectrum
!
!

subroutine compute_flux(istep,dir,u,v,w,n,ng,np,lhet,dims,plan_r2c, plan_c2r)
  implicit none
  real(rp),     intent(in), dimension(1:,1:,1:)                           :: u,v,w
  real(rp),     intent(in)                                                :: lhet
  integer(8),   intent(in)                                                :: plan_r2c, plan_c2r
  integer,      intent(in), dimension(2)                                  :: dims
  integer,      intent(in), dimension(3)                                  :: ng, np, n
  integer,      intent(in)                                                :: istep, dir
  complex(rp),              dimension(np(1),np(2),np(3))                  :: u_ft, v_ft, w_ft
  complex(rp),              dimension(np(1),np(2),np(3))                  :: uv_ft, uw_ft, vw_ft, uu_ft, vv_ft, ww_ft
  real(rp),                 dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1)         :: u_l, v_l, w_l
  real(rp),                 dimension(n(1),n(2),n(3))                     :: uv_l, uw_l, vw_l, uu_l, vv_l, ww_l, uc, vc, wc
  real(rp),                 dimension(3)                                  :: dudi_l, dvdi_l, dwdi_l
  integer                                                                 :: i, j, k, li, ii, nhom1,nhom2, nhet
  integer                                                                 :: ip, jp, kp, im, jm, km
  real(rp),                 dimension(np(dir))                            :: G_l, wavenum      
  real(rp),                 dimension(n(1),n(2),n(3),nflux)               :: Pi_l
  real(rp),                 dimension(nflux)                              :: Pi_lsum
  real(rp),     parameter                                                 :: pi = acos(-1.)
  real(rp)                                                                :: deltak
  character(len=1)                                                        :: idirchar
  character(len=1)                                                        :: listr
  character(len=7)                                                        :: fldnum


  do k=1,n(3)
    km = k - 1
    do j=1,n(2)
      jm = j - 1
      do i=1,n(1)
        im = i - 1
        uc(i,j,k) = 0.5*(u(i,j,k) + u(im,j,k))
        vc(i,j,k) = 0.5*(v(i,j,k) + v(i,jm,k))
        wc(i,j,k) = 0.5*(w(i,j,k) + w(i,j,km))
      enddo
    enddo
  enddo


  do li=1,nflux
    l_flux(li) = li/(nflux*1.0)*lhet
    select case (dir)
      case(1)
        nhom1 = n(2)
        nhom2 = n(3)
        nhet  = n(1)

        deltak = 2*pi/lhet
        do ii=1,nhet/2+1
         wavenum(ii) = (ii-1)*deltak
        enddo
          
          G_l = exp(-wavenum**2 * l_flux(li)**2/2.)
          do j = 1,n(2)
            do k = 1,n(3)
                call fftr2c(nhet, uc(1:nhet,j,k)                    ,u_ft  (:,j,k),plan_r2c)
                call fftr2c(nhet, vc(1:nhet,j,k)                    ,v_ft  (:,j,k),plan_r2c)
                call fftr2c(nhet, wc(1:nhet,j,k)                    ,w_ft  (:,j,k),plan_r2c)
                call fftr2c(nhet, uc(1:nhet,j,k)*vc(1:nhet,j,k)     ,uv_ft (:,j,k),plan_r2c)
                call fftr2c(nhet, uc(1:nhet,j,k)*wc(1:nhet,j,k)     ,uw_ft (:,j,k),plan_r2c)
                call fftr2c(nhet, wc(1:nhet,j,k)*vc(1:nhet,j,k)     ,vw_ft (:,j,k),plan_r2c)
                call fftr2c(nhet, uc(1:nhet,j,k)*uc(1:nhet,j,k)     ,uu_ft (:,j,k),plan_r2c)
                call fftr2c(nhet, vc(1:nhet,j,k)*vc(1:nhet,j,k)     ,vv_ft (:,j,k),plan_r2c)
                call fftr2c(nhet, wc(1:nhet,j,k)*wc(1:nhet,j,k)     ,ww_ft (:,j,k),plan_r2c)
                u_ft (:,j,k)  = u_ft  (:,j,k)*G_l
                v_ft (:,j,k)  = v_ft  (:,j,k)*G_l
                w_ft (:,j,k)  = w_ft  (:,j,k)*G_l
                uv_ft (:,j,k) = uv_ft (:,j,k)*G_l
                uw_ft (:,j,k) = uw_ft (:,j,k)*G_l
                vw_ft (:,j,k) = vw_ft (:,j,k)*G_l
                uu_ft (:,j,k) = uu_ft (:,j,k)*G_l
                vv_ft (:,j,k) = vv_ft (:,j,k)*G_l
                ww_ft (:,j,k) = ww_ft (:,j,k)*G_l
               call fftc2r(nhet, u_ft(:,j,k)  ,u_l (1:nhet,j,k) ,plan_c2r)
               call fftc2r(nhet, v_ft(:,j,k)  ,v_l (1:nhet,j,k) ,plan_c2r)
               call fftc2r(nhet, w_ft(:,j,k)  ,w_l (1:nhet,j,k) ,plan_c2r)
               call fftc2r(nhet, uv_ft(:,j,k) ,uv_l(1:nhet,j,k) ,plan_c2r)
               call fftc2r(nhet, uw_ft(:,j,k) ,uw_l(1:nhet,j,k) ,plan_c2r)
               call fftc2r(nhet, vw_ft(:,j,k) ,vw_l(1:nhet,j,k) ,plan_c2r)
               call fftc2r(nhet, uu_ft(:,j,k) ,uu_l(1:nhet,j,k) ,plan_c2r)
               call fftc2r(nhet, vv_ft(:,j,k) ,vv_l(1:nhet,j,k) ,plan_c2r)
               call fftc2r(nhet, ww_ft(:,j,k) ,ww_l(1:nhet,j,k) ,plan_c2r)
            enddo
          enddo

      case(2)
        nhom1 = n(1)
        nhom2 = n(3)
        nhet  = n(2)
        deltak = 2*pi/lhet
        do ii=1,nhet/2+1
         wavenum(ii) = (ii-1)*deltak
        enddo
        
          G_l = exp(-wavenum**2 *l_flux(li)**2/2.)
          do k = 1,n(3)
            do i = 1,n(1)
                call fftr2c(nhet, uc(i,1:nhet,k)                   ,u_ft (i,:,k),plan_r2c)
                call fftr2c(nhet, vc(i,1:nhet,k)                   ,v_ft (i,:,k),plan_r2c)
                call fftr2c(nhet, wc(i,1:nhet,k)                   ,w_ft (i,:,k),plan_r2c)
                call fftr2c(nhet, uc(i,1:nhet,k)*vc(i,1:nhet,k)    ,uv_ft (i,:,k),plan_r2c)
                call fftr2c(nhet, uc(i,1:nhet,k)*wc(i,1:nhet,k)    ,uw_ft (i,:,k),plan_r2c)
                call fftr2c(nhet, wc(i,1:nhet,k)*vc(i,1:nhet,k)    ,vw_ft (i,:,k),plan_r2c)
                call fftr2c(nhet, uc(i,1:nhet,k)*uc(i,1:nhet,k)    ,uu_ft (i,:,k),plan_r2c)
                call fftr2c(nhet, vc(i,1:nhet,k)*vc(i,1:nhet,k)    ,vv_ft (i,:,k),plan_r2c)
                call fftr2c(nhet, wc(i,1:nhet,k)*wc(i,1:nhet,k)    ,ww_ft (i,:,k),plan_r2c)
                u_ft (i,:,k)  = u_ft  (i,:,k)*G_l
                v_ft (i,:,k)  = v_ft  (i,:,k)*G_l
                w_ft (i,:,k)  = w_ft  (i,:,k)*G_l
                uv_ft (i,:,k) = uv_ft (i,:,k)*G_l
                uw_ft (i,:,k) = uw_ft (i,:,k)*G_l
                vw_ft (i,:,k) = vw_ft (i,:,k)*G_l
                uu_ft (i,:,k) = uu_ft (i,:,k)*G_l
                vv_ft (i,:,k) = vv_ft (i,:,k)*G_l
                ww_ft (i,:,k) = ww_ft (i,:,k)*G_l
               call fftc2r(nhet, u_ft(i,:,k)  ,u_l (i,1:nhet,k) ,plan_c2r)
               call fftc2r(nhet, v_ft(i,:,k)  ,v_l (i,1:nhet,k) ,plan_c2r)
               call fftc2r(nhet, w_ft(i,:,k)  ,w_l (i,1:nhet,k) ,plan_c2r)
               call fftc2r(nhet, uv_ft(i,:,k) ,uv_l(i,1:nhet,k) ,plan_c2r)
               call fftc2r(nhet, uw_ft(i,:,k) ,uw_l(i,1:nhet,k) ,plan_c2r)
               call fftc2r(nhet, vw_ft(i,:,k) ,vw_l(i,1:nhet,k) ,plan_c2r)
               call fftc2r(nhet, uu_ft(i,:,k) ,uu_l(i,1:nhet,k) ,plan_c2r)
               call fftc2r(nhet, vv_ft(i,:,k) ,vv_l(i,1:nhet,k) ,plan_c2r)
               call fftc2r(nhet, ww_ft(i,:,k) ,ww_l(i,1:nhet,k) ,plan_c2r)
            enddo
          enddo

      case(3)
      nhom1 = n(1)
      nhom2 = n(2)
      nhet  = n(3)
      deltak = 2*pi/lhet
      do ii=1,nhet/2+1
       wavenum(ii) = (ii-1)*deltak
      enddo
      
        G_l = exp(-wavenum**2 *l_flux(li)**2/2.)

        
        do j = 1,n(2)
          do i = 1,n(1)
              call fftr2c(nhet, uc(i,j,1:nhet)                    ,u_ft (i,j,:),plan_r2c)
              call fftr2c(nhet, vc(i,j,1:nhet)                    ,v_ft (i,j,:),plan_r2c)
              call fftr2c(nhet, wc(i,j,1:nhet)                    ,w_ft (i,j,:),plan_r2c)
              call fftr2c(nhet, uc(i,j,1:nhet)*vc(i,j,1:nhet)     ,uv_ft(i,j,:),plan_r2c)
              call fftr2c(nhet, uc(i,j,1:nhet)*wc(i,j,1:nhet)     ,uw_ft(i,j,:),plan_r2c)
              call fftr2c(nhet, wc(i,j,1:nhet)*vc(i,j,1:nhet)     ,vw_ft(i,j,:),plan_r2c)
              call fftr2c(nhet, uc(i,j,1:nhet)*uc(i,j,1:nhet)     ,uu_ft(i,j,:),plan_r2c)
              call fftr2c(nhet, vc(i,j,1:nhet)*vc(i,j,1:nhet)     ,vv_ft(i,j,:),plan_r2c)
              call fftr2c(nhet, wc(i,j,1:nhet)*wc(i,j,1:nhet)     ,ww_ft(i,j,:),plan_r2c)
              u_ft (i,j,:)  = u_ft  (i,j,:)*G_l
              v_ft (i,j,:)  = v_ft  (i,j,:)*G_l
              w_ft (i,j,:)  = w_ft  (i,j,:)*G_l
              uv_ft (i,j,:) = uv_ft (i,j,:)*G_l
              uw_ft (i,j,:) = uw_ft (i,j,:)*G_l
              vw_ft (i,j,:) = vw_ft (i,j,:)*G_l
              uu_ft (i,j,:) = uu_ft (i,j,:)*G_l
              vv_ft (i,j,:) = vv_ft (i,j,:)*G_l
              ww_ft (i,j,:) = ww_ft (i,j,:)*G_l
             call fftc2r(nhet, u_ft (i,j,:)  ,u_l (i,j,1:nhet) ,plan_c2r)
             call fftc2r(nhet, v_ft (i,j,:)  ,v_l (i,j,1:nhet) ,plan_c2r)
             call fftc2r(nhet, w_ft (i,j,:)  ,w_l (i,j,1:nhet) ,plan_c2r)
             call fftc2r(nhet, uv_ft(i,j,:)  ,uv_l(i,j,1:nhet) ,plan_c2r)
             call fftc2r(nhet, uw_ft(i,j,:)  ,uw_l(i,j,1:nhet) ,plan_c2r)
             call fftc2r(nhet, vw_ft(i,j,:)  ,vw_l(i,j,1:nhet) ,plan_c2r)
             call fftc2r(nhet, uu_ft(i,j,:)  ,uu_l(i,j,1:nhet) ,plan_c2r)
             call fftc2r(nhet, vv_ft(i,j,:)  ,vv_l(i,j,1:nhet) ,plan_c2r)
             call fftc2r(nhet, ww_ft(i,j,:)  ,ww_l(i,j,1:nhet) ,plan_c2r)
          enddo
        enddo

    end select

  

    call updthalo(n,1,u_l ,dir)
    call updthalo(n,2,u_l ,dir)
    call updthalo(n,3,u_l ,dir)
    call updthalo(n,1,v_l ,dir)
    call updthalo(n,2,v_l ,dir)
    call updthalo(n,3,v_l ,dir)
    call updthalo(n,1,w_l ,dir)
    call updthalo(n,2,w_l ,dir)
    call updthalo(n,3,w_l ,dir)

    select case (dir)
      case(1)
        do j=0,n(2)+1
          do k=0,n(3)+1
            u_l(0,j,k)         = u_l(n(1),j,k)
            u_l(n(1)+1,j,k)    = u_l(1,j,k)
            v_l(0,j,k)         = v_l(n(1),j,k)
            v_l(n(1)+1,j,k)    = v_l(1,j,k)
            w_l(0,j,k)         = w_l(n(1),j,k)
            w_l(n(1)+1,j,k)    = w_l(1,j,k)
          enddo
        enddo
      case(2)
        do i=0,n(1)+1
          do k=0,n(3)+1
            u_l(i,0,k)         = u_l(i,n(2),k)
            u_l(i,n(2)+1,k)    = u_l(i,1,k)
            v_l(i,0,k)         = v_l(i,n(2),k)
            v_l(i,n(2)+1,k)    = v_l(i,1,k)
            w_l(i,0,k)         = w_l(i,n(2),k)
            w_l(i,n(2)+1,k)    = w_l(i,1,k)
          enddo
        enddo
      case(3)
        do i=0,n(1)+1
          do j=0,n(2)+1
            u_l(i,j,0)         = u_l(i,j,n(3))
            u_l(i,j,n(3)+1)    = u_l(i,j,1)
            v_l(i,j,0)         = v_l(i,j,n(3))
            v_l(i,j,n(3)+1)    = v_l(i,j,1)
            w_l(i,j,0)         = w_l(i,j,n(3))
            w_l(i,j,n(3)+1)    = w_l(i,j,1)
          enddo
        enddo
    end select



    Pi_lsum(li) = 0.0

    do k=1,n(3)
      kp = k + 1
      km = k - 1
      do j=1,n(2)
        jp = j + 1
        jm = j - 1
        do i=1,n(1)

          ip = i + 1
          im = i - 1
          
          dudi_l(1) = ( u_l(i,j,k) - u_l(im,j,k) )*dli(1)
          dudi_l(2) = 0.25*( u_l(i,j,k) + u_l(i,jp,k) + u_l(im,j,k) + u_l(im,jp,k) )*dli(2) - &
                      0.25*( u_l(i,j,k) + u_l(i,jm,k) + u_l(im,j,k) + u_l(im,jm,k) )*dli(2)
          dudi_l(3) = 0.25*( u_l(i,j,k) + u_l(i,j,kp) + u_l(im,j,k) + u_l(im,j,kp) )*dli(3) - &
                      0.25*( u_l(i,j,k) + u_l(i,j,km) + u_l(im,j,k) + u_l(im,j,km) )*dli(3)

          dvdi_l(1) = 0.25*( v_l(i,j,k) + v_l(ip,j,k) + v_l(i,jm,k) + v_l(ip,jm,k) )*dli(1) - &
                      0.25*( v_l(i,j,k) + v_l(im,j,k) + v_l(i,jm,k) + v_l(im,jm,k) )*dli(1)
          dvdi_l(2) = ( v_l(i,j,k) - v_l(i,jm,k) )*dli(2)
          dvdi_l(3) = 0.25*( v_l(i,j,k) + v_l(i,j,kp) + v_l(i,jm,k) + v_l(i,jm,kp) )*dli(3) - &
                      0.25*( v_l(i,j,k) + v_l(i,j,km) + v_l(i,jm,k) + v_l(i,jm,km) )*dli(3)

          dwdi_l(1) = 0.25*( w_l(i,j,k) + w_l(ip,j,k) + w_l(i,j,km) + w_l(ip,j,km) )*dli(1) - &
                      0.25*( w_l(i,j,k) + w_l(im,j,k) + w_l(i,j,km) + w_l(im,j,km) )*dli(1)
          dwdi_l(2) = 0.25*( w_l(i,j,k) + w_l(i,jp,k) + w_l(i,j,km) + w_l(i,jp,km) )*dli(2) - &
                      0.25*( w_l(i,j,k) + w_l(i,jm,k) + w_l(i,j,km) + w_l(i,jm,km) )*dli(2)
          dwdi_l(3) = ( w_l(i,j,k) - w_l(i,j,km) )*dli(3)


          Pi_l(i,j,k,li)  =  dudi_l(1)            * (uu_l(i,j,k) - (u_l(i,j,k)*u_l(i,j,k))) +&
                            (dudi_l(2)+dvdi_l(1)) * (uv_l(i,j,k) - (u_l(i,j,k)*v_l(i,j,k))) +&
                            (dudi_l(3)+dwdi_l(1)) * (uw_l(i,j,k) - (u_l(i,j,k)*w_l(i,j,k))) +&
                             dvdi_l(2)            * (vv_l(i,j,k) - (v_l(i,j,k)*v_l(i,j,k))) +&
                            (dvdi_l(3)+dwdi_l(2)) * (vw_l(i,j,k) - (v_l(i,j,k)*w_l(i,j,k))) +&
                             dwdi_l(3)            * (ww_l(i,j,k) - (w_l(i,j,k)*w_l(i,j,k)))      
          ! Pi_l(i,j,k,li)  =  dudi_l(1)             +&
          !                   (dudi_l(2)+dvdi_l(1))  +&
          !                   (dudi_l(3)+dwdi_l(1))  +&
          !                    dvdi_l(2)             +&
          !                   (dvdi_l(3)+dwdi_l(2))  +&
          !                    dwdi_l(3)                  

          ! Pi_l(i,j,k,li)  =  (uu_l(i,j,k) - (u_l(i,j,k)*u_l(i,j,k))) +&
          !                    (uv_l(i,j,k) - (u_l(i,j,k)*v_l(i,j,k))) +&
          !                    (uw_l(i,j,k) - (u_l(i,j,k)*w_l(i,j,k))) +&
          !                    (vv_l(i,j,k) - (v_l(i,j,k)*v_l(i,j,k))) +&
          !                    (vw_l(i,j,k) - (v_l(i,j,k)*w_l(i,j,k))) +&
          !                    (ww_l(i,j,k) - (w_l(i,j,k)*w_l(i,j,k))) 

          ! Pi_l(i,j,k,li)  =  uu_l(i,j,k)+& ! - (u_l(i,j,k)*u_l(i,j,k))) +&
          !                    uv_l(i,j,k)+& ! - (u_l(i,j,k)*v_l(i,j,k))) +&
          !                    uw_l(i,j,k)+&! - (u_l(i,j,k)*w_l(i,j,k))) +&
          !                    vv_l(i,j,k)+&! - (v_l(i,j,k)*v_l(i,j,k))) +&
          !                    vw_l(i,j,k)+&! - (v_l(i,j,k)*w_l(i,j,k))) +&
          !                    ww_l(i,j,k)! - (w_l(i,j,k)*w_l(i,j,k))) 
          ! Pi_l(i,j,k,li) = uw_l(i,j,k)
          Pi_lsum(li) = Pi_lsum(li) + Pi_l(i,j,k,li)


        enddo
      enddo
    enddo
    Pi_lsum(li) = Pi_lsum(li)/product(n)
    call mpi_allreduce(MPI_IN_PLACE,Pi_lsum(li) ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    Pi_lsum(li) = Pi_lsum(li)/product(dims)
    ! call compHist3d(Pi_l(:,:,:,li),n,nbin, 0,istep, varname)
  enddo
  
  if(myid.eq.0) then
    write(fldnum,'(i7.7)') istep
    write(idirchar,'(i1.1)') dir
    open(20, file='data/post/flux/pi'//idirchar//'fld'//fldnum//'.out')
    ! write(20,'(A,I1,A)') '# E'//index//'(k',idir,')'
    write(20,'(A)') '# l, pi'
    wavenum = 0.
    ! Eab_av_sum = sum(Eab_av(:))

    !write(*,*) 'small = ', var - Eab_av_sum/(1.*Nhet)

    do li=1,nflux
      write(20,'(2E15.7)') l_flux(li),Pi_lsum(li)
    enddo
    close(20)
  endif

  ! deallocate(G_l, wavenum)
  call mpi_barrier(comm_cart,ierr)
  return
end subroutine compute_flux
  

subroutine updthalo(n,idir,p,stencil_dir)
    use mpi
    use mod_common_mpi, only: myid, halo
    use mod_types
    use decomp_2d
    implicit none
    integer , dimension(3), intent(in) :: n
    integer , intent(in) :: idir, stencil_dir
    real(rp), dimension(0:,0:,0:), intent(inout) :: p
    integer, dimension(0:1,3) :: nb
    integer, dimension(3) :: nh
    integer               :: ierr,comm_cart1
    integer                           :: status(MPI_STATUS_SIZE)
    !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  this subroutine updates the haloalos that store info
    !  from the neighboring computational sub-domain
    !
    select case(stencil_dir)
    case(1)
      comm_cart1 = DECOMP_2D_COMM_CART_X
      nb(0,1) = myid
      nb(1,1) = myid
      call MPI_CART_SHIFT(comm_cart1,0,1,nb(0,2),nb(1,2),ierr)
      call MPI_CART_SHIFT(comm_cart1,1,1,nb(0,3),nb(1,3),ierr)
    case(2)
      comm_cart1 = DECOMP_2D_COMM_CART_Y
      call MPI_CART_SHIFT(comm_cart1,0,1,nb(0,1),nb(1,1),ierr)
      nb(0,2) = myid
      nb(1,2) = myid
      call MPI_CART_SHIFT(comm_cart1,1,1,nb(0,3),nb(1,3),ierr)
    case(3)
      comm_cart1 = DECOMP_2D_COMM_CART_Z
      call MPI_CART_SHIFT(comm_cart1,0,1,nb(0,1),nb(1,1),ierr)
      call MPI_CART_SHIFT(comm_cart1,1,1,nb(0,2),nb(1,2),ierr)
      nb(0,3) = myid
      nb(1,3) = myid
    end select

    
    select case(idir)
    case(1) ! x direction
    if (stencil_dir.ne.1) then
          call MPI_SENDRECV(p(1     ,0,0),1,halo(stencil_dir,idir),nb(0,idir),0, &
                        p(n(1)+1,0,0),1,halo(stencil_dir,idir),nb(1,idir),0, &
                        comm_cart1,status,ierr)
          call MPI_SENDRECV(p(n(1)  ,0,0),1,halo(stencil_dir,idir),nb(1,idir),0, &
                        p(0     ,0,0),1,halo(stencil_dir,idir),nb(0,idir),0, &
                        comm_cart1,status,ierr)
    endif

    case(2) ! y direction
    if (stencil_dir.ne.2) then 
      call MPI_SENDRECV(p(0,1     ,0),1,halo(stencil_dir,idir),nb(0,idir),0, &
                        p(0,n(2)+1,0),1,halo(stencil_dir,idir),nb(1,idir),0, &
                        comm_cart1,status,ierr)
      call MPI_SENDRECV(p(0,n(2)  ,0),1,halo(stencil_dir,idir),nb(1,idir),0, &
                        p(0,0     ,0),1,halo(stencil_dir,idir),nb(0,idir),0, &
                        comm_cart1,status,ierr)
    endif
    case(3) ! z direction
    if (stencil_dir.ne.3) then
      call MPI_SENDRECV(p(0,0,1     ),1,halo(stencil_dir,idir),nb(0,idir),0, &
                        p(0,0,n(3)+1),1,halo(stencil_dir,idir),nb(1,idir),0, &
                        comm_cart1,status,ierr)
      call MPI_SENDRECV(p(0,0,n(3)  ),1,halo(stencil_dir,idir),nb(1,idir),0, &
                        p(0,0,0     ),1,halo(stencil_dir,idir),nb(0,idir),0, &
                        comm_cart1,status,ierr)
    endif
    end select
    return
  end subroutine updthalo

subroutine structureFunction(istep, dir, u,v,w, ng,n,dims)
  implicit none
  real(rp),     intent(in), dimension(1:,1:,1:)                           :: u,v,w
  integer,      intent(in), dimension(2)                                  :: dims
  integer,      intent(in), dimension(3)                                  :: ng, n
  integer,      intent(in)                                                :: istep, dir
  real(rp),                 dimension(2*ng(dir))                          :: uExt, vExt, wExt
  real(rp)                                                                :: du,dv,dw
  real(rp),                 dimension(nmom_sf)                            :: Su, Sv, Sw
  integer                                                                 :: nhet, nhom1, nhom2
  integer                                                                 :: i1,i2,i3, il, ex
  character(len=1)                                                        :: idirchar,exstr
  character(len=7)                                                        :: fldnum

  if(myid.eq.0) then
    write(fldnum,'(i7.7)') istep
    write(idirchar,'(i1.1)') dir
    do ex=1,nmom_sf
      write(exstr,'(i1.1)') ex
      open(20+ex, file='data/post/sf/sf'//idirchar//'mom'//exstr//'fld'//fldnum//'.out')
      write(20+ex,'(A)') '# Structure function moment '//exstr
      write(20+ex,'(A)') '# Su^n, Sv^n, Sw^n'
    enddo
  endif
  select case(dir)
    case(1)
      nhom1 = n(2)
      nhom2 = n(3)
      nhet  = n(1)
    case(2)
      nhom1 = n(1)
      nhom2 = n(3)
      nhet  = n(2)
    case(3)
      nhom1 = n(1)
      nhom2 = n(2)
      nhet  = n(3)
  end select

  do il=1,nl_sf
    Su = 0.0
    Sv = 0.0
    Sw = 0.0
    do i1=1,nhom1
      do i2=1,nhom2
        select case(dir)
          case(1)
            uExt(1:nhet)  = u(:,i1,i2)
            uExt(nhet+1:) = u(:,i1,i2)
            vExt(1:nhet)  = v(:,i1,i2)
            vExt(nhet+1:) = v(:,i1,i2)
            wExt(1:nhet)  = w(:,i1,i2)
            wExt(nhet+1:) = w(:,i1,i2)
          case(2)
            uExt(1:nhet)  = u(i1,:,i2)
            uExt(nhet+1:) = u(i1,:,i2)
            vExt(1:nhet)  = v(i1,:,i2)
            vExt(nhet+1:) = v(i1,:,i2)
            wExt(1:nhet)  = w(i1,:,i2)
            wExt(nhet+1:) = w(i1,:,i2)
          case(3)
            uExt(1:nhet)  = u(i1,i2,:)
            uExt(nhet+1:) = u(i1,i2,:)
            vExt(1:nhet)  = v(i1,i2,:)
            vExt(nhet+1:) = v(i1,i2,:)
            wExt(1:nhet)  = w(i1,i2,:)
            wExt(nhet+1:) = w(i1,i2,:)
        end select
          do i3=1,nhet
            du = uExt(i3) - uExt(i3+il*dl_sf)
            dv = vExt(i3) - vExt(i3+il*dl_sf)
            dw = wExt(i3) - wExt(i3+il*dl_sf)
            do ex=1,nmom_sf
              Su(ex) = Su(ex) + du**ex
              Sv(ex) = Sv(ex) + dv**ex
              Sw(ex) = Sw(ex) + dw**ex
            enddo
          enddo
        enddo
      enddo
        do ex=1,nmom_sf
          Su(ex) = Su(ex)/product(n)
          Sv(ex) = Sv(ex)/product(n)
          Sw(ex) = Sw(ex)/product(n)
          call mpi_allreduce(MPI_IN_PLACE,Su(ex) ,1,mpi_real8,mpi_sum,comm_cart,ierr)
          call mpi_allreduce(MPI_IN_PLACE,Sv(ex) ,1,mpi_real8,mpi_sum,comm_cart,ierr)
          call mpi_allreduce(MPI_IN_PLACE,Sw(ex) ,1,mpi_real8,mpi_sum,comm_cart,ierr)
          Su(ex) = Su(ex)/product(dims)
          Sv(ex) = Sv(ex)/product(dims)
          Sw(ex) = Sw(ex)/product(dims)
          if (myid.eq.0) write(20+ex,'(3E15.7)') Su(ex), Sv(ex), Sw(ex)
        enddo
    enddo

    do ex=1,nmom_sf
      close(20+ex)
    enddo
    
end subroutine structureFunction

subroutine decoupledAutocorr(istep, dir, u,v,w,vof, ng,n,dims)
  implicit none
  real(rp),     intent(in), dimension(1:,1:,1:)                           :: u,v,w,vof
  integer,      intent(in), dimension(2)                                  :: dims
  integer,      intent(in), dimension(3)                                  :: ng, n
  integer,      intent(in)                                                :: istep, dir
  real(rp),                 dimension(2*ng(dir))                          :: uExt, vExt, wExt, vofExt
  real(rp)                                                                :: du,dv,dw
  real(rp),                 dimension(ng(dir))                            :: Ru1, Rv1, Rw1, Ru0, Rv0, Rw0
  integer                                                                 :: nhet, nhom1, nhom2
  integer                                                                 :: i1,i2,i3, dr, count1, count0
  character(len=1)                                                        :: idirchar
  character(len=7)                                                        :: fldnum

  select case(dir)
    case(1)
      nhom1 = n(2)
      nhom2 = n(3)
      nhet  = n(1)
    case(2)
      nhom1 = n(1)
      nhom2 = n(3)
      nhet  = n(2)
    case(3)
      nhom1 = n(1)
      nhom2 = n(2)
      nhet  = n(3)
  end select

  Ru1 = 0.0
  Rv1 = 0.0
  Rw1 = 0.0
  Ru0 = 0.0
  Rv0 = 0.0
  Rw0 = 0.0
  do dr=1,nhet
    count1 = 0
    count0 = 0
    do i1=1,nhom1
      do i2=1,nhom2
        select case(dir)
          case(1)
            uExt(1:nhet)  = u(:,i1,i2)
            uExt(nhet+1:) = u(:,i1,i2)
            vExt(1:nhet)  = v(:,i1,i2)
            vExt(nhet+1:) = v(:,i1,i2)
            wExt(1:nhet)  = w(:,i1,i2)
            wExt(nhet+1:) = w(:,i1,i2)
            vofExt(1:nhet)  = vof(:,i1,i2)
            vofExt(nhet+1:) = vof(:,i1,i2)
          case(2)
            uExt(1:nhet)  = u(i1,:,i2)
            uExt(nhet+1:) = u(i1,:,i2)
            vExt(1:nhet)  = v(i1,:,i2)
            vExt(nhet+1:) = v(i1,:,i2)
            wExt(1:nhet)  = w(i1,:,i2)
            wExt(nhet+1:) = w(i1,:,i2)
            vofExt(1:nhet)  = vof(i1,:,i2)
            vofExt(nhet+1:) = vof(i1,:,i2)
          case(3)
            uExt(1:nhet)  = u(i1,i2,:)
            uExt(nhet+1:) = u(i1,i2,:)
            vExt(1:nhet)  = v(i1,i2,:)
            vExt(nhet+1:) = v(i1,i2,:)
            wExt(1:nhet)  = w(i1,i2,:)
            wExt(nhet+1:) = w(i1,i2,:)
            vofExt(1:nhet)  = vof(i1,i2,:)
            vofExt(nhet+1:) = vof(i1,i2,:)
        end select
          do i3=1,nhet
            if ((vofExt(i3).ge.0.51) .and. (vofExt(i3+dr-1).ge.0.51))  then
              Ru1(dr) = Ru1(dr) + uExt(i3)*uExt(i3+dr-1)
              Rv1(dr) = Rv1(dr) + vExt(i3)*vExt(i3+dr-1)
              Rw1(dr) = Rw1(dr) + wExt(i3)*wExt(i3+dr-1)
              count1  = count1  + 1
            else if ((vofExt(i3).le.0.49) .and. (vofExt(i3+dr-1).le.0.49)) then
              Ru0(dr) = Ru0(dr) + uExt(i3)*uExt(i3+dr-1)
              Rv0(dr) = Rv0(dr) + vExt(i3)*vExt(i3+dr-1)
              Rw0(dr) = Rw0(dr) + wExt(i3)*wExt(i3+dr-1)
              count0  = count0  + 1
            endif
          enddo
      enddo
    enddo
    if (count1.eq.0) count1=1
    if (count0.eq.0) count0=1
    Ru1(dr) = Ru1(dr)/(1.0*count1)
    Rv1(dr) = Rv1(dr)/(1.0*count1)
    Rw1(dr) = Rw1(dr)/(1.0*count1)
    Ru0(dr) = Ru0(dr)/(1.0*count0)
    Rv0(dr) = Rv0(dr)/(1.0*count0)
    Rw0(dr) = Rw0(dr)/(1.0*count0)
  enddo
  ! Ru1 = Ru1/(nhom1*nhom2)
  ! Rv1 = Rv1/(nhom1*nhom2)
  ! Rw1 = Rw1/(nhom1*nhom2)
  ! Ru0 = Ru0/(nhom1*nhom2)
  ! Rv0 = Rv0/(nhom1*nhom2)
  ! Rw0 = Rw0/(nhom1*nhom2)
  call mpi_allreduce(MPI_IN_PLACE,Ru1 ,nhet,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,Rv1 ,nhet,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,Rw1 ,nhet,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,Ru0 ,nhet,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,Rv0 ,nhet,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,Rw0 ,nhet,mpi_real8,mpi_sum,comm_cart,ierr)
  Ru1 = Ru1/product(dims)
  Rv1 = Rv1/product(dims)
  Rw1 = Rw1/product(dims)
  Ru0 = Ru0/product(dims)
  Rv0 = Rv0/product(dims)
  Rw0 = Rw0/product(dims)

  if(myid.eq.0) then
    write(fldnum,'(i7.7)') istep
    write(idirchar,'(i1.1)') dir
    open(20, file='data/post/RR/RR'//idirchar//'fld'//fldnum//'.out')
    ! write(20,'(A,I1,A)') '# E'//index//'(k',idir,')'
    write(20,'(A)') '# Decoupled autocorrelation for vof=1 and vof=0'
    write(20,'(A)') '# Ru0, Ru1, Rv0, Rv1, Rw0, Rw1'

    do dr=1,nhet
      write(20,'(6E15.7)') Ru0(dr),Ru1(dr),Rv0(dr),Rv1(dr),Rw0(dr),Rw1(dr)
    enddo
  close(20)
  endif
end subroutine decoupledAutocorr
end module mod_spectra


