!
! SPDX-License-Identifier: MIT
!
module mod_output
  !
  use mpi
  use decomp_2d_io
  use mod_common_mpi, only: ierr,myid,ipencil,ijk_start
  use mod_types
  !
  implicit none
  !
  private
  public  :: out0d,out1d,out1d_2,out2d,out3d
#if defined(_TURB_FORCING)
  public  :: budget
#endif
  !
  contains
  !
  subroutine out0d(fname,n,var)
    !
    ! appends the first n entries of an array
    ! var to a file
    ! fname -> name of the file
    ! n     -> number of entries
    ! var   -> input array of real values
    !
    implicit none
    !
    character(len=*), intent(in)               :: fname
    integer         , intent(in)               :: n
    real(rp)        , intent(in), dimension(:) :: var
    !
    character(len=30) :: cfmt
    integer :: iunit
    integer :: i
    !
    write(cfmt,'(A,I3,A)') '(',n,'E15.7)'
    iunit = 10
    if(myid.eq.0) then
      open(iunit,file=fname,position='append')
      write(iunit,trim(cfmt)) (var(i),i=1,n) 
      close(iunit)
    endif
    !
    return
  end subroutine out0d
  !
  subroutine out1d(fname,n,dims,dl,idir,nh_d,nh_p,z,dzlzi,p)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions
    !
    ! fname -> name of the file
    ! n     -> size of the input array
    ! idir  -> direction of the profile
    ! z     -> z coordinate (grid is non-uniform in z)
    ! dzlzi -> dz/lz weight of a grid cell for averaging over z
    ! p     -> 3D input scalar field
    !
    implicit none
    !
    character(len=*), intent(in)                                     :: fname
    integer         , intent(in), dimension(3)                       :: n,dims
    real(rp)        , intent(in), dimension(3)                       :: dl
    integer         , intent(in)                                     :: idir
    integer         , intent(in)                                     :: nh_d,nh_p
    real(rp)        , intent(in), dimension(1-nh_d:)                 :: z,dzlzi
    real(rp)        , intent(in), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    !
    real(rp), allocatable, dimension(:) :: p1d
    integer, dimension(3) :: ng
    integer :: i,j,k,ii,jj,kk
    integer :: iunit
    !
    ng(:) = n(:)*dims(:)
    iunit = 10
    select case(idir)
    case(3)
      allocate(p1d(ng(3)))
      p1d(:) = 0._rp
      do k=1,n(3)
        kk = ijk_start(3) + k
        p1d(kk) = 0._rp
        do j=1,n(2)
          do i=1,n(1)
            p1d(kk) = p1d(kk) + p(i,j,k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,p1d(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d(:) = p1d(:)/(1.*ng(1)*ng(2))
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do kk=1,ng(3)
          write(iunit,'(2E15.7)') z(kk),p1d(kk)
        enddo
        close(iunit)
      endif
    case(2)
      allocate(p1d(ng(2)))
      p1d(:) = 0._rp
      do j=1,n(2)
        jj = ijk_start(2) + j
        p1d(jj) = 0._rp
        do k=1,n(3)
          do i=1,n(1)
            p1d(jj) = p1d(jj) + p(i,j,k)*dzlzi(k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,p1d(1),ng(2),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d(:) = p1d(:)/(1.*ng(1))
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do jj=1,ng(2)
          write(iunit,'(2E15.7)') (jj-0.5_rp)*dl(2),p1d(jj)
        enddo
        close(iunit)
      endif
    case(1)
      allocate(p1d(ng(1)))
      p1d(:) = 0._rp
      do i=1,n(1)
        ii = ijk_start(1) + i
        p1d(ii) = 0._rp
        do k=1,n(3)
          do j=1,n(2)
            p1d(ii) = p1d(ii) + p(i,j,k)*dzlzi(k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,p1d(1),ng(1),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d(:) = p1d(:)/(1.*ng(2))
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do ii=1,ng(1)
          write(iunit,'(2E15.7)') (ii-0.5_rp)*dl(1),p1d(ii)
        enddo
        close(iunit)
      endif
    end select
    deallocate(p1d)
    !
    return
  end subroutine out1d
  !
  subroutine out2d(fname,inorm,islice,p)
    !
    ! saves a planar slice of a scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! inorm  -> plane is perpendicular to direction
    !           inorm (1,2,3)
    ! islice -> plane is of constant index islice 
    !           in direction inorm
    ! p      -> 3D input scalar field
    !
    implicit none
    !
    character(len=*), intent(in)                   :: fname
    integer         , intent(in)                   :: inorm,islice
    real(rp)        , intent(in), dimension(:,:,:) :: p
    !
    select case(inorm)
    case(1) !normal to x --> yz plane
      call decomp_2d_write_plane(ipencil,p,inorm,islice,fname)
    case(2) !normal to y --> zx plane
      call decomp_2d_write_plane(ipencil,p,inorm,islice,fname)
    case(3) !normal to z --> xy plane
      call decomp_2d_write_plane(ipencil,p,inorm,islice,fname)
    end select
    !
    return
  end subroutine out2d
  !
  subroutine out3d(fname,nskip,p)
    !
    ! saves a 3D scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! nskip  -> array with the step size for which the
    !           field is written; i.e.: (/1,1,1/) 
    !           writes the full field 
    ! p      -> 3D input scalar field
    !
    implicit none
    !
    character(len=*), intent(in)                   :: fname
    integer         , intent(in), dimension(3)     :: nskip
    real(rp)        , intent(in), dimension(:,:,:) :: p
    !
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    integer :: fh
    !
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fname, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_write_every(ipencil,p,nskip(1),nskip(2),nskip(3),fname,.true.)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine out3d
  !
  subroutine out1d_2(fname,n,dims,idir,nh_d,nh_u,z,u,v,w) ! e.g. for a channel with streamwise dir in x
    !
    implicit none
    !
    character(len=*), intent(in)                                     :: fname
    integer         , intent(in), dimension(3)                       :: n,dims
    integer         , intent(in)                                     :: idir
    integer         , intent(in)                                     :: nh_d,nh_u
    real(rp)        , intent(in), dimension(1-nh_d:)                 :: z
    real(rp)        , intent(in), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp), allocatable, dimension(:) :: um,vm,wm,u2,v2,w2,uw
    integer, dimension(3) :: ng
    integer :: i,j,k,kk
    integer :: iunit
    integer :: q
    !
    ng(:) = n(:)*dims(:)
    iunit = 10
    select case(idir)
    case(3)
      q = ng(3)
      allocate(um(0:q+1),vm(0:q+1),wm(0:q+1),u2(0:q+1),v2(0:q+1),w2(0:q+1),uw(0:q+1))
      um(:) = 0._rp
      vm(:) = 0._rp
      wm(:) = 0._rp
      u2(:) = 0._rp
      v2(:) = 0._rp
      w2(:) = 0._rp
      uw(:) = 0._rp
      do k=1,n(3)
        kk = ijk_start(3) + k
        um(kk) = 0._rp
        vm(kk) = 0._rp
        wm(kk) = 0._rp
        u2(kk) = 0._rp
        v2(kk) = 0._rp
        w2(kk) = 0._rp
        uw(kk) = 0._rp
        do j=1,n(2)
          do i=1,n(1)
            um(kk) = um(kk) + u(i,j,k)
            vm(kk) = vm(kk) + v(i,j,k)
            wm(kk) = wm(kk) + 0.50_rp*(w(i,j,k-1) + w(i,j,k))
            u2(kk) = u2(kk) + u(i,j,k)**2
            v2(kk) = v2(kk) + v(i,j,k)**2
            w2(kk) = w2(kk) + 0.50_rp*(w(i,j,k)**2+w(i,j,k-1)**2)
            uw(kk) = uw(kk) + 0.25_rp*(u(i-1,j,k) + u(i,j,k))* &
                                   (w(i,j,k-1) + w(i,j,k))
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,um(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,wm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,u2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,v2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,w2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,uw(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:) = um(:)/(1.*ng(1)*ng(2))
      vm(:) = vm(:)/(1.*ng(1)*ng(2))
      wm(:) = wm(:)/(1.*ng(1)*ng(2))
      u2(:) = sqrt(u2(:)/(1.*ng(1)*ng(2)) - um(:)**2)
      v2(:) = sqrt(v2(:)/(1.*ng(1)*ng(2)) - vm(:)**2)
      w2(:) = sqrt(w2(:)/(1.*ng(1)*ng(2)) - wm(:)**2)
      uw(:) = uw(:)/(1.*ng(1)*ng(2)) - um(:)*wm(:)
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do kk=1,ng(3)
          write(iunit,'(8E15.7)') z(kk),um(kk),vm(kk),wm(kk), &
                                        u2(kk),v2(kk),w2(kk), &
                                        uw(kk)
        enddo
        close(iunit)
      endif
      deallocate(um,vm,wm,u2,v2,w2,uw)
    case(2)
    case(1)
    end select
    !
    return
  end subroutine out1d_2
  !
  subroutine out2d_2(fname,n,dims,dl,idir,nh_d,nh_u,z,u,v,w) ! e.g. for a duct
    !
    implicit none
    !
    character(len=*), intent(in)                                     :: fname
    integer         , intent(in), dimension(3)                       :: n,dims
    real(rp)        , intent(in), dimension(3      )                 :: dl
    integer         , intent(in)                                     :: idir
    integer         , intent(in)                                     :: nh_d,nh_u
    real(rp)        , intent(in), dimension(1-nh_d:)                 :: z
    real(rp)        , intent(in), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp), allocatable, dimension(:,:) :: um,vm,wm,u2,v2,w2,uv,uw,vw
    real(rp) :: x,y
    integer, dimension(3) :: ng
    integer  :: i,j,k,ii,jj,kk
    integer  :: iunit
    integer  :: p,q
    !
    ng(:) = n(:)*dims(:)
    iunit = 10
    select case(idir) ! streamwise direction
    case(3)
    case(2)
      p = ng(1)
      q = ng(3)
      allocate(um(p,q),vm(p,q),wm(p,q),u2(p,q),v2(p,q),w2(p,q),uv(p,q),vw(p,q))
      !
      um(:,:) = 0._rp
      vm(:,:) = 0._rp
      wm(:,:) = 0._rp
      u2(:,:) = 0._rp
      v2(:,:) = 0._rp
      w2(:,:) = 0._rp
      uv(:,:) = 0._rp
      vw(:,:) = 0._rp
      do k=1,n(3)
        kk = ijk_start(3) + k
        do i=1,n(1)
          ii = ijk_start(1) + i
          um(ii,kk) = 0._rp
          vm(ii,kk) = 0._rp
          wm(ii,kk) = 0._rp
          u2(ii,kk) = 0._rp
          v2(ii,kk) = 0._rp
          w2(ii,kk) = 0._rp
          vw(ii,kk) = 0._rp
          uv(ii,kk) = 0._rp
          do j=1,n(2)
            um(ii,kk) = um(ii,kk) + 0.50_rp*(u(i-1,j,k)+u(i,j,k))
            vm(ii,kk) = vm(ii,kk) + v(i,j,k)
            wm(ii,kk) = wm(ii,kk) + 0.50_rp*(w(i,j,k-1)+w(i,j,k))
            u2(ii,kk) = u2(ii,kk) + 0.50_rp*(u(i-1,j,k)**2+u(i,j,k)**2)
            v2(ii,kk) = v2(ii,kk) + v(i,j,k)**2
            w2(ii,kk) = w2(ii,kk) + 0.50_rp*(w(i,j,k-1)**2+w(i,j,k)**2)
            vw(ii,kk) = vw(ii,kk) + 0.25_rp*(v(i,j-1,k) + v(i,j,k))* &
                                            (w(i,j,k-1) + w(i,j,k))
            uv(ii,kk) = uv(ii,kk) + 0.25_rp*(u(i-1,j,k) + u(i,j,k))* &
                                            (v(i,j-1,k) + v(i,j,k))
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,um(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,wm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,u2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,v2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,w2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vw(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,uv(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:,:) =      um(:,:)/(1.*ng(2))
      vm(:,:) =      vm(:,:)/(1.*ng(2))
      wm(:,:) =      wm(:,:)/(1.*ng(2))
      u2(:,:) = sqrt(u2(:,:)/(1.*ng(2)) - um(:,:)**2)
      v2(:,:) = sqrt(v2(:,:)/(1.*ng(2)) - vm(:,:)**2)
      w2(:,:) = sqrt(w2(:,:)/(1.*ng(2)) - wm(:,:)**2)
      vw(:,:) =      vw(:,:)/(1.*ng(2)) - vm(:,:)*wm(:,:)
      uv(:,:) =      uv(:,:)/(1.*ng(2)) - um(:,:)*vm(:,:)
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do kk=1,ng(3)
          do ii=1,ng(1)
            x = (ii-0.50_rp)*dl(1)
            write(iunit,'(10E15.7)') x,z(kk),um(ii,kk),vm(ii,kk),wm(ii,kk), &
                                             u2(ii,kk),v2(ii,kk),w2(ii,kk), &
                                             vw(ii,kk),uv(ii,kk)
          enddo
        enddo
        close(iunit)
      endif
      deallocate(um,vm,wm,u2,v2,w2,vw,uv)
    case(1)
      p = ng(2)
      q = ng(3)
      allocate(um(p,q),vm(p,q),wm(p,q),u2(p,q),v2(p,q),w2(p,q),uv(p,q),uw(p,q))
      !
      um(:,:) = 0._rp
      vm(:,:) = 0._rp
      wm(:,:) = 0._rp
      u2(:,:) = 0._rp
      v2(:,:) = 0._rp
      w2(:,:) = 0._rp
      uv(:,:) = 0._rp
      uw(:,:) = 0._rp
      do k=1,n(3)
        kk = ijk_start(3) + k
        do j=1,n(2)
          jj = ijk_start(2) + j
          um(jj,kk) = 0._rp
          vm(jj,kk) = 0._rp
          wm(jj,kk) = 0._rp
          u2(jj,kk) = 0._rp
          v2(jj,kk) = 0._rp
          w2(jj,kk) = 0._rp
          uv(jj,kk) = 0._rp
          uw(jj,kk) = 0._rp
          do i=1,n(1)
            um(jj,kk) = um(jj,kk) + u(i,j,k)
            vm(jj,kk) = vm(jj,kk) + 0.50_rp*(v(i,j-1,k)+v(i,j,k))
            wm(jj,kk) = wm(jj,kk) + 0.50_rp*(w(i,j,k-1)+w(i,j,k))
            u2(jj,kk) = u2(jj,kk) + u(i,j,k)**2
            v2(jj,kk) = v2(jj,kk) + 0.50_rp*(v(i,j-1,k)**2+v(i,j,k)**2)
            w2(jj,kk) = w2(jj,kk) + 0.50_rp*(w(i,j,k-1)**2+w(i,j,k)**2)
            uv(jj,kk) = uv(jj,kk) + 0.25_rp*(u(i-1,j,k) + u(i,j,k))* &
                                            (v(i,j-1,k) + v(i,j,k))
            uw(jj,kk) = uw(jj,kk) + 0.25_rp*(u(i-1,j,k) + u(i,j,k))* &
                                            (w(i,j,k-1) + w(i,j,k))
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,um(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,wm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,u2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,v2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,w2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,uv(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,uw(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:,:) =      um(:,:)/(1.*ng(2))
      vm(:,:) =      vm(:,:)/(1.*ng(2))
      wm(:,:) =      wm(:,:)/(1.*ng(2))
      u2(:,:) = sqrt(u2(:,:)/(1.*ng(2)) - um(:,:)**2)
      v2(:,:) = sqrt(v2(:,:)/(1.*ng(2)) - vm(:,:)**2)
      w2(:,:) = sqrt(w2(:,:)/(1.*ng(2)) - wm(:,:)**2)
      uv(:,:) =      uv(:,:)/(1.*ng(2)) - um(:,:)*vm(:,:)
      uw(:,:) =      uw(:,:)/(1.*ng(2)) - um(:,:)*wm(:,:)
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do kk=1,ng(3)
          do jj=1,ng(2)
            y = (jj-0.50_rp)*dl(2)
            write(iunit,'(10E15.7)') y,z(kk),um(jj,kk),vm(jj,kk),wm(jj,kk), &
                                             u2(jj,kk),v2(jj,kk),w2(jj,kk), &
                                             uv(jj,kk),uw(jj,kk)
          enddo
        enddo
        close(iunit)
      endif
      deallocate(um,vm,wm,u2,v2,w2,vw,uv)
    end select
    !
    return
  end subroutine out2d_2
  !
  subroutine budget(datadir,nx,ny,nz,ngx,ngy,ngz,rho1,rho2,dx,dy,dz, &
                    nh_u,u,v,w,rho,vof,time,istep)
    !
    implicit none
    !
    character(len=100), intent(in )                                     :: datadir
    integer           , intent(in )                                     :: nx ,ny ,nz
    integer           , intent(in )                                     :: ngx,ngy,ngz
    real(rp)          , intent(in )                                     :: rho1,rho2
    real(rp)          , intent(in )                                     :: dx,dy,dz
    integer           , intent(in )                                     :: nh_u
    real(rp)          , intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp)          , intent(in ), dimension(     0:,     0:,     0:) :: rho
    real(rp)          , intent(in ), dimension(     0:,     0:,     0:) :: vof
    real(rp)          , intent(in )                                     :: time
    integer           , intent(in )                                     :: istep
    !
    real(rp) :: ke_t,ke_1,ke_2,dvol,volt,vol_t1,vol_t2
    integer  :: i,j,k,ip,jp,kp,im,jm,km
    !
    real(rp), parameter :: eps = real(1e-12,rp)
    ! 
    volt = 1._rp*ngx*ngy*ngz
    dvol = dx*dy*dz
    !
    vol_t1 = 0._rp
    vol_t2 = 0._rp
    !
    ke_t = 0._rp
    ke_1 = 0._rp
    ke_2 = 0._rp
    !
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
          ! 1. volume
          ! 
          vol_t1 = vol_t1 + vof(i,j,k)
          vol_t2 = vol_t2 + (1._rp-vof(i,j,k))
          !
          ! 2. kinetic energy
          ! 
          ke_t = ke_t + rho(i,j,k)*&
          0.5_rp*(0.25_rp*(u(i,j,k)+u(im,j,k))**2 + 0.25_rp*(v(i,j,k)+v(i,jm,k))**2 + 0.25_rp*(w(i,j,k)+w(i,j,km))**2)
          ke_1 = ke_1 + rho1*vof(i,j,k)*&
          0.5_rp*(0.25_rp*(u(i,j,k)+u(im,j,k))**2 + 0.25_rp*(v(i,j,k)+v(i,jm,k))**2 + 0.25_rp*(w(i,j,k)+w(i,j,km))**2)
          ke_2 = ke_2 + rho2*(1._rp-vof(i,j,k))*&
          0.5_rp*(0.25_rp*(u(i,j,k)+u(im,j,k))**2 + 0.25_rp*(v(i,j,k)+v(i,jm,k))**2 + 0.25_rp*(w(i,j,k)+w(i,j,km))**2)
          !
        enddo
      enddo
    enddo
    !
    call mpi_allreduce(MPI_IN_PLACE,vol_t1,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,vol_t2,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    call mpi_allreduce(MPI_IN_PLACE,ke_t  ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,ke_1  ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,ke_2  ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    ke_t = ke_t/volt
    ke_1 = ke_1/(vol_t1+eps)
    ke_2 = ke_2/(vol_t2+eps)
    !
    if(myid.eq.0) then
      !
      ! a. post-processing one fluid
      !
      open(92,file=trim(datadir)//'ke_t.out',position='append')
      write(92,'(3E15.7)') 1._rp*istep,time,ke_t
      close(92)
      ! 
      ! b1. phase 1
      ! 
      open(93,file=trim(datadir)//'ke_1.out',position='append')
      write(93,'(4E15.7)') 1._rp*istep,time,vol_t1,ke_1
      close(93)
      ! 
      ! b2. phase 2
      !
      open(94,file=trim(datadir)//'ke_2.out',position='append')
      write(94,'(4E15.7)') 1._rp*istep,time,vol_t2,ke_2 
      close(94)
      !
    endif
    !
    return
  end subroutine budget
  !
end module mod_output
