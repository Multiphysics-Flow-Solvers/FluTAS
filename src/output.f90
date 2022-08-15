!
! SPDX-License-Identifier: MIT
!
module mod_output
  !
  ! note: this output module serves for including generic input/output 
  !       subroutines that are usefull for all the different applications:
  !          --> plot one observable (out0d);
  !          --> plot one observable over a line (out1d);
  !          --> plot one observable over a plane (out2d);
  !          --> plot one observable over the entire computational domain (out3d);
  !       
  ! We recommend to include the postprocessing routines in the following location:
  !   --> If the routine is application-dependent include in apps/<YOUR_APPLICATION>/postp.<YOUR_APPLICATION>
  !   --> If it can serve many applications include in postprocessing/post.f90
  !
  ! Remember to add a comment on what the subroutine is doing
  !
  use mpi
  use decomp_2d_io
  use mod_common_mpi, only: ierr,myid,ipencil,ijk_start
  use mod_types
  !@cuf use mod_common_mpi, only: mydev
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: out0d,out1d,out2d,out3d,write_visu_2d,write_visu_3d
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
    ! n     -> number of points in each computational subdomain
    ! dims  -> number of computational subdomains along the parallelized directions
    ! dl    -> spacing along x, y and z
    ! idir  -> direction along which we print the profile (the other two set are averaged directions)
    ! nh_d  -> number of halo points for the non uniform spacing along z
    ! nh_p  -> number of halo points for the generic 3D input scalar field
    ! z     -> global z coordinate (grid is non-uniform in z)
    ! dzlzi -> dz/lz weight of a grid cell for averaging over z (it is sufficient to provide the local spacing)
    ! p     -> generic 3D input scalar field which we want to average
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
    integer :: ngx,ngy,ngz,nx,ny,nz,start,ng_idir
    integer :: i,j,k,ii,jj,kk,mm
    integer :: iunit
    !@cuf integer :: istat
    !@cuf attributes(managed) :: z,dzlzi,p 
    !
    iunit   = 10
    ngx     = n(1)*dims(1)
    ngy     = n(2)*dims(2)
    ngz     = n(3)*dims(3)
    ng_idir = n(idir)*dims(idir)
    nx      = n(1)
    ny      = n(2)
    nz      = n(3)
    start   = ijk_start(idir)
    !
    allocate(p1d(ng_idir))
    !
    ! first, set to zero
    !
    !$acc kernels
    do mm=1,ng_idir
      p1d(mm) = 0._rp
    enddo
    !$acc end kernels 
    !
    select case(idir)
    case(3)
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            mm = start + k
            !
            p1d(mm) = p1d(mm) + p(i,j,k)
            !
          enddo
        enddo
      enddo
      !$acc end kernels 
      call mpi_allreduce(MPI_IN_PLACE,p1d(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do kk=1,ngz
          write(iunit,'(2E15.7)') z(kk),p1d(kk)/(1._rp*ngx*ngy)
        enddo
        close(iunit)
      endif
    case(2)
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            mm = start + j
            !
            p1d(mm) = p1d(mm) + p(i,j,k)*dzlzi(k)
            !
          enddo
        enddo
      enddo
      !$acc end kernels
      call mpi_allreduce(MPI_IN_PLACE,p1d(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do jj=1,ngy
          write(iunit,'(2E15.7)') (jj-0.5_rp)*dl(2),p1d(jj)/ngy
        enddo
        close(iunit)
      endif
    case(1)
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            mm = start + i
            !
            p1d(mm) = p1d(mm) + p(i,j,k)*dzlzi(k)
            !
          enddo
        enddo
      enddo
      !$acc end kernels 
      call mpi_allreduce(MPI_IN_PLACE,p1d(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do ii=1,ngx
          write(iunit,'(2E15.7)') (ii-0.5_rp)*dl(1),p1d(ii)/ngy
        enddo
        close(iunit)
      endif
    end select
    !
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
  subroutine write_log_output(fname,fname_fld,varname,nmin,nmax,nskip,time,istep)
    !
    ! appends information about a saved binary file to a log file
    ! this file is used to generate a xdmf file for visualization of field data
    !
    ! fname     -> name of the output log file
    ! fname_fld -> name of the saved binary file (excluding the directory)
    ! varname   -> name of the variable that is saved
    ! nmin      -> first element of the field that is saved in each direction, e.g. [1,1,1]
    ! nmax      -> last  element of the field that is saved in each direction, e.g. [ng(1),ng(2),ng(3)]
    ! nskip     -> step size between nmin and nmax, e.g. [1,1,1] if the whole array is saved
    ! time      -> physical time
    ! istep     -> time step number
    !
    implicit none
    !
    character(len=*), intent(in)               :: fname,fname_fld,varname
    integer         , intent(in), dimension(3) :: nmin,nmax,nskip
    real(rp)        , intent(in)               :: time
    integer         , intent(in)               :: istep
    !
    character(len=100) :: cfmt
    integer :: iunit
    !
    write(cfmt, '(A)') '(A,A,A,9i5,E15.7,i7)'
    if(myid  ==  0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,trim(cfmt)) trim(fname_fld),' ',trim(varname),nmin,nmax,nskip,time,istep
      close(iunit)
    endif
    !
    return
  end subroutine write_log_output
  !
  subroutine write_visu_3d(datadir,fname_bin,fname_log,varname,nmin,nmax,nskip,time,istep,p)
    !
    ! wraps the calls of out3d and write_log_output into the same subroutine
    !
    implicit none
    !
    character(len=*), intent(in)                   :: datadir,fname_bin,fname_log,varname
    integer         , intent(in), dimension(3)     :: nmin,nmax,nskip
    real(rp)        , intent(in)                   :: time
    integer         , intent(in)                   :: istep
    real(rp)        , intent(in), dimension(:,:,:) :: p
    !
    call out3d(trim(datadir)//trim(fname_bin),nskip,p)
    call write_log_output(trim(datadir)//trim(fname_log),trim(fname_bin),trim(varname),nmin,nmax,nskip,time,istep)
    !
    return
  end subroutine write_visu_3d
  !
  subroutine write_visu_2d(datadir,fname_bin,fname_log,varname,inorm,nslice,ng,time,istep,p)
    !
    ! wraps the calls of out2d and write-log_output into the same subroutine
    !
    implicit none
    !
    character(len=*), intent(in)                   :: datadir,fname_bin,fname_log,varname
    integer         , intent(in)                   :: inorm,nslice
    integer         , intent(in), dimension(3)     :: ng
    real(rp)        , intent(in)                   :: time
    integer         , intent(in)                   :: istep
    real(rp)        , intent(in), dimension(:,:,:) :: p
    integer         , dimension(3)                 :: nmin_2d,nmax_2d
    !
    call out2d(trim(datadir)//trim(fname_bin),inorm,nslice,p)
    select case(inorm)
    case(1)
      nmin_2d(:) = [nslice,1    ,1    ]
      nmax_2d(:) = [nslice,ng(2),ng(3)]
    case(2)
      nmin_2d(:) = [1    ,nslice,1    ]
      nmax_2d(:) = [ng(1),nslice,ng(3)]
    case(3)
      nmin_2d(:) = [1    ,1    ,nslice]
      nmax_2d(:) = [ng(1),ng(2),nslice]
    end select
    call write_log_output(trim(datadir)//trim(fname_log),trim(fname_bin),trim(varname),nmin_2d,nmax_2d,[1,1,1],time,istep)
    !
    return
  end subroutine write_visu_2d
  !
end module mod_output
