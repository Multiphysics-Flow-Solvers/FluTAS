!
! SPDX-License-Identifier: MIT
!
module mod_load
  !
  use mpi
  use mod_common_mpi, only: ierr,myid,ipencil
  use mod_param     , only: ng
  use mod_types     , only: rp
  use decomp_2d
  use decomp_2d_io
  !
  implicit none
  !
  private
  public  :: load,load_ns,load_scalar
  !
  contains
  !
  subroutine load(io,filename,n,fld)
    !
    ! reads/writes a restart file
    !
    implicit none
    !
    character(len=1), intent(in   )                            :: io
    character(len=*), intent(in   )                            :: filename
    integer         , intent(in   ), dimension(3)              :: n
    real(rp)        , intent(inout), dimension(n(1),n(2),n(3)) :: fld ! generic field to be read/written
    !
    integer(MPI_OFFSET_KIND) :: filesize,disp,good
    integer :: lenr,fh
    !
    select case(io)
    case('r')
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      lenr = storage_size(fld(1,1,1))/8
      !good = int(ng(1)*ng(2)*ng(3),MPI_OFFSET_KIND)*lenr
      good = int(lenr,MPI_OFFSET_KIND)*ng(1)*ng(2)*ng(3)
      !
      if(filesize.ne.good) then
        if(myid.eq.0) print*, ''
        if(myid.eq.0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
        if(myid.eq.0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call decomp_2d_finalize
        call MPI_FINALIZE(ierr)
        call exit
      endif
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_read_var(fh,disp,ipencil,fld )
      call MPI_FILE_CLOSE(fh,ierr)
      ! 
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_write_var(fh,disp,ipencil,fld)
      call MPI_FILE_CLOSE(fh,ierr)
      !
    end select
    !
    return
  end subroutine load
  !
  subroutine load_scalar(io,filename,time,istep,dto)
    !
    implicit none
    !
    character(len=1), intent(in   ) :: io
    character(len=*), intent(in   ) :: filename
    integer         , intent(inout) :: istep
    real(rp)        , intent(inout) :: time, dto
    !
    integer :: fh
    !
    select case(io)
    case('r')
      open(88,file=filename,status='old',action='read')
      read(88,*) time, dto, istep
      close(88)
    case('w')
      if(myid.eq.0) then
        open (88, file=filename)
        write(88,'(2E15.7, 1I9.8)') time, dto, istep
        close(88)
      endif
    end select
    !
    return
  end subroutine
  !
  subroutine load_ns(io,filename,n,u,v,w,p,time,istep)
    !
    ! reads/writes a restart file
    !
    implicit none
    !
    character(len=1), intent(in   )                            :: io
    character(len=*), intent(in   )                            :: filename
    integer         , intent(in   ), dimension(3)              :: n
    real(rp)        , intent(inout), dimension(n(1),n(2),n(3)) :: u,v,w,p
    real(rp)        , intent(inout)                            :: time,istep
    !
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
    real(rp)  , dimension(2) :: fldinfo
    integer(8), dimension(3) :: ng
    integer(8) :: lenr
    integer    :: fh
    !
    select case(io)
    case('r')
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      !ng(:)   = n(:)
      !ng(1:3) = ng(1:3)*dims(1:3)
      lenr = storage_size(time)/8
      good = int(product(ng)*4+2,MPI_OFFSET_KIND)*lenr
      if(filesize.ne.good) then
        if(myid.eq.0) print*, ''
        if(myid.eq.0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
        if(myid.eq.0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call decomp_2d_finalize
        call MPI_FINALIZE(ierr)
        call exit
      endif
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_read_var(fh,disp,ipencil,u   )
      call decomp_2d_read_var(fh,disp,ipencil,v   )
      call decomp_2d_read_var(fh,disp,ipencil,w   )
      call decomp_2d_read_var(fh,disp,ipencil,p   )
      call decomp_2d_read_scalar(fh,disp,2,fldinfo)
      time  = fldinfo(1)
      istep = fldinfo(2)
      call MPI_FILE_CLOSE(fh,ierr)
      !
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_write_var(fh,disp,ipencil,u  )
      call decomp_2d_write_var(fh,disp,ipencil,v  )
      call decomp_2d_write_var(fh,disp,ipencil,w  )
      call decomp_2d_write_var(fh,disp,ipencil,p  )
      fldinfo = (/time,istep/)
      call decomp_2d_write_scalar(fh,disp,2,fldinfo)
      call MPI_FILE_CLOSE(fh,ierr)
      !
    end select
    !
    return
  end subroutine load_ns
  !
end module mod_load
