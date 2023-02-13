!
! SPDX-License-Identifier: MIT
!
module mod_load
  !
  use mpi
  use mod_common_mpi, only: ierr,myid,ipencil
  use mod_param     , only: ng
  use mod_types     , only: rp
  use mod_sanity    , only: flutas_error
  use decomp_2d
  use decomp_2d_io
  !
  implicit none
  !
  private
  public  :: load,load_ns,load_scalar,cmpt_it_chkpt
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
    logical :: is_fld
    integer :: lenr,fh
    !
    select case(io)
    case('r')
      !
      ! check if the field <<fld>> exists
      !
      inquire(file=filename,exist=is_fld)
      if(.not.is_fld.and.myid.eq.0) then
        print*, ''
        print*, '*** The restarting field ', filename, 'does not exist! ***'
        print*, 'Please ensure that the field ', filename, 'has been printed in the correct restarting directory!'
        call flutas_error
      endif
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      !
      ! check its size 
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
      ! read it
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
    logical :: is_file
    integer :: fh
    !
    select case(io)
    case('r')
      !
      ! check if the file exists
      !
      inquire(file=filename,exist=is_file)
      if(.not.is_file.and.myid.eq.0) then
        print*, ''
        print*, '*** The restarting field ', filename, 'does not exist! ***'
        print*, 'Please ensure that the field ', filename, 'has been printed in the correct restarting directory!'
        call flutas_error
      endif
      !
      ! read it
      !
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
    ! reads/writes a restart file in the CaNS's output format (see: https://github.com/CaNS-World/CaNS)
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
    logical :: is_fld
    integer :: fh
    !
    select case(io)
    case('r')
      !
      ! check if the field <<fld>> exists
      !
      inquire(file=filename,exist=is_fld)
      if(.not.is_fld.and.myid.eq.0) then
        print*, ''
        print*, '*** The restarting field ', filename, 'does not exist! ***'
        print*, 'Please ensure that ', filename, 'has been printed in the correct restarting directory!'
        call flutas_error
      endif
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      !
      ! check its size
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
      ! read it
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
  subroutine cmpt_it_chkpt(restart_dir,it_chkpt)
    !
    ! determine the checkpoint from which to restart the simulation
    !
    implicit none
    !
    character(len=*), intent(in ) :: restart_dir
    integer         , intent(out) :: it_chkpt
    !
    real(rp), allocatable, dimension(:,:) :: mat_chkpt
    integer :: i,nrow
    logical :: is_file
    !
    ! check if the file "restart_checkpoints" exists 
    !
    inquire(file=trim(restart_dir)//'restart_checkpoints.out',exist=is_file)
    if(.not.is_file.and.myid.eq.0) call flutas_error('The file <<restart_checkpoints.out>> is absent. Please, provide it')
    if(.not.is_file.and.myid.eq.0) call flutas_error('In this case, the restarting files are probably missing. Check dns.in')
    !
    if(myid.eq.0) then ! only the first task
      !
      ! count the number line of the file "restart_checkpoints.out"
      !
      call execute_command_line('rm -f '//trim(restart_dir)//'line_of_restart.out') ! remove possible old files
      call execute_command_line('cat '//trim(restart_dir)//'restart_checkpoints.out | & 
                                 wc -l>> '//trim(restart_dir)//'line_of_restart.out')
      !
      ! get the number of line of the restarting file
      !
      open(10,file=trim(restart_dir)//'line_of_restart.out',status='old',action="read",iostat=ierr)
      read(10,*) nrow
      close(10)
      !
      ! go to the last line, third column, where the latest checkpoint is written
      !
      allocate(mat_chkpt(nrow,3))
      open(10,file=trim(restart_dir)//'restart_checkpoints.out',status='old',action='read',iostat=ierr)
      do i=1,nrow
        read(10,*) mat_chkpt(i,1),mat_chkpt(i,2),mat_chkpt(i,3)
      enddo
      close(10)
      it_chkpt = int(mat_chkpt(nrow,3))
      deallocate(mat_chkpt)
      !
    endif
    !
    ! broadcast "it_chkpt" to all the other processors
    !
    call MPI_BCAST(it_chkpt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    !
    return
  end subroutine cmpt_it_chkpt
  !
end module mod_load
