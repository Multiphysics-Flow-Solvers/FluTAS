!
! SPDX-License-Identifier: MIT
!
! ** PROFILER (v0.3) ** 
! This module combines the ability to track time using MPI_WTIME and trigger 
! NVTX. A small report is printed on stdout. It works under the assumption that
! there are always matching and consistent pairs of 'start' + 'stop'.
! I am going to extend it a bit further and create a standalone piece 
! of code that can be reused in other codes.
! -- Filippo Spiga
!
module profiler
  !
  use iso_fortran_env, only: real64, stdout => output_unit
  use mpi
  use iso_c_binding
  !
  implicit none
  save
  !
  integer, parameter :: max_clocks = 16
  integer, parameter :: labels_length_max = 16
  real(KIND=real64), parameter :: not_running = -1.0_real64
  !
  integer :: my_rank
  integer :: num_clocks = 0
  character(len=labels_length_max) :: labels(max_clocks)
  real(KIND=real64) :: stime(max_clocks)
  real(KIND=real64) :: ttime(max_clocks)
  integer :: how_many(max_clocks)
  logical :: tagged(max_clocks)
  !
  ! Color palette
  integer(C_INT), parameter :: COLOR_RED = int(z'00ff0000',4)
  integer(C_INT), parameter :: COLOR_GREEN = int(z'0000ff00',4) 
  integer(C_INT), parameter :: COLOR_BLUE = int(z'000000ff',4)
  integer(C_INT), parameter :: COLOR_YELLOW = int(z'00ffff00',4) 
  integer(C_INT), parameter :: COLOR_CYAN = int(z'0000ffff',4)
  integer(C_INT), parameter :: COLOR_PURPLE = int(z'00ff00ff',4)
  integer(C_INT), parameter :: COLOR_WHITE = int(z'00ffffff',4)
  !
#if defined(_USE_NVTX)
  character(len=256), private :: tempName
  !
  type, bind(C) :: nvtxEventAttributes
    integer(C_INT16_T) :: version = 1
    integer(C_INT16_T) :: size = 48 !
    integer(C_INT) :: category = 0
    integer(C_INT) :: colorType = 1 ! NVTX_COLOR_ARGB = 1
    integer(C_INT) :: color
    integer(C_INT) :: payloadType = 0 ! NVTX_PAYLOAD_UNKNOWN = 0
    integer(C_INT) :: reserved0
    integer(C_INT64_T) :: payload   ! union uint,int,double
    integer(C_INT) :: messageType = 1  ! NVTX_MESSAGE_TYPE_ASCII     = 1
    type(C_PTR) :: message  ! ascii char
  end type nvtxEventAttributes
  !
  interface nvtxRangePush
    ! push range with custom label and standard color
    subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
        use iso_c_binding
        character(kind=C_CHAR, len=*) :: name
    end subroutine nvtxRangePushA
    !
    ! push range with custom label and custom color
    subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
        use iso_c_binding
        import :: nvtxEventAttributes
        type(nvtxEventAttributes) :: event
    end subroutine nvtxRangePushEx
  end interface nvtxRangePush
  !
  interface nvtxRangePop
    subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
    end subroutine nvtxRangePop
  end interface nvtxRangePop
#endif
  !
  public :: profiler_init, profiler_start, profiler_stop, profiler_report
  !
contains
  !
  subroutine profiler_init()
    implicit none  
    !
    integer :: id, ierr
    !
    num_clocks = 0
    !
    do id = 1, max_clocks
      how_many(id) = 0
      ttime(id) = 0.0_real64
      stime(id) = not_running
      labels(id) = ' '
      tagged(id) = .false.
    end do
    !
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    !
    return
  end subroutine profiler_init
  !
  subroutine profiler_start(label, tag, tag_color)
    implicit none
    !
    character(len=*), intent(in) :: label
    logical, intent(in), optional :: tag
    integer, intent(in), optional :: tag_color
    !
    character(len=labels_length_max) :: label_
    integer :: id
    !
    label_ = trim(label)
    !
    do id = 1, num_clocks
      !
      if (labels(id) == label_) then
        !
        stime(id) = MPI_WTIME()
        !
        if (tagged(id)) then
#if defined(_USE_NVTX)
          if (present(tag_color)) then
            call nvtxStartRange(label_, tag_color)
          else
            call nvtxStartRange(label_)
          end if
#endif
        end if
        !
        ! Early exit
        return
      end if
      !
    end do
    !
    ! First time clock is encountered, add it in the pool
    num_clocks = num_clocks + 1
    labels(num_clocks) = label_
    stime(num_clocks) = MPI_WTIME()
    !
    if (present(tag)) then
      if (tag) then
        tagged(id) = .true.
#if defined(_USE_NVTX)
        if (present(tag_color)) then
          call nvtxStartRange(label_, tag_color)
        else
          call nvtxStartRange(label_)
        end if
#endif
      end if
    end if
    !
    return
  end subroutine profiler_start
  !
  subroutine profiler_stop(label)
    implicit none
    !
    character(len=*), intent(in) :: label
    character(len=labels_length_max) :: label_
    integer :: id
    !
    label_ = trim(label)
    !
    do id = 1, num_clocks
      !
      if (labels(id) == label_) then
        !
        ttime(id) = ttime(id) + (MPI_WTIME() - stime(id))
        stime(id) = not_running
        how_many(id) = how_many(id) + 1
        if (tagged(id)) then
#if defined(_USE_NVTX)
          call nvtxEndRange
#endif
        end if
        !
        ! Early exit
        return
      end if
      !
    end do
    !
    return
  end subroutine profiler_stop
  !
  subroutine profiler_report()
    implicit none
    !
    character(len=labels_length_max) :: label_
    integer :: id, ierr
    !
    ! Print all clocks
    write (stdout, *)
    !
    do id = 1, num_clocks
      !
      if (my_rank .eq. 0) write (stdout, '(1X,A12," : ",F10.3,"s ( ", I8," calls)")') labels(id), ttime(id), how_many(id)
      !
    end do
    !
    return
  end subroutine profiler_report
  !
#if defined(_USE_NVTX)
  ! 
  subroutine nvtxStartRange(name, id)
    character(kind=c_char, len=*) :: name
    integer, optional :: id
    type(nvtxEventAttributes) :: event
    integer :: istat
    !
    tempName = trim(name)//c_null_char
    !
    if (.not. present(id)) then
        call nvtxRangePush(tempName)
    else
        event%color = id
        event%message = c_loc(tempName)
        call nvtxRangePushEx(event)
    end if
  end subroutine nvtxStartRange
  !
  subroutine nvtxEndRange
    integer :: istat
    call nvtxRangePop
  end subroutine nvtxEndRange
  !
  subroutine nvtxStartRangeAsync(name, id)
    character(kind=c_char, len=*) :: name
    integer, optional :: id
    type(nvtxEventAttributes) :: event
    !
    tempName = trim(name)//c_null_char
    if (present(id)) then
        event%color = id
    end if
    !
    event%message = c_loc(tempName)
  end subroutine nvtxStartRangeAsync
  !
  subroutine nvtxEndRangeAsync
    call nvtxRangePop
  end subroutine nvtxEndRangeAsync
#endif
  !
end module profiler
