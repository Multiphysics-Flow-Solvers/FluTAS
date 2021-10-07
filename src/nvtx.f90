!
! Copyright (c) 2016, NVIDIA CORPORATION. All rights reserved.
!
!
! Permission is hereby granted, free of charge, to any person obtaining a
! copy of this software and associated documentation files (the "Software"),
! to deal in the Software without restriction, including without limitation
! the rights to use, copy, modify, merge, publish, distribute, sublicense,
! and/or sell copies of the Software, and to permit persons to whom the
! Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.
!

! Interface for NVTX markers
module nvtx
   use iso_c_binding
#if defined(_OPENACC)
   use cudafor
   implicit none

   integer, private :: col(7) = [Z'0000ff00', Z'000000ff', Z'00ffff00', Z'00ff00ff', Z'0000ffff', &
                                 Z'00ff0000', Z'00ffffff']
   character(len=256), private :: tempName

   logical, save :: isopen = .false.

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

#if defined(_USE_NVTX)
   interface nvtxRangePush
      ! push range with custom label and standard color
      subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
         use iso_c_binding
         character(kind=C_CHAR, len=*) :: name
      end subroutine nvtxRangePushA

      ! push range with custom label and custom color
      subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
         use iso_c_binding
         import :: nvtxEventAttributes
         type(nvtxEventAttributes) :: event
      end subroutine nvtxRangePushEx
   end interface nvtxRangePush

   interface nvtxRangePop
      subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
      end subroutine nvtxRangePop
   end interface nvtxRangePop
#endif
#endif

contains

   subroutine nvtxStartRange(name, id)
      character(kind=c_char, len=*) :: name
      integer, optional :: id
#if defined(_OPENACC)
#if defined(_USE_NVTX)
      type(nvtxEventAttributes) :: event
      integer :: istat
      istat = cudaDeviceSynchronize()

      tempName = trim(name)//c_null_char

      if (.not. present(id)) then
         call nvtxRangePush(tempName)
      else
         event%color = col(mod(id, 7) + 1)
         event%message = c_loc(tempName)
         call nvtxRangePushEx(event)
      end if
#endif
#endif
   end subroutine nvtxStartRange

   subroutine nvtxEndRange
#if defined(_OPENACC)
#if defined(_USE_NVTX)
      integer :: istat
      istat = cudaDeviceSynchronize()
      call nvtxRangePop
#endif
#endif
   end subroutine nvtxEndRange

   subroutine nvtxStartRangeAsync(name, id)
      character(kind=c_char, len=*) :: name
      integer, optional :: id
#if defined(_OPENACC)
#if defined(_USE_NVTX)
      type(nvtxEventAttributes) :: event
      tempName = trim(name)//c_null_char
      if (.not. present(id)) then
         call nvtxRangePush(tempName)
      else
         event%color = col(mod(id, 7) + 1)
         event%message = c_loc(tempName)
         call nvtxRangePushEx(event)
      end if
#endif
#endif
   end subroutine nvtxStartRangeAsync

   subroutine nvtxEndRangeAsync
#if defined(_OPENACC)
#if defined(_USE_NVTX)
      call nvtxRangePop
#endif
#endif
   end subroutine nvtxEndRangeAsync

end module nvtx
