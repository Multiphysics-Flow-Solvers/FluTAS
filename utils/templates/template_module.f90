module mod_template
  !
  use mod_sub1, only: what_to_import 
  !
  implicit none 
  !
  private ! what is not public is private by default, no need to decleare.
  public  :: sub_a,sub_b
  !
  contains
  !
  subroutine sub_a(a1,a2,a3,nh_1,nh_2,a4,a5,a6,a7,a8)
    !
    ! Few lines of comments (please add them)
    !
    implicit none 
    !
    character(len=*), intent(in )                                     :: a1
    integer         , intent(in ), dimension(3)                       :: a2 
    integer         , intent(in )                                     :: a3 
    integer         , intent(in )                                     :: nh_1,nh_2
    real(rp)        , intent(in ), dimension(3)                       :: a4 
    logical         , intent(in )                                     :: a5 
    real(rp)        , intent(in ), dimension(1-nh_1:,1-nh_1:,1-nh_1:) :: a6 
    real(rp)        , intent(in ), dimension(1-nh_2:,1-nh_2:,1-nh_2:) :: a7 
    real(rp)        , intent(out), dimension(     0:,     0:,     0:) :: a8 
    !
    real(rp), dimension(0:,0:,0:) :: b1
    !
    ! Body of the subroutine
    !
    return
  end subroutine sub_a
  !
  subroutine sub_b(nx,ny,nz,a1,a2,a3,a4,nh_1,nh_2,a5,a6,a7)
    !
    ! Few lines of comments (please add them)
    !
    implicit none 
    !
    integer , intent(in )                                     :: nx,ny,nz
    integer , intent(in ), dimension(3)                       :: a1  
    real(rp), intent(in ), dimension(3)                       :: a2 
    integer , intent(in )                                     :: a3
    logical , intent(in )                                     :: a4 
    integer , intent(in )                                     :: nh_1,nh_2
    real(rp), intent(in ), dimension(1-nh_1:,1-nh_1:,1-nh_1:) :: a5 
    real(rp), intent(in ), dimension(1-nh_2:,1-nh_2:,1-nh_2:) :: a6 
    real(rp), intent(out), dimension(     0:,     0:,     0:) :: a7 
    !                                                             
    real(rp) :: b1 
    real(rp) :: b2 
    integer  :: b3 
    !
    ! Body of the subroutine
    !
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ! ...
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine sub_b
  !
end module mod_template
