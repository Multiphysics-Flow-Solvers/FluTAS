module mod_tagging
  use mod_common_mpi, only: myid, comm_Cart,ierr, coord,status, &
  xsl_buf, xrl_buf, ysl_buf, yrl_buf, xsr_buf, xrr_buf, ysr_buf, yrr_buf, &
  left,right,front,back,xhalo,yhalo
  use mod_types
  use mpi
  use mod_param, only:  dl, dims
  ! use mod_bound, only: updthalo
  real(rp) :: vof_th = 0.5
  contains
  ! Import Threshold
  subroutine dropletTagging(n,vof,u,v,w, &
#if defined(_HEAT_TRANSFER)
                            tmp,         &
#endif
                            istep)
    implicit none
    integer,    intent(in),    dimension(3)                              :: n
    real(rp),   intent(in),    dimension(0:,0:,0:)                       :: vof,u,v,w
#if defined(_HEAT_TRANSFER)
    real(rp),   intent(in),    dimension(0:,0:,0:)                       :: tmp
#endif
    integer                                                              :: istep
    integer(8),                dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1)     :: vofTag
    ! integer(8),    intent(out),            dimension(0:,0:,0:)     :: vofTag0
    real(rp),                dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1)       :: vof_f
    real(rp)                                                             :: rdnInit
    integer                                                              :: i,j,k,b, nd, cellcount
    integer(4)                                                           :: id, idst
    integer(4), allocatable,    dimension(:,:)                           :: procShare, share,idShare, orderField
    integer(4), allocatable,    dimension(:)                             :: myid_out, drid
    ! integer(8), allocatable,    dimension(:,:)                           :: idShare
    integer(8), allocatable,    dimension(:)                             :: idloc
    ! integer,                    dimension(6)                             :: share
    real(8),    allocatable,    dimension(:)                             :: xd,yd,zd, ud, vd, wd, vold
#if defined(_HEAT_TRANSFER)
    real(8),    allocatable,    dimension(:)                             :: td
#endif
    character(len=3) ::proidch
    
    call random_number(rdnInit)
    ! idst = floor(rdnInit*100000+100000*myid) ! avoid possible ID duplication between processors. May fail for n_d >>10000???
    idst = 0 ! avoid possible ID duplication between processors. May fail for n_d >>10000???
    id = idst+1
    vof_f = vof
    vofTag = 0
    nd = 0
    ! main loop for tagging droplets 
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if (vof_f(i,j,k).ge.vof_th) then
            vof_f(i,j,k)  = 0.0
            vofTag(i,j,k) = id
            call recursiveTagging(vof_f,vofTag,id, i,j,k, n)
            id = id+1
            nd = nd+1
          endif
        enddo
      enddo
    enddo

    ! Update halo
    call updthalo((/n(1),n(2)/),1,vofTag)
    call updthalo((/n(1),n(2)/),2,vofTag)

    ! write(proidch,'(i3.3)') myid
    ! open(90+myid,file='proc'//proidch)
    ! do i=0,n(1)+1
    !   do j=0,n(2)+1
    !     do k=0,n(3)+1
    !       write(90+myid,'(4I5.4)') i,j,k,vofTag(i,j,k)
    !     enddo
    !   enddo
    ! enddo
    ! close(92)

    ! print*,'max vof', maxval(vof_f), myid

    ! allocate droplet properties
    allocate(xd(nd),yd(nd),zd(nd),ud(nd),vd(nd),wd(nd),vold(nd), & 
#if defined(_HEAT_TRANSFER)
             td(nd), &
#endif
             idloc(nd),share(nd,6),myid_out(nd),procShare(nd,6), &
             idShare(nd,6), orderField(nd,6),drid(nd))
    xd = 0.
    yd = 0.
    zd = 0.
    ud = 0.
    vd = 0.
    wd = 0.
    vold = 0.
#if defined(_HEAT_TRANSFER)
    td =0.
#endif
    id = 0
    cellcount = 0
    idShare = 0
    procShare = -1
    share = 0

    ! Loop for properties characterization
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if (vofTag(i,j,k).gt.0) then
            id       = vofTag(i,j,k)-idst
            drid(id) = id
            idloc(id)= vofTag(i,j,k) 
            xd(id)   = xd(id) + vof(i,j,k)*(i+coord(1)*n(1)-0.5)*dl(1)
            yd(id)   = yd(id) + vof(i,j,k)*(j+coord(2)*n(2)-0.5)*dl(2)
            zd(id)   = zd(id) + vof(i,j,k)*(k              -0.5)*dl(3)
            ud(id)   = ud(id) + u(i,j,k)*vof(i,j,k)
            vd(id)   = vd(id) + v(i,j,k)*vof(i,j,k)
            wd(id)   = wd(id) + w(i,j,k)*vof(i,j,k)
#if defined(_HEAT_TRANSFER)
            td(id)   = td(id) + tmp(i,j,k)*vof(i,j,k)
#endif
            vold(id) = vold(id) + vof(i,j,k)
            ! print*, vof(i,j,k), i,j,k
            cellcount = 1
            if ((i.eq.1).and. share(id,1).eq.0 .and. vofTag(0,j,k).ne.0) then
              ! print*, 'i0i0i0i0i'
              share(id,1)     = 1
              procShare(id,1) = left
              idShare(id,1)   = vofTag(0,j,k)
              ! if (vofTag(0,j,k).ne.0)      idShare(id,1)   = vofTag(0,j,k)
            else if ((i.eq.n(1)).and. share(id,2).eq.0 .and.vofTag(n(1)+1,j,k).ne.0) then
              ! print*, 'n1n1n1n1n1'
              share(id,2)     = 1
              procShare(id,2) = right
              idShare(id,2)   = vofTag(n(1)+1,j,k)
            else if ((j.eq.1) .and. share(id,3).eq.0 .and. vofTag(i,0,k).ne.0) then
              share(id,3)     = 1
              procShare(id,3) = front
              idShare(id,3)   = vofTag(i,0,k)
              ! print*, 'j1', idShare(id,:)
            else if ((j.eq.n(2)) .and. share(id,4).eq.0 .and. vofTag(i,n(2)+1,k).ne.0) then
              share(id,4)        = 1
              procShare(id,4) = back
              idShare(id,4)   = vofTag(i,n(2)+1,k)
              ! print*, 'n2', idShare(id,:)
            else if ((k.eq.1) .and. share(id,5).eq.0) then
              share(id,5)        = 1
              procShare(id,5) = myid
              idShare(id,5)   = vofTag(i,j,0)
            else if ((k.eq.n(3)).and. share(id,6).eq.0) then
              share(id,6)        = 1
              procShare(id,6) = myid
              idShare(id,6)   = vofTag(i,j,n(3)+1)
            endif
          endif
        enddo
      enddo
    enddo 

    do i=1,nd
      do j=1,6
        orderField(i,j) = myid*100+i*10+J
        ! print*,'proc', procShare(i,:), myid
        ! print*, nd, myid
        ! print*,'id  ', idShare(i,:), myid
      enddo

      ! print*,'id  ', idShare(i,:), myid
      ! print*,'proc', procShare(i,:), myid
    enddo

    myid_out = 0
    do id=1,nd
      myid_out(id) = myid
      xd(id)   = xd(id)/vold(id)    
      yd(id)   = yd(id)/vold(id)    
      zd(id)   = zd(id)/vold(id)    
      ud(id)   = ud(id)/vold(id)    
      vd(id)   = vd(id)/vold(id)    
      wd(id)   = wd(id)/vold(id)   
#if defined(_HEAT_TRANSFER)
      td(id)   = td(id)/vold(id)
#endif
      ! print*, xd(id),yd(id),zd(id),vold(id),id,myid 
      vold(id) = vold(id)*dl(1)*dl(2)*dl(3)
    enddo 

    ! Writing results
    call write_output_r8(istep,'xpos',xd  , nd)
    call write_output_r8(istep,'ypos',yd  , nd)
    call write_output_r8(istep,'zpos',zd  , nd)
    call write_output_r8(istep,'uvel',ud  , nd)
    call write_output_r8(istep,'vvel',vd  , nd)
    call write_output_r8(istep,'wvel',wd  , nd)
#if defined(_HEAT_TRANSFER)
    call write_output_r8(istep,'tmpd',td  , nd)
#endif
    call write_output_r8(istep,'vold',vold, nd)
    call write_output_i1(istep,'drid',drid, nd)
    call write_output_i1(istep,'idmy',myid_out, nd)
    call write_output_i6(istep,'proc',transpose(procShare), nd)
    call write_output_i6(istep,'orde',transpose(orderField), nd)
    call write_output_i6(istep,'idsh',transpose(idShare), nd)





    deallocate(xd,yd,zd,ud,vd,wd,vold, &
     procShare, idShare)
    return
  end subroutine dropletTagging

recursive subroutine recursiveTagging(vof_f,vofTag, id,i,j,k, n)
  implicit none
  real(rp),   intent(inout),     dimension(0:,0:,0:)         :: vof_f
  integer(8) ,   intent(inout),     dimension(0:,0:,0:)         :: vofTag
  integer (4),   intent(inout)                                  :: id
  integer ,   intent(in)                                  :: i,j,k
  integer ,   intent(in),        dimension(3)                :: n
  integer                                                    :: ii,jj,kk

  do kk=k-1,k+1
    do jj=j-1,j+1
      do ii=i-1,i+1
        if ((vof_f(ii,jj,kk).gt.vof_th).and.               &
          (ii.ge.1).and.(jj.ge.1).and.(kk.ge.1) .and.  &
          (ii.le.n(1)).and.(jj.le.n(2)).and.(kk.le.n(3))) then
          vof_f(ii,jj,kk) = 0.0
          vofTag(ii,jj,kk) = id
          ! print*, ii,jj,kk, vof_f(ii,jj,kk), vofTag(ii,jj,kk)
          call recursiveTagging(vof_f,vofTag,id, ii,jj,kk, n)
        endif
      enddo
    enddo
  enddo
  return
end subroutine recursiveTagging



! recursive subroutine recursiveProperties(i,j,k,id,n,vofTag,vof,u,v,w,xd,yd,zd,ud,vd,wd,vold,procShare,idShare,share,cellcount)
  ! integer,     intent(in)                           :: i,j,k
  ! integer,     intent(inout)                           :: cellcount
  ! integer,     intent(in),     dimension(3)         :: n
  ! integer(8),  intent(in)                           :: id
  ! integer(8),  intent(inout),  dimension(0:,0:,0:)  :: vofTag
  ! real(rp),    intent(in),     dimension(0:,0:,0:)  :: u,v,w,vof
  ! real(rp),       intent(inout),  dimension(1:)        :: xd,yd,zd,ud,vd,wd,vold
  ! integer(4),    intent(inout),  dimension(1:,1:)   :: procShare
  ! integer(8),    intent(inout),  dimension(1:,1:)   :: idShare
  ! integer(4),    intent(inout),  dimension(6)         :: share
  ! integer                                           :: ii,jj,kk
  ! do kk=k-1,k+1
  !   do jj=j-1,j+1
  !     do ii=i-1,i+1
  !       ! print*, ii,jj,kk
  !       if ((vofTag(ii,jj,kk).eq.id) .and. &
  !         (ii.ne.0)    .and.(jj.ne.0)    .and.(kk.ne.0) .and.  &
  !         (ii.ne.n(1)).and.(jj.ne.n(2)).and.(kk.ne.n(3))) then
  !         vofTag(ii,jj,kk) = 0
  !         xd(id)   = xd(id)   +vof(ii,jj,kk)*(ii+coord(1)*n(1)-0.5)*dl(1)
  !         yd(id)   = yd(id)   +vof(ii,jj,kk)*(jj+coord(2)*n(2)-0.5)*dl(2)
  !         zd(id)   = zd(id)   +vof(ii,jj,kk)*(kk              -0.5)*dl(3)
  !         ud(id)   = ud(id)   +u(ii,jj,kk)  *vof(ii,jj,kk)
  !         vd(id)   = vd(id)   +v(ii,jj,kk)  *vof(ii,jj,kk)
  !         wd(id)   = wd(id)   +w(ii,jj,kk)  *vof(ii,jj,kk)
  !         vold(id) = vold(id) +vof(ii,jj,kk)
  !         cellcount = cellcount +1
  !         if ((ii.eq.1) .and. (share(1).eq.0)) then
  !           share(1)        = 1
  !           procShare(id,1) = left
  !           idShare(id,1)   = vofTag(ii-1,jj,kk)
  !         elseif ((ii.eq.n(1)) .and. (share(2).eq.0)) then
  !           share(2)        = 1
  !           procShare(id,2) = right
  !           idShare(id,2)   = vofTag(ii+1,jj,kk)
  !         elseif ((jj.eq.1) .and. (share(3).eq.0)) then
  !           share(3)        = 1
  !           procShare(id,3) = back
  !           idShare(id,3)   = vofTag(ii,jj-1,kk)
  !         elseif ((jj.eq.n(2)) .and. (share(4).eq.0)) then
  !           share(4)        = 1
  !           procShare(id,4) = front
  !           idShare(id,4)   = vofTag(ii,jj+1,kk)
  !         elseif ((kk.eq.1) .and. (share(5).eq.0)) then
  !           share(5)        = 1
  !           procShare(id,5) = myid
  !           idShare(id,5)   = vofTag(ii,jj,kk-1)
  !         elseif ((jj.eq.n(2)) .and. (share(4).eq.0)) then
  !           share(6)        = 1
  !           procShare(id,6) = myid
  !           idShare(id,6)   = vofTag(ii,jj,kk)+1
  !         endif
  !         call recursiveProperties(i,j,k,id,n,vofTag,vof,u,v,w,xd,yd,zd,ud,vd,wd,vold,procShare,idShare,share,cellcount)


  !       endif
  !       ! endif
  !     enddo
  !   enddo
  ! enddo
  ! return

! end subroutine recursiveProperties


subroutine write_output_r8(istep,varname,var, nd)
  implicit none
  integer,            intent(in)                    :: istep, nd
  character(len=4),   intent(in)                    :: varname
  real(rp),           intent(in), dimension(1:)     :: var
  integer                                           :: offset, fh
  character(7)                                      :: fldnum
  integer(kind=MPI_OFFSET_KIND)                     :: filesize,disp,nreals,nreals_myid
  character(len=1024)                               :: fname

  write(fldnum,'(i7.7)') istep
  disp = 0
  fname = 'data/post/tagging/'//varname//'fld'//fldnum//'.bin'
  call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  disp = disp-nd
  call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
  filesize = 0_MPI_OFFSET_KIND
  call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
  call MPI_FILE_SET_VIEW(fh,disp*8, MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION, 'native', &
       MPI_INFO_NULL, ierr)
  call MPI_FILE_WRITE(fh,var,nd,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  call MPI_FILE_CLOSE(fh,ierr)

end subroutine write_output_r8


! subroutine write_output_i6(istep,varname,var, nd)
!   implicit none
!   integer,            intent(in)                    :: istep, nd
!   character(len=4),   intent(in)                    :: varname
!   integer(4),         intent(in), dimension(1:,1:)  :: var
!   integer                                           :: offset, fh, ndtot,nelem,disp, subar, ierr
!   character(7)                                      :: fldnum
!   integer(kind=MPI_OFFSET_KIND)                     :: filesize, displacement
!   character(len=1024)                               :: fname
!   integer, dimension(2) :: sizes,subsizes,starts

!   nelem = 6
!   write(fldnum,'(i7.7)') istep
!   disp = 0
!   fname = 'data/post/tagging/'//varname//'fld'//fldnum//'.bin'
!   ! if (nd.gt.0) then
!   call MPI_allreduce(nd,ndtot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!   call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!   disp = disp - nd
!   print*, ndtot, disp
!   sizes = (/ndtot,nelem/)
!   subsizes = (/nd,nelem/)
!   starts = (/disp,0/)
!   Call MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_INTEGER, subar ,ierr)
!   ! Call MPI_Type_create_subarray(2,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_real,subarray,ierr)
!   Call MPI_Type_commit(subar,ierr)
!   ! endif
!   call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
!   filesize = 0_MPI_OFFSET_KIND
!   ! displacement= disp*nelem
!   displacement= 0
!   call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
!   call MPI_FILE_SET_VIEW(fh,displacement, MPI_INTEGER, subar, 'native', &
!                           MPI_INFO_NULL, ierr)
!   ! displacement = disp-np
!   ! call MPI_FILE_SET_VIEW(fh,displacement, MPI_INT, subar, 'native', &
!   !      MPI_INFO_NULL, ierr)
!   ! print*, 'sub', subarray
!   call MPI_FILE_WRITE_all(fh,var,nelem*nd,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
!   call MPI_FILE_CLOSE(fh,ierr)

! end subroutine write_output_i6


subroutine write_output_i6(istep,varname,var, nd)
  implicit none
  integer,            intent(in)                    :: istep, nd
  character(len=4),   intent(in)                    :: varname
  integer(4),           intent(in), dimension(1:,1:)  :: var
  integer                                           :: offset, fh, ndtot,nelem
  integer, save   ::subarray
  character(7)                                      :: fldnum
  integer(kind=MPI_OFFSET_KIND)                     :: filesize,nreals,nreals_myid, displacement,disp
  character(len=1024)                               :: fname

  nelem = 6
  write(fldnum,'(i7.7)') istep
  disp = 0
  fname = 'data/post/tagging/'//varname//'fld'//fldnum//'.bin'
  call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  disp = disp-nd
  call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
  filesize = 0_MPI_OFFSET_KIND
  call MPI_FILE_SET_SIZE(fh,filesize,ierr)    
  call MPI_FILE_SET_VIEW(fh,disp*nelem*4, MPI_INTEGER, MPI_INTEGER, 'native', &
       MPI_INFO_NULL, ierr)
  call MPI_FILE_WRITE_all(fh,var,nelem*nd,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
  call MPI_FILE_CLOSE(fh,ierr)

end subroutine write_output_i6

subroutine write_output_i1(istep,varname, var, nd)
  implicit none
  integer,            intent(in)                    :: istep, nd
  character(len=4),   intent(in)                    :: varname
  integer(4),           intent(in), dimension(1:)  :: var
  integer                                           :: offset, fh, nelem
  character(7)                                      :: fldnum
  integer(kind=MPI_OFFSET_KIND)                     :: filesize,disp,nreals,nreals_myid
  character(len=1024)                               :: fname

  nelem = 1
  write(fldnum,'(i7.7)') istep
  disp = 0
  fname = 'data/post/tagging/'//varname//'fld'//fldnum//'.bin'
  call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  disp = disp-nd
  call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
  filesize = 0_MPI_OFFSET_KIND
  call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
  call MPI_FILE_SET_VIEW(fh,disp*4*nelem, MPI_INTEGER,MPI_INTEGER, 'native', &
       MPI_INFO_NULL, ierr)
  call MPI_FILE_WRITE(fh,var,nd,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
  call MPI_FILE_CLOSE(fh,ierr)

end subroutine write_output_i1

subroutine updthalo(n,idir,p)
  implicit none
  integer , dimension(2), intent(in) :: n
  integer , intent(in) :: idir
  integer(rp), dimension(0:,0:,0:), intent(inout) :: p
#if defined(_OPENACC)
  attributes(managed) :: p
  integer :: istat
#endif
  integer :: i,j,k,n_1,n_2
  !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
  !
  !  this subroutine updates the halos that store info
  !  from the neighboring computational sub-domain
  !
  select case(idir)
  case(1) ! x direction
    !if( .false. ) then
    if(dims(1) .eq.  1) then
#if defined(_OPENACC)
      n_1=n(1)
      !$cuf kernel do(2) <<<*,*>>>
      do k=lbound(p,3),ubound(p,3)
        do j=lbound(p,2),ubound(p,2)
          p(n_1+1 ,j,k) = p(  1,j,k)
          p(0     ,j,k) = p(n_1,j,k)
        enddo
      enddo
#else
      !$OMP WORKSHARE
      p(n(1)+1,:,:) = p(   1,:,:)
      p(0     ,:,:) = p(n(1),:,:)
      !$OMP END WORKSHARE
#endif
    else
      n_1=n(1)
#if defined(_OPENACC)
      !$cuf kernel do(2) <<<*,*>>>
#endif
      do k=lbound(p,3),ubound(p,3)
        do j=lbound(p,2),ubound(p,2)
          xsl_buf(j,k) = p(  1,j,k)
          xsr_buf(j,k) = p(n_1,j,k)
        enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize()
      !
      call MPI_SENDRECV(xsl_buf(0,0), size( xsl_buf ),MPI_REAL_RP,left ,0, &
        xrr_buf(0,0), size( xrr_buf ),MPI_REAL_RP,right,0, &
        comm_cart,status,ierr)
      call MPI_SENDRECV(xsr_buf(0,0), size( xsr_buf ),MPI_REAL_RP,right,0, &
        xrl_buf(0,0), size( xrl_buf ),MPI_REAL_RP,left ,0, &
        comm_cart,status,ierr)
#if defined(_OPENACC)
      !$cuf kernel do(2) <<<*,*>>>
#endif
      do k=lbound(p,3),ubound(p,3)
        do j=lbound(p,2),ubound(p,2)
          p(n_1+1,j,k) = xrr_buf(j,k)
          p(0    ,j,k) = xrl_buf(j,k)
        enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize()
      !
      !call MPI_SENDRECV(p(1     ,0,0),1,xhalo,left ,0, &
      !                  p(n(1)+1,0,0),1,xhalo,right,0, &
      !                  comm_cart,status,ierr)
      !call MPI_SENDRECV(p(n(1),0,0),1,xhalo,right,0, &
      !                  p(0   ,0,0),1,xhalo,left ,0, &
      !                  comm_cart,status,ierr)
      !!call MPI_IRECV(p(0     ,0,0),1,xhalo,left ,1, &
      !!               comm_cart,requests(2),error)
      !!call MPI_IRECV(p(n(1)+1,0,0),1,xhalo,right,0, &
      !!               comm_cart,requests(1),error)
      !!call MPI_ISSEND(p(n(1),0,0),1,xhalo,right,1, &
      !!               comm_cart,requests(4),error)
      !!call MPI_ISSEND(p(1   ,0,0),1,xhalo,left ,0, &
      !!               comm_cart,requests(3),error)
      !!call MPI_WAITALL(4, requests, statuses, error)
    endif
  case(2) ! y direction
    !if( .false. ) then
    if( dims(2) .eq.  1 ) then
#if defined(_OPENACC)
      n_2=n(2)
      !$cuf kernel do(2) <<<*,*>>>
      do k=lbound(p,3),ubound(p,3)
        do i=lbound(p,1),ubound(p,1)
          p(i,n_2+1,k) = p(i,  1,k)
          p(i,    0,k) = p(i,n_2,k)
        enddo
      enddo
#else
      !$OMP WORKSHARE
      p(:,n(2)+1,:)  = p(:, 1  ,:)
      p(:, 0     ,:) = p(:,n(2),:)
      !$OMP END WORKSHARE
#endif
    else
      n_2=n(2)
#if defined(_OPENACC)
      !$cuf kernel do(2) <<<*,*>>>
#endif
      do k=lbound(p,3),ubound(p,3)
        do i=lbound(p,1),ubound(p,1)
          ysl_buf(i,k) = p(i,  1,k)
          ysr_buf(i,k) = p(i,n_2,k)
        enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize()
      !
      call MPI_SENDRECV(ysl_buf(0,0), size( ysl_buf ),MPI_REAL_RP,front,0, &
        yrr_buf(0,0), size( yrr_buf ),MPI_REAL_RP,back ,0, &
        comm_cart,status,ierr)
      call MPI_SENDRECV(ysr_buf(0,0), size( ysr_buf ),MPI_REAL_RP,back ,0, &
        yrl_buf(0,0), size( yrl_buf ),MPI_REAL_RP,front,0, &
        comm_cart,status,ierr)
#if defined(_OPENACC)
      !$cuf kernel do(2) <<<*,*>>>
#endif
      do k=lbound(p,3),ubound(p,3)
        do i=lbound(p,1),ubound(p,1)
          p(i,n_2+1,k) = yrr_buf(i,k)
          p(i,    0,k) = yrl_buf(i,k)
        enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize()
      !call MPI_SENDRECV(p(0,1     ,0),1,yhalo,front,0, &
      !                  p(0,n(2)+1,0),1,yhalo,back ,0, &
      !                  comm_cart,status,ierr)
      !call MPI_SENDRECV(p(0,n(2),0),1,yhalo,back ,0, &
      !                  p(0,0   ,0),1,yhalo,front,0, &
      !                  comm_cart,status,ierr)
      !!call MPI_IRECV(p(0,n(2)+1,0),1,yhalo,back ,0, &
      !!               comm_cart,requests(1),error)
      !!call MPI_IRECV(p(0,0     ,0),1,yhalo,front,1, &
      !!               comm_cart,requests(2),error)
      !!call MPI_ISSEND(p(0,1   ,0),1,yhalo,front,0, &
      !!               comm_cart,requests(3),error)
      !!call MPI_ISSEND(p(0,n(2),0),1,yhalo,back ,1, &
      !!               comm_cart,requests(4),error)
      !!call MPI_WAITALL(4, requests, statuses, error)
    endif
  end select

  return
end subroutine updthalo

end module
