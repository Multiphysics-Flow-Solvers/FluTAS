module mod_tagging
  !
  ! note: not yet on GPU
  ! 
  use mpi
  use mod_bound     , only: boundp
  use mod_common_mpi, only: myid,comm_cart,ierr,ijk_start, &
                            left,right,front,back,top,bottom
  use mod_param     , only: cbcvof,bcvof
  use mod_types     , only: rp
  !
  implicit none 
  !
  real(rp), parameter :: vof_th = 0.5_rp
  !
  private
  public  :: droplet_tagging
  !
  contains
  !
  subroutine droplet_tagging(n,dims,datadir_ta,dl,nh_d,nh_v,nh_u,halo_v,dzc,dzf,vof,u,v,w,istep,time)
    !
    implicit none
    !
    integer         , intent(in), dimension(3)                       :: n
    integer         , intent(in), dimension(3)                       :: dims
    character(len=*), intent(in)                                     :: datadir_ta
    real(rp)        , intent(in), dimension(3)                       :: dl
    integer         , intent(in)                                     :: nh_d,nh_v,nh_u
    integer         , intent(in), dimension(3)                       :: halo_v
    real(rp)        , intent(in), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp)        , intent(in), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: vof
    real(rp)        , intent(in), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    integer         , intent(in)                                     :: istep
    real(rp)        , intent(in)                                     :: time
    !
    integer(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: vofTag
    real(rp)  , dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: vof_f,vofTag_r
    real(rp) :: rdnInit
    integer  :: i,j,k,b, nd, cellcount
    integer(4) :: id, idst
    integer(4), allocatable, dimension(:,:) :: procShare, share,idShare, orderField
    integer(4), allocatable, dimension(:)   :: myid_out, drid
    integer(8), allocatable, dimension(:)   :: idloc
    real(8),    allocatable, dimension(:)   :: xd,yd,zd,ud,vd,wd,vold
    character(len=3) :: proidch
    !@cuf attributes(managed) :: dzc, dzf, vofTag_r
    ! 
    !call random_number(rdnInit)
    ! idst = floor(rdnInit*100000+100000*myid) ! avoid possible ID duplication between processors. May fail for n_d >>10000???
    idst = 0 ! avoid possible ID duplication between processors. May fail for n_d >>10000???
    id   = idst+1
    nd   = 0
    !
    ! initialize vofTag and vof_f
    !
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          vofTag(i,j,k) = 0
          vof_f(i,j,k)  = vof(i,j,k)
        enddo
      enddo
    enddo
    !
    ! main loop for tagging droplets 
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(vof_f(i,j,k).ge.vof_th) then
            vof_f(i,j,k)  = 0._rp
            vofTag(i,j,k) = id
            call recursive_tagging(n,i,j,k,id,vof_f,vofTag)
            id = id+1
            nd = nd+1
          endif
        enddo
      enddo
    enddo
    !
    ! conversion to real(rp) to set_bc
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          vofTag_r(i,j,k) = 1._rp*vofTag(i,j,k)
        enddo
      enddo
    enddo
    call boundp(cbcvof,n,bcvof,nh_d,nh_v,halo_v,dl,dzc,dzf,vofTag_r) ! meaningful bc: N or P
    !
    ! re-conversion to integer
    !
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          vofTag(i,j,k) = int(vofTag_r(i,j,k),8)
        enddo
      enddo
    enddo
    ! 
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
    !
    ! print*,'max vof', maxval(vof_f), myid
    !
    ! allocate droplet properties
    !
    allocate(xd(nd),yd(nd),zd(nd),ud(nd),vd(nd),wd(nd),vold(nd), & 
             idloc(nd),share(nd,6),myid_out(nd),procShare(nd,6), &
             idShare(nd,6), orderField(nd,6),drid(nd))
    xd = 0._rp
    yd = 0._rp
    zd = 0._rp
    ud = 0._rp
    vd = 0._rp
    wd = 0._rp
    vold = 0._rp
    id = 0
    cellcount = 0
    idShare = 0
    procShare = -1
    share = 0
    !
    ! Loop for properties characterization
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          if(vofTag(i,j,k).gt.0) then
            !
            id        = vofTag(i,j,k)-idst
            drid(id)  = id
            idloc(id) = vofTag(i,j,k) 
            xd(id)    = xd(id) + vof(i,j,k)*(i+ijk_start(1)-0.5_rp)*dl(1)
            yd(id)    = yd(id) + vof(i,j,k)*(j+ijk_start(2)-0.5_rp)*dl(2)
            zd(id)    = zd(id) + vof(i,j,k)*(k+ijk_start(3)-0.5_rp)*dl(3)
            ud(id)    = ud(id) + u(i,j,k)*vof(i,j,k)
            vd(id)    = vd(id) + v(i,j,k)*vof(i,j,k)
            wd(id)    = wd(id) + w(i,j,k)*vof(i,j,k)
            vold(id)  = vold(id) + vof(i,j,k)
            !
            ! print*, vof(i,j,k), i,j,k
            cellcount = 1
            if((i.eq.1).and.share(id,1).eq.0.and.vofTag(0,j,k).ne.0) then
              ! print*, 'i0i0i0i0i'
              share(id,1)     = 1
              procShare(id,1) = left
              idShare(id,1)   = vofTag(0,j,k)
              ! if (vofTag(0,j,k).ne.0)      idShare(id,1)   = vofTag(0,j,k)
            elseif((i.eq.n(1)).and.share(id,2).eq.0.and.vofTag(n(1)+1,j,k).ne.0) then
              ! print*, 'n1n1n1n1n1'
              share(id,2)     = 1
              procShare(id,2) = right
              idShare(id,2)   = vofTag(n(1)+1,j,k)
            elseif((j.eq.1).and.share(id,3).eq.0.and.vofTag(i,0,k).ne.0) then
              share(id,3)     = 1
              procShare(id,3) = front
              idShare(id,3)   = vofTag(i,0,k)
              ! print*, 'j1', idShare(id,:)
            elseif((j.eq.n(2)).and.share(id,4).eq.0.and.vofTag(i,n(2)+1,k).ne.0) then
              share(id,4)     = 1
              procShare(id,4) = back
              idShare(id,4)   = vofTag(i,n(2)+1,k)
              ! print*, 'n2', idShare(id,:)
            elseif((k.eq.1).and.share(id,5).eq.0.and.vofTag(i,j,0).ne.0) then
              share(id,5)     = 1
              procShare(id,5) = bottom
              idShare(id,5)   = vofTag(i,j,0)
            elseif((k.eq.n(3)).and.share(id,6).eq.0.and.vofTag(i,j,n(3)+1).ne.0) then
              share(id,6)     = 1
              procShare(id,6) = top
              idShare(id,6)   = vofTag(i,j,n(3)+1)
            endif
          endif
        enddo
      enddo
    enddo 
    !
    do i=1,nd
      do j=1,6
        orderField(i,j) = myid*100+i*10+J
        ! print*,'proc', procShare(i,:), myid
        ! print*, nd, myid
        ! print*,'id  ', idShare(i,:), myid
      enddo
    enddo
    !
    myid_out = 0
    do id=1,nd
      myid_out(id) = myid
      xd(id)   = xd(id)/vold(id)    
      yd(id)   = yd(id)/vold(id)    
      zd(id)   = zd(id)/vold(id)    
      ud(id)   = ud(id)/vold(id)    
      vd(id)   = vd(id)/vold(id)    
      wd(id)   = wd(id)/vold(id)   
      vold(id) = vold(id)*dl(1)*dl(2)*dl(3)
    enddo 
    !
    ! Writing results
    !
    call write_output_r8(istep,datadir_ta,'xpos',xd  , nd)
    call write_output_r8(istep,datadir_ta,'ypos',yd  , nd)
    call write_output_r8(istep,datadir_ta,'zpos',zd  , nd)
    call write_output_r8(istep,datadir_ta,'uvel',ud  , nd)
    call write_output_r8(istep,datadir_ta,'vvel',vd  , nd)
    call write_output_r8(istep,datadir_ta,'wvel',wd  , nd)
    call write_output_r8(istep,datadir_ta,'vold',vold, nd)
    !
    call write_output_i1(istep,datadir_ta,'drid',drid, nd)
    call write_output_i1(istep,datadir_ta,'idmy',myid_out, nd)
    call write_output_i6(istep,datadir_ta,'proc',transpose(procShare), nd)
    call write_output_i6(istep,datadir_ta,'orde',transpose(orderField), nd)
    call write_output_i6(istep,datadir_ta,'idsh',transpose(idShare), nd)
    !
    ! deallocate the droplet variables
    !
    deallocate(xd,yd,zd,ud,vd,wd,vold,procShare,idShare)
    !
    return
  end subroutine droplet_tagging
  !
  recursive subroutine recursive_tagging(n,i,j,k,id,vof_f,vofTag)
    !
    implicit none
    !
    integer   , intent(in   ), dimension(3)        :: n
    integer   , intent(in   )                      :: i,j,k
    integer(4), intent(inout)                      :: id
    real(rp)  , intent(inout), dimension(0:,0:,0:) :: vof_f
    integer(8), intent(inout), dimension(0:,0:,0:) :: vofTag
    !
    integer :: ii,jj,kk
    ! 
    do kk=k-1,k+1
      do jj=j-1,j+1
        do ii=i-1,i+1
          !
          if((vof_f(ii,jj,kk).gt.vof_th).and.(ii.ge.1   ).and.(jj.ge.1   ).and.(kk.ge.1    ).and. &
                                             (ii.le.n(1)).and.(jj.le.n(2)).and.(kk.le.n(3))) then
            vof_f(ii,jj,kk)  = 0._rp
            vofTag(ii,jj,kk) = id
            ! print*, ii,jj,kk, vof_f(ii,jj,kk), vofTag(ii,jj,kk)
            call recursive_tagging(n,ii,jj,kk,id,vof_f,vofTag)
            !
          endif
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine recursive_tagging
  !
  subroutine write_output_r8(istep,datadir_ta,varname,var,nd)
    !
    implicit none
    !
    integer,          intent(in)                :: istep
    character(len=*), intent(in)                :: datadir_ta
    character(len=*), intent(in)                :: varname
    real(rp),         intent(in), dimension(1:) :: var
    integer,          intent(in)                :: nd
    !
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,nreals,nreals_myid
    integer :: offset,fh
    character(9) :: fldnum
    character(len=1024) :: fname
    !
    write(fldnum,'(i9.9)') istep
    disp = 0
    fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    !
    call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    disp = disp-nd
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
    call MPI_FILE_SET_VIEW(fh,disp*8, MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,var,nd,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine write_output_r8
  !
  subroutine write_output_i6(istep,datadir_ta,varname,var,nd)
    !
    implicit none
    !
    integer,          intent(in)                   :: istep, nd
    character(len=*), intent(in)                   :: datadir_ta
    character(len=4), intent(in)                   :: varname
    integer(4)      , intent(in), dimension(1:,1:) :: var
    !
    integer(kind=MPI_OFFSET_KIND) :: filesize,nreals,nreals_myid, displacement,disp
    integer :: offset, fh, ndtot,nelem
    integer, save :: subarray
    character(9) :: fldnum
    character(len=1024) :: fname
    ! 
    nelem = 6
    write(fldnum,'(i9.9)') istep
    disp = 0
    fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    disp = disp-nd
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)    
    call MPI_FILE_SET_VIEW(fh,disp*nelem*4, MPI_INTEGER, MPI_INTEGER, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE_all(fh,var,nelem*nd,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine write_output_i6
  !
  subroutine write_output_i1(istep,datadir_ta,varname,var,nd)
    !
    implicit none
    !
    integer,          intent(in)                :: istep, nd
    character(len=*), intent(in)                :: datadir_ta
    character(len=4), intent(in)                :: varname
    integer(4),       intent(in), dimension(1:) :: var
    !
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,nreals,nreals_myid
    integer :: offset,fh,nelem
    character(9) :: fldnum
    character(len=1024) :: fname
    ! 
    nelem = 1
    write(fldnum,'(i9.9)') istep
    disp = 0
    fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    disp = disp-nd
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
    call MPI_FILE_SET_VIEW(fh,disp*4*nelem, MPI_INTEGER,MPI_INTEGER, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,var,nd,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine write_output_i1
  !
end module
