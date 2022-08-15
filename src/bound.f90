!
! SPDX-License-Identifier: MIT
!
module mod_bound
  !
  use mpi
  use mod_types
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: bounduvw, boundp, updt_rhs_b
  !
  contains
  !
  subroutine bounduvw(cbc,n,bc,nh_d,nh_u,halo,isoutflow,dl,dzc,dzf,u,v,w)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    ! imposes velocity boundary conditions 
    ! on an arbitrary stecil size
    !
    implicit none
    !
    character(len=1), intent(in   ), dimension(0:1,3,3)                 :: cbc
    integer         , intent(in   ), dimension(3)                       :: n 
    real(rp)        , intent(in   ), dimension(0:1,3,3)                 :: bc
    integer         , intent(in   )                                     :: nh_d,nh_u
    integer         , intent(in   ), dimension(3)                       :: halo
    logical         , intent(in   ), dimension(0:1,3)                   :: isoutflow
    real(rp)        , intent(in   ), dimension(3)                       :: dl
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp), dimension(0:nh_u-1) :: dr
    real(rp) :: dx,dy,dz
    integer  :: q,idir,sgn,ioutflowdir,qmin
    integer  :: ind1,ind2
    !@cuf attributes(managed) :: u, v, w, dzc, dzf
    !
    dx = dl(1)
    dy = dl(2)
    dz = dl(3)
    !
    qmin = abs(1-nh_u)
    !
#if defined(_DECOMP_X)
    ind1 = 2
    ind2 = 3
#elif _DECOMP_Y
    ind1 = 1
    ind2 = 3
#else
    ind1 = 1
    ind2 = 2
#endif
    !
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,u)
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,u)
    !
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,v)
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,v)
    !
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,w)
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,w)
    !
    ! along x
    !
    do q = 0,qmin
      dr(q) = dl(1)
    enddo
    if(left .eq.MPI_PROC_NULL) then ! x - bottom boundary
      call set_bc(n(1),n(2),n(3),cbc(0,1,1),0,1,.false.,bc(0,1,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(0,1,2),0,1,.true. ,bc(0,1,2),qmin,nh_u,dr,v)
      call set_bc(n(1),n(2),n(3),cbc(0,1,3),0,1,.true. ,bc(0,1,3),qmin,nh_u,dr,w)
    endif
    if(right.eq.MPI_PROC_NULL) then ! x - top boundary
      call set_bc(n(1),n(2),n(3),cbc(1,1,1),1,1,.false.,bc(1,1,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(1,1,2),1,1,.true. ,bc(1,1,2),qmin,nh_u,dr,v)
      call set_bc(n(1),n(2),n(3),cbc(1,1,3),1,1,.true. ,bc(1,1,3),qmin,nh_u,dr,w)
    endif
    !
    ! along y
    !
    do q = 0,qmin
      dr(q) = dl(2)
    enddo
    if(front.eq.MPI_PROC_NULL) then ! y - bottom boundary
      call set_bc(n(1),n(2),n(3),cbc(0,2,1),0,2,.true. ,bc(0,2,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(0,2,2),0,2,.false.,bc(0,2,2),qmin,nh_u,dr,v)
      call set_bc(n(1),n(2),n(3),cbc(0,2,3),0,2,.true. ,bc(0,2,3),qmin,nh_u,dr,w)
    endif
    if(back .eq.MPI_PROC_NULL) then ! y - top boundary
      call set_bc(n(1),n(2),n(3),cbc(1,2,1),1,2,.true. ,bc(1,2,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(1,2,2),1,2,.false.,bc(1,2,2),qmin,nh_u,dr,v)
      call set_bc(n(1),n(2),n(3),cbc(1,2,3),1,2,.true. ,bc(1,2,3),qmin,nh_u,dr,w)
    endif
    !
    ! along z
    !
    do q = 0,qmin
      dr(q) = dzc(-q)
    enddo
    if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
      call set_bc(n(1),n(2),n(3),cbc(0,3,1),0,3,.true. ,bc(0,3,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(0,3,2),0,3,.true. ,bc(0,3,2),qmin,nh_u,dr,v)
    endif
    do q = 0,qmin
      dr(q) = dzf(-q)
    enddo
    if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
      call set_bc(n(1),n(2),n(3),cbc(0,3,3),0,3,.false.,bc(0,3,3),qmin,nh_u,dr,w)
    endif
    !
    do q = 0,qmin 
      dr(q) = dzc(n(3)+q)
    enddo
    if(top.eq.MPI_PROC_NULL) then ! z - top boundary 
      call set_bc(n(1),n(2),n(3),cbc(1,3,1),1,3,.true. ,bc(1,3,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(1,3,2),1,3,.true. ,bc(1,3,2),qmin,nh_u,dr,v)
    endif
    do q = 0,qmin 
      dr(q) = dzf(n(3)+q)
    enddo
    if(top.eq.MPI_PROC_NULL) then ! z - top boundary
      call set_bc(n(1),n(2),n(3),cbc(1,3,3),1,3,.false.,bc(1,3,3),qmin,nh_u,dr,w)
    endif
    !
    ! outflow bc.
    !
    do q = 1,3
      do idir = 0,1
        if(isoutflow(idir,q)) then
          if(idir.eq.0) sgn = -1
          if(idir.eq.1) sgn = +1
          ioutflowdir = q*sgn
          call outflow(n(1),n(2),n(3),ioutflowdir,nh_d,nh_u,dx,dy,dz,dzf,u,v,w)
        endif
      enddo
    enddo
    !
    return
  end subroutine bounduvw
  !
  subroutine boundp(cbc,n,bc,nh_d,nh_p,halo,dl,dzc,dzf,p)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    ! imposes pressure and scalar boundary conditions
    ! on an arbitrary stencil size
    !
    implicit none
    !
    character(len=1), intent(in   ), dimension(0:1,3)                   :: cbc
    integer         , intent(in   ), dimension(3)                       :: n
    real(rp)        , intent(in   ), dimension(0:1,3)                   :: bc
    integer         , intent(in   )                                     :: nh_d,nh_p 
    integer         , intent(in   ), dimension(3)                       :: halo 
    real(rp)        , intent(in   ), dimension(3)                       :: dl
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp)        , intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    !
    real(rp), dimension(0:nh_p-1) :: dr
    integer :: q,qmin
    integer :: ind1,ind2
    !@cuf attributes(managed) :: p, dzc, dzf
    !
    qmin = abs(1-nh_p)
    !
#if defined(_DECOMP_X)
    ind1 = 2
    ind2 = 3
#elif _DECOMP_Y
    ind1 = 1
    ind2 = 3
#else
    ind1 = 1
    ind2 = 2
#endif
    !
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind1),ind1,p)
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind2),ind2,p)
    !
    ! along x
    !
    do q = 0,qmin
      dr(q) = dl(1)
    enddo
    if(left .eq.MPI_PROC_NULL) then ! x - bottom part
      call set_bc(n(1),n(2),n(3),cbc(0,1),0,1,.true.,bc(0,1),qmin,nh_p,dr,p)
    endif
    if(right.eq.MPI_PROC_NULL) then ! x - upper part
      call set_bc(n(1),n(2),n(3),cbc(1,1),1,1,.true.,bc(1,1),qmin,nh_p,dr,p)
    endif
    !
    ! along y
    !
    do q = 0,qmin
      dr(q) = dl(2)
    enddo
    if(front.eq.MPI_PROC_NULL) then ! y - bottom part
      call set_bc(n(1),n(2),n(3),cbc(0,2),0,2,.true.,bc(0,2),qmin,nh_p,dr,p)
     endif
    if(back .eq.MPI_PROC_NULL) then ! y - upper part
      call set_bc(n(1),n(2),n(3),cbc(1,2),1,2,.true.,bc(1,2),qmin,nh_p,dr,p)
    endif
    !
    ! along z
    !
    do q=0,qmin 
      dr(q) = dzc(-q)
    enddo
    if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
      call set_bc(n(1),n(2),n(3),cbc(0,3),0,3,.true.,bc(0,3),qmin,nh_p,dr,p)
    endif
    do q=0,qmin 
      dr(q) = dzc(n(3)+q)
    enddo
    if(top.eq.MPI_PROC_NULL) then ! z - top boundary
      call set_bc(n(1),n(2),n(3),cbc(1,3),1,3,.true.,bc(1,3),qmin,nh_p,dr,p)
    endif
    !
    return
  end subroutine boundp
  !
  subroutine set_bc(nx,ny,nz,ctype,ibound,idir,centered,rvalue,qq_d,nh_p,dr,p)
    !
    ! note: to be generalized with non-homogeneous Neumann
    !
    implicit none
    !
    integer         , intent(in   )                                     :: nx,ny,nz
    character(len=1), intent(in   )                                     :: ctype
    integer         , intent(in   )                                     :: ibound,idir
    logical         , intent(in   )                                     :: centered
    real(rp)        , intent(in   )                                     :: rvalue
    integer         , intent(in   )                                     :: qq_d,nh_p
    real(rp)        , intent(in   ), dimension(0:qq_d)                  :: dr
    real(rp)        , intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    !
    real(rp), dimension(0:qq_d) :: factor
    real(rp) :: factor_value, sgn
    integer  :: i,j,k,q
    !@cuf attributes(managed) :: p
    !
    do q=0,qq_d
      factor(q) = rvalue
    enddo
    sgn = 0._rp
    !
    if(ctype.eq.'D'.and.centered) then
      do q=0,qq_d
        factor(q) = +2._rp*factor(q)
      enddo
      sgn = -1._rp
    endif
    if(ctype.eq.'N') then
      if(    ibound.eq.0) then
        do q=0,qq_d
          factor(q) = -dr(q)*factor(q)
        enddo
      elseif(ibound.eq.1) then
        do q=0,qq_d
          factor(q) =  dr(q)*factor(q)
        enddo
      endif
      sgn = +1._rp
    endif
    !
    select case(ctype)
    case('P')
      select case(idir)
      case(1)
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k,q,factor,sgn) &
        !$OMP SHARED(ny,nz,nh_p,p)
        !$acc parallel loop collapse(3)
        do k=1-nh_p,nz+nh_p
          do j=1-nh_p,ny+nh_p
            do q=0,nh_p-1
              p(0-q   ,j,k) = p(nx-q,j,k)
              p(nx+1+q,j,k) = p(1+q ,j,k)
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$OMP END PARALLEL DO
      case(2)
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k,q,factor,sgn) &
        !$OMP SHARED(nx,nz,nh_p,p)
        !$acc parallel loop collapse(3)
        do k=1-nh_p,nz+nh_p
          do q=0,nh_p-1
            do i=1-nh_p,nx+nh_p
              p(i,0-q   ,k) = p(i,ny-q,k)
              p(i,ny+1+q,k) = p(i,1+q ,k)
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$OMP END PARALLEL DO
      case(3)
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,j,q,factor,sgn) &
        !$OMP SHARED(nx,ny,nz,nh_p,p)
        !$acc parallel loop collapse(3)
        do q=0,nh_p-1
          do j=1-nh_p,ny+nh_p
            do i=1-nh_p,nx+nh_p
              p(i,j,0-q   ) = p(i,j,nz-q)
              p(i,j,nz+1+q) = p(i,j, 1+q)
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$OMP END PARALLEL DO
      end select
    case('D','N')
      if(centered) then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(0-q,j,k) = factor_value+sgn*p(1+q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(nx+1+q,j,k) = factor_value+sgn*p(nx-q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,0-q  ,k) = factor_value+sgn*p(i,1+q,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,ny+1+q,k) = factor_value+sgn*p(i,ny-q,k)
                enddo
              enddo 
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,0-q ) = factor_value+sgn*p(i,j,1+q)
                enddo 
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,nz+1+q) = factor_value+sgn*p(i,j,nz-q)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'D') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(0-q,j,k) = factor_value
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(nx+q  ,j,k) = factor_value
                  p(nx+1+q,j,k) = p(nx-1-q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,0-q,k) = factor_value
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,ny+q  ,k) = factor_value
                  p(i,ny+1+q,k) = p(i,ny-1-q,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,0-q) = factor_value 
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,nz+q  ) = factor_value
                  p(i,j,nz+1+q) = p(i,j,nz-1-q)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'N') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(0-q,j,k) = 1.0_rp*factor_value + p(1+q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(nx+q  ,j,k) = 1.0_rp*factor_value + p(nx-1-q,j,k)
                  p(nx+1+q,j,k) = 2.0_rp*factor_value + p(nx-1-q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,0-q  ,k) = 1.0_rp*factor_value + p(i,1+q  ,k) 
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,ny+q  ,k) = 1.0_rp*factor_value + p(i,ny-1-q,k)
                  p(i,ny+1+q,k) = 2.0_rp*factor_value + p(i,ny-1-q,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,0-q  ) = 1.0_rp*factor_value + p(i,j,1+q  )
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,nz+q  ) = 1.0_rp*factor_value + p(i,j,nz-1-q)
                  p(i,j,nz+1+q) = 2.0_rp*factor_value + p(i,j,nz-1-q)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        end select
      endif
    end select
    !
    return
  end subroutine set_bc
  !
  ! TODO: Add GPU here? 
  subroutine outflow(nx,ny,nz,idir,nh_d,nh_u,dx,dy,dz,dzf,u,v,w)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    implicit none
    !
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   )                                     :: idir
    integer , intent(in   )                                     :: nh_d,nh_u
    real(rp), intent(in   )                                     :: dx,dy,dz
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzf
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp), dimension(1-nh_d:nz+nh_d) :: dzfi
    real(rp) :: dxi,dyi
    integer  :: i,j,k,q,qmin
    !
    qmin = abs(1-nh_u)
    !
    dxi  = dx**(-1)
    dyi  = dy**(-1)
    dzfi = dzf**(-1)
    !
    ! determine face velocity from zero divergence
    ! Note: here we assume that in outflow boundaries, 
    !       there is not phase change or any other volume sources.
    !
    select case(idir)
    case(1) ! x direction, right
      if(right.eq.MPI_PROC_NULL) then
        !i = nx + 1
        i = nx + 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k,q) &
        !$OMP SHARED(nz,ny,i,qmin,u,v,w,dx,dyi,dzfi)
        do k=1,nz
          do j=1,ny
            do q=0,qmin+1 ! note that the value at u(i+qmin+1,j,k) is not used
              u(i+q,j,k) = u(i-1-q,j,k) - dx*((v(i+q,j,k)-v(i+q,j-1,k))*dyi+(w(i+q,j,k)-w(i+q,j,k-1))*dzfi(k))
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
    case(2) ! y direction, back
      if(back.eq.MPI_PROC_NULL) then
        !j = ny + 1
        j = ny + 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k,q) &
        !$OMP SHARED(nx,nz,j,qmin,u,v,w,dy,dxi,dzfi)
        do k=1,nz
          do q=0,qmin+1 ! note that the value at v(i,j+qmin+1,k) is not used
            do i=1,nx
              v(i,j+q,k) = v(i,j-1-q,k) - dy*((u(i,j+q,k)-u(i-1,j+q,k))*dxi+(w(i,j+q,k)-w(i,j+q,k-1))*dzfi(k))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(3) ! z direction, top
      if(top.eq.MPI_PROC_NULL) then
        !k = n(3) + 1
        k = nz + 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,j,q) &
        !$OMP SHARED(nx,ny,k,qmin,u,v,w,dzf,dxi,dyi)
        do q=0,qmin+1 ! note that the value at w(i,j,k+qmin+1) is not used
          do j=1,ny
            do i=1,nx
              w(i,j,k+q) = w(i,j,k-1-q) - dzf(k+q)*((u(i,j,k+q)-u(i-1,j,k+q))*dxi+(v(i,j,k+q)-v(i,j-1,k+q))*dyi)
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-1) ! x direction, left
      if(left.eq.MPI_PROC_NULL) then
        i = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k,q) &
        !$OMP SHARED(ny,nz,i,qmin,u,v,w,dx,dyi,dzfi)
        do k=1,nz
          do j=1,ny
            do q=0,qmin
              u(i-q,j,k) = u(i+1+q,j,k) + dx*((v(i+1+q,j,k)-v(i+1+q,j-1,k))*dyi+(w(i+1+q,j,k)-w(i+1+q,j,k-1))*dzfi(k))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-2) ! y direction, front
      if(front.eq.MPI_PROC_NULL) then
        j = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k,q) &
        !$OMP SHARED(nx,nz,j,qmin,u,v,w,dy,dxi,dzfi)
        do k=1,nz
          do q=0,qmin
            do i=1,nx
              v(i,j-q,k) = v(i,j+1+q,k) + dy*((u(i,j+1+q,k)-u(i-1,j+1+q,k))*dxi+(w(i,j+1+q,k)-w(i,j+1+q,k-1))*dzfi(k))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-3) ! z direction, bottom
      if(bottom.eq.MPI_PROC_NULL) then
        k = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,j,q) &
        !$OMP SHARED(nx,ny,k,qmin,u,v,w,dzf,dxi,dyi)
        do q=0,qmin
          do j=1,ny
            do i=1,nx
              w(i,j,k-q) = w(i,j,k+1+q) + dzf(k-q)*((u(i,j,k+1+q)-u(i-1,j,k+1+q))*dxi+(v(i,j,k+1+q)-v(i,j-1,k+1+q))*dyi)
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    end select
    !
    return
  end subroutine outflow
  !
  ! TODO: add GPU here?
  subroutine inflow(nx,ny,nz,idir,nh_u,vel2d,u,v,w)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    implicit none
    !
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   )                                     :: idir
    integer , intent(in   )                                     :: nh_u
    real(rp), intent(in   ), dimension(     0:,    0:)          :: vel2d
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    integer :: i,j,k,q
    !
    select case(idir)
    case(1) ! x direction
      i = 0
      if(left.eq.MPI_PROC_NULL) then
        do k=1,nz
          do j=1,ny
            do q=0,nh_u-1
              u(i-q,j,k) = vel2d(j,k)
            enddo
          enddo
        enddo 
      endif
    case(2) ! y direction
      j = 0
      if(front.eq.MPI_PROC_NULL) then
        do k=1,nz
          do q=0,nh_u-1
            do i=1,nx
              v(i,j-q,k) = vel2d(i,k)
            enddo
          enddo
        enddo 
      endif
    case(3) ! z direction
      k = 0
      if(bottom.eq.MPI_PROC_NULL) then
        do q=0,nh_u-1
          do j=1,ny
            do i=1,nx
              w(i,j,k-q) = vel2d(i,j)
            enddo
          enddo
        enddo
      endif 
    end select
    !
    return
  end subroutine inflow
  !
  subroutine updt_rhs_b(nx,ny,nz,c_or_f,cbc,nh_p,rhsbx,rhsby,rhsbz,p)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    implicit none
    !
    integer         , intent(in   )                                     :: nx,ny,nz
    character       , intent(in   ), dimension(3)                       :: c_or_f
    character(len=1), intent(in   ), dimension(0:1,3)                   :: cbc
    integer         , intent(in   )                                     :: nh_p
    real(rp)        , intent(in   ), dimension(      :,      :,     0:) :: rhsbx,rhsby,rhsbz
    real(rp)        , intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    !
    integer :: q1, q2, q3
    integer :: i,j,k,idir
    !@cuf attributes(managed) :: p, rhsbx, rhsby, rhsbz
    !
    q1 = 0
    q2 = 0
    q3 = 0
    if(c_or_f(1).eq.'f'.and.cbc(1,1).eq.'D') q1 = 1
    if(c_or_f(2).eq.'f'.and.cbc(1,2).eq.'D') q2 = 1
    if(c_or_f(3).eq.'f'.and.cbc(1,3).eq.'D') q3 = 1
    !
    ! along x
    !
    if(left.eq.MPI_PROC_NULL) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(j,k) &
      !$OMP SHARED(ny,nz,p,rhsbx)
      !$acc kernels
      do k=1,nz
        do j=1,ny
          p(1,j,k) = p(1,j,k) + rhsbx(j,k,0)
        enddo
      enddo
      !$acc end kernels
      !$OMP END PARALLEL DO
    endif  
    !
    if(right.eq.MPI_PROC_NULL) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(j,k,q) &
      !$OMP SHARED(ny,nx,nz,p,rhsbx)
      !$acc kernels
      do k=1,nz
        do j=1,ny
          p(nx-q1,j,k) = p(nx-q1,j,k) + rhsbx(j,k,1)
        enddo
      enddo
      !$acc end kernels
      !$OMP END PARALLEL DO
    endif
    !
    ! along y
    !
    if(front.eq.MPI_PROC_NULL) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,k) &
      !$OMP SHARED(nx,nz,p,rhsby)
      !$acc kernels
      do k=1,nz
        do i=1,nx
          p(i,1,k) = p(i,1,k) + rhsby(i,k,0)
        enddo
      enddo
      !$acc end kernels
      !$OMP END PARALLEL DO
    endif
    !
    if(back.eq.MPI_PROC_NULL) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,k,q) &
      !$OMP SHARED(nx,ny,nz,p,rhsby)
      !$acc kernels
      do k=1,nz
        do i=1,nx
          p(i,ny-q2,k) = p(i,ny-q2,k) + rhsby(i,k,1)
        enddo
      enddo
      !$acc end kernels
      !$OMP END PARALLEL DO
    endif
    !
    ! along z
    !
    if(bottom.eq.MPI_PROC_NULL) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,k,q) &
      !$OMP SHARED(nx,ny,p,rhsbz)
      !$acc kernels
      do j=1,ny
        do i=1,nx
          p(i,j,1      ) = p(i,j,1      ) + rhsbz(i,j,0)
        enddo
      enddo
      !$acc end kernels
      !$OMP END PARALLEL DO
    endif
    !
    if(top.eq.MPI_PROC_NULL) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,k,q) &
      !$OMP SHARED(nx,ny,nz,p,rhsbz)
      !$acc kernels
      do j=1,ny
        do i=1,nx
          p(i,j,nz-q3) = p(i,j,nz-q3) + rhsbz(i,j,1)
        enddo
      enddo
      !$acc end kernels
      !$OMP END PARALLEL DO
    endif
    !
    return
  end subroutine updt_rhs_b
  !
  subroutine updthalo(nx,ny,nz,nh,halo,idir,p)
    !
#if defined(_OPENACC)
    use mod_common_mpi, only: xsl_buf, xrl_buf, xsr_buf, xrr_buf, &
                              ysr_buf, yrr_buf, ysl_buf, yrl_buf, &
                              zsr_buf, zrr_buf, zsl_buf, zrl_buf
#endif
    use mod_common_mpi, only: left,right,front,back,top,bottom, &
                              comm_cart,status,ierr
    !
    implicit none
    !
    integer , intent(in   )                               :: nx,ny,nz
    integer , intent(in   )                               :: nh,halo,idir
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    !@cuf attributes(managed) :: p
    !
    integer :: i,j,k, q, lb1, lb2
    !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  this subroutine updates the halos that store info
    !  from the neighboring computational sub-domain
    !
    select case(idir)
    case(1) ! x direction
#if defined(_OPENACC)
      !$acc kernels present(xsl_buf,xrl_buf,xsr_buf,xrr_buf)
      xsl_buf(:,:) = - huge(1._rp)
      xrl_buf(:,:) = - huge(1._rp)
      xsr_buf(:,:) = - huge(1._rp)
      xrr_buf(:,:) = - huge(1._rp)
      !$acc end kernels
      !
      do q=0,nh-1
        !$acc kernels present(xsl_buf,xsr_buf)
        do k=1-nh,nz+nh
          do j=1-nh,ny+nh
            xsl_buf(j,k) = p(1+q ,j,k)
            xsr_buf(j,k) = p(nx-q,j,k)
          enddo
        enddo
        !$acc end kernels
        !
        lb1 = lbound(xsl_buf,1) 
        lb2 = lbound(xsl_buf,2) 
        call MPI_SENDRECV(xsl_buf(lb1,lb2), size( xsl_buf ),MPI_REAL_RP,left ,0, &
                          xrr_buf(lb1,lb2), size( xrr_buf ),MPI_REAL_RP,right,0, &
                          comm_cart,status,ierr)
        call MPI_SENDRECV(xsr_buf(lb1,lb2), size( xsr_buf ),MPI_REAL_RP,right,0, &
                          xrl_buf(lb1,lb2), size( xrl_buf ),MPI_REAL_RP,left ,0, &
                          comm_cart,status,ierr)
        !
        !$acc kernels present(xrl_buf,xrr_buf)                  
        do k=1-nh,nz+nh
          do j=1-nh,ny+nh
            p(nx+1+q  ,j,k) = xrr_buf(j,k)
            p(0-q     ,j,k) = xrl_buf(j,k)
          enddo
        enddo
        !$acc end kernels
      enddo
#else
      call MPI_SENDRECV(p(1      ,1-nh,1-nh),1,halo,left ,0, &
                        p(nx+1   ,1-nh,1-nh),1,halo,right,0, &
                        comm_cart,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(nx+1-nh,1-nh,1-nh),1,halo,right,0, &
                        p(1-nh   ,1-nh,1-nh),1,halo,left ,0, &
                        comm_cart,MPI_STATUS_IGNORE,ierr)
#endif
    case(2) ! y direction
#if defined(_OPENACC)
      !$acc kernels present(ysl_buf,yrl_buf,ysr_buf,yrr_buf)
      ysl_buf(:,:) = - huge(1._rp)
      yrl_buf(:,:) = - huge(1._rp)
      ysr_buf(:,:) = - huge(1._rp)
      yrr_buf(:,:) = - huge(1._rp)
      !$acc end kernels
      !
      do q=0,nh-1
        !$acc kernels present(ysl_buf,ysr_buf)
        do k=1-nh,nz+nh 
          do i=1-nh,nx+nh
            ysl_buf(i,k) = p(i, 1+q,k)
            ysr_buf(i,k) = p(i,ny-q,k)
          enddo
        enddo
        !$acc end kernels
        !
        lb1 = lbound(ysl_buf,1) 
        lb2 = lbound(ysl_buf,2) 
        call MPI_SENDRECV(ysl_buf(lb1,lb2), size( ysl_buf ),MPI_REAL_RP,front,0, &
                          yrr_buf(lb1,lb2), size( yrr_buf ),MPI_REAL_RP,back ,0, &
                          comm_cart,status,ierr)
        call MPI_SENDRECV(ysr_buf(lb1,lb2), size( ysr_buf ),MPI_REAL_RP,back ,0, &
                          yrl_buf(lb1,lb2), size( yrl_buf ),MPI_REAL_RP,front,0, &
                          comm_cart,status,ierr)
        !
        !$acc kernels present(yrr_buf,yrl_buf)
        do k=1-nh,nz+nh 
          do i=1-nh,nx+nh
            p(i,ny+1+q,k) = yrr_buf(i,k)
            p(i,   0-q,k) = yrl_buf(i,k)
          enddo
        enddo
        !$acc end kernels
      enddo
#else
      call MPI_SENDRECV(p(1-nh,1      ,1-nh),1,halo,front,0, &
                        p(1-nh,ny+1   ,1-nh),1,halo,back ,0, &
                        comm_cart,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(1-nh,ny+1-nh,1-nh),1,halo,back ,0, &
                        p(1-nh,1-nh   ,1-nh),1,halo,front,0, &
                        comm_cart,MPI_STATUS_IGNORE,ierr)
#endif
    case(3) ! z direction
#if defined(_OPENACC)
      !$acc kernels present(zsl_buf,zrl_buf,zsr_buf,zrr_buf)
      zsl_buf(:,:) = - huge(1._rp)
      zrl_buf(:,:) = - huge(1._rp)
      zsr_buf(:,:) = - huge(1._rp)
      zrr_buf(:,:) = - huge(1._rp)
      !$acc end kernels
      !
      do q=0,nh-1
        !$acc kernels present(zsl_buf,zsr_buf)
        do j=1-nh,ny+nh 
          do i=1-nh,nx+nh
            zsl_buf(i,j) = p(i,j,  1+q)
            zsr_buf(i,j) = p(i,j, nz-q)
          enddo
        enddo
        !$acc end kernels
        !
        lb1 = lbound(zsl_buf,1) 
        lb2 = lbound(zsl_buf,2) 
        call MPI_SENDRECV(zsl_buf(lb1,lb2), size( zsl_buf ),MPI_REAL_RP,bottom,0, &
                          zrr_buf(lb1,lb2), size( zrr_buf ),MPI_REAL_RP,top   ,0, &
                          comm_cart,status,ierr)
        call MPI_SENDRECV(zsr_buf(lb1,lb2), size( zsr_buf ),MPI_REAL_RP,top   ,0, &
                          zrl_buf(lb1,lb2), size( zrl_buf ),MPI_REAL_RP,bottom,0, &
                          comm_cart,status,ierr)
        !
        !$acc kernels present(zrr_buf,zrl_buf)
        do j=1-nh,ny+nh 
          do i=1-nh,nx+nh
            p(i,j,nz+1+q) = zrr_buf(i,j)
            p(i,j,   0-q) = zrl_buf(i,j)
          enddo
        enddo
        !$acc end kernels
      enddo
#else
     call MPI_SENDRECV(p(1-nh,1-nh   ,1   ),1,halo,bottom,0, &
                       p(1-nh,1-nh   ,nz+1),1,halo,top   ,0, &
                       comm_cart,MPI_STATUS_IGNORE,ierr)
     call MPI_SENDRECV(p(1-nh,1-nh,nz+1-nh),1,halo,top   ,0, &
                       p(1-nh,1-nh   ,1-nh),1,halo,bottom,0, &
                       comm_cart,MPI_STATUS_IGNORE,ierr)
#endif
    !
    end select
    !
    return
    !
  end subroutine updthalo
  !
end module mod_bound
