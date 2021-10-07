module mod_fftw_spectra
! use mod_param
! use mod_common_post
use mod_common_mpi
use mod_types
! use mpi
implicit none
!
! FFTW 
!
integer FFTW_FORWARD,FFTW_BACKWARD
parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

integer FFTW_ESTIMATE,FFTW_MEASURE
parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)

integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
parameter (FFTW_OUT_OF_PLACE=0)
parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

integer FFTW_THREADSAFE
parameter (FFTW_THREADSAFE=128)

private
public init_fft,fftr2c,fftc2r!,clean_fft
contains

! subroutine init_fft(nx,ny,nz)
! 	implicit none
! 	integer, intent(in) :: nx,ny,nz
! 	real(rp), dimension(nx) :: varxr
! 	real(rp), dimension(ny) :: varyr
! 	real(rp), dimension(nz) :: varzr
! 	real(rp), dimension(nx,ny,nz) :: var
! 	complex(rp), dimension(nx/2+1,ny,nz) :: varc
! 	complex(rp), dimension(nx/2+1) :: varxc
! 	complex(rp), dimension(ny/2+1) :: varyc
! 	complex(rp), dimension(nz/2+1) :: varzc
! 	call dfftw_plan_dft_r2c_1d(plan_x_d,nx,varxr,varxc,FFTW_ESTIMATE)
! 	call dfftw_plan_dft_c2r_1d(plan_x_i,nx,varxc,varxr,FFTW_ESTIMATE)
! 	call dfftw_plan_dft_r2c_1d(plan_y_d,ny,varyr,varyc,FFTW_ESTIMATE)
! 	call dfftw_plan_dft_c2r_1d(plan_y_i,ny,varyc,varyr,FFTW_ESTIMATE)
! 	call dfftw_plan_dft_r2c_1d(plan_z_d,nz,varzr,varzc,FFTW_ESTIMATE)
! 	call dfftw_plan_dft_c2r_1d(plan_z_d,nz,varzc,varzr,FFTW_ESTIMATE)
! 	!
! 	! call dfftw_plan_dft_1d(plan_x_d_c2c,nx/2+1,varxc,varxc,FFTW_FORWARD,FFTW_ESTIMATE)
! 	! call dfftw_plan_dft_1d(plan_x_i_c2c,nx/2+1,varxc,varxc,FFTW_BACKWARD,FFTW_ESTIMATE)
! 	! call dfftw_plan_dft_1d(plan_y_d_c2c,ny/2+1,varyc,varyc,FFTW_FORWARD,FFTW_ESTIMATE)
! 	! call dfftw_plan_dft_1d(plan_y_i_c2c,ny/2+1,varyc,varyc,FFTW_BACKWARD,FFTW_ESTIMATE)
! 	! call dfftw_plan_dft_1d(plan_z_d_c2c,nz/2+1,varzc,varzc,FFTW_FORWARD,FFTW_ESTIMATE)
! 	! call dfftw_plan_dft_1d(plan_z_i_c2c,nz/2+1,varzc,varzc,FFTW_BACKWARD,FFTW_ESTIMATE)
! 	! call dfftw_plan_dft_c2r_3d(plan_3d_r2c,nx,ny,nz,var,varc,FFTW_ESTIMATE)
! 	! call dfftw_plan_dft_r2c_3d(plan_3d_c2r,nx,ny,nz,varc,var,FFTW_ESTIMATE)

! end subroutine init_fft
subroutine init_fft(n, plan_r2c, plan_c2r)
implicit none
integer, intent(in) :: n
integer(8), intent(out) :: plan_r2c, plan_c2r
real, dimension(n) :: var
complex, dimension(n/2+1) :: varc
!
call dfftw_plan_dft_r2c_1d(plan_r2c,n,var,varc,FFTW_ESTIMATE)
call dfftw_plan_dft_c2r_1d(plan_c2r,n,varc,var,FFTW_ESTIMATE)
!
return
end subroutine init_fft



subroutine fftr2c(n,varin,varout,plan_r2c)
implicit none
integer, intent(in) :: n
integer(8), intent(in) ::plan_r2c
real(rp), intent(in), dimension(n) :: varin
complex(rp), intent(out), dimension(n/2+1) :: varout
!
call dfftw_execute_dft_r2c(plan_r2c,varin,varout)
!
return
end subroutine fftr2c
!
subroutine fftc2r(n,varin,varout,plan_c2r)
implicit none
integer, intent(in) :: n
integer(8), intent(in) :: plan_c2r
complex(rp), intent(in), dimension(n/2+1) :: varin
real(rp), intent(out), dimension(n) :: varout
!
call dfftw_execute_dft_c2r(plan_c2r,varin,varout)
varout(:) = varout(:)/n
!
return
end subroutine fftc2r

! subroutine clean_fft
! implicit none
! !
! call dfftw_destroy_plan(plan_r2c)
! call dfftw_destroy_plan(plan_c2r)
! !
! return
! end subroutine clean_fft

end module mod_fftw_spectra
