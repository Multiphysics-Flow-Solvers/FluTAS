module mod_common_post
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

	! integer:: plan_x_d, plan_x_i, plan_y_d, plan_y_i, plan_z_d, plan_z_i
	! integer :: plan_r2c, plan_c2r
end module mod_common_post