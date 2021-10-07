### Enabling input files
Depending on the pre-processor flags considered different input files are required. In the following, descriptions of the input files are provided. NOTE: line position is important in input files, hence use the profided templates.

***dns.in***, always available

* `itot,jtot,ktot` : grid points in each direction
* `lx,ly,lz`       : domain dimensions in each direction
* `gr`             : stretching parameter (UNTESTED)
* `cfl,dt_input`   : cfl and fixed-dt used
* `constant_dt`    : whether to use constant CFL or constant dt
* `visc`           : fluid viscosity (active only if single-phase)
* `inivel`         : type of initialized velocity field
* `is_wallturb`    : initializes velocity conditions for faster turbulence transition in wall-bounded flows
* `nstep, time_max,tw_max`   : stopping criteria, i.e. maximum number of timestep, maximum simulation time, maximum wall-time
* `stop_type(1),stop_type(2),stop_type(3)`  : enables stopping criteria for the simulation
* `restart`                               : restart simulation from last dump file
* `icheck,iout0d,iout1d,iout2d,iout3d,isave` : set number of timesteps between saved output
* `cbcvel(0,1,1),cbcvel(1,1,1),cbcvel(0,2,1),cbcvel(1,2,1),cbcvel(0,3,1),cbcvel(1,3,1)`: U velocity BC type
* `cbcvel(0,1,2),cbcvel(1,1,2),cbcvel(0,2,2),cbcvel(1,2,2),cbcvel(0,3,2),cbcvel(1,3,2)`: V velocity BC type
* `cbcvel(0,1,3),cbcvel(1,1,3),cbcvel(0,2,3),cbcvel(1,2,3),cbcvel(0,3,3),cbcvel(1,3,3)`: W velocity BC type
* `cbcpre(0,1  ),cbcpre(1,1  ),cbcpre(0,2  ),cbcpre(1,2  ),cbcpre(0,3  ),cbcpre(1,3  )`: pressure BC type
* ` bcvel(0,1,1), bcvel(1,1,1), bcvel(0,2,1), bcvel(1,2,1), bcvel(0,3,1), bcvel(1,3,1)`: U velocity BC value
* ` bcvel(0,1,2), bcvel(1,1,2), bcvel(0,2,2), bcvel(1,2,2), bcvel(0,3,2), bcvel(1,3,2)`: V velocity BC value
* ` bcvel(0,1,3), bcvel(1,1,3), bcvel(0,2,3), bcvel(1,2,3), bcvel(0,3,3), bcvel(1,3,3)`: W velocity BC value
* ` bcpre(0,1  ), bcpre(1,1  ), bcpre(0,2  ), bcpre(1,2  ), bcpre(0,3  ), bcpre(1,3  )`: pressure BC value
* ` bforce(1),bforce(2),bforce(3)` : body force source term
* ` is_forced(1),is_forced(2),is_forced(3)` : if forced mean velocity in one direction 
* ` bulk_velx,bulk_vely,bulk_velz` : value of the forced velocity
* ` is_outflow(0,1),is_outflow(1,1),is_outflow(0,2),is_outflow(1,2),is_outflow(0,3),is_outflow(1,3)` : defines in outflow BC are used
* `dims_in(1),dims_in(2)`: number of cores/GPU per direction (note that GPU should be parallelized in slabs) 
* `nthreadsmax`: maximum number of threads for OpenMP (UNTESTED)       

***vof.in***, enabled when `_USE_VOF` is used
* `rho1, rho2, mu1,mu2` : density and dynamic viscosity for fluid 1 and 2
* `inivof` : type of vof initialization
* `nbub` : number of bubbles/droplet initialized if `inivof='bub'`. Phase 1 is initialized inside the dispersed phase. If `nbub>1` the file bub.in is used, where all the coordinates and radii are set (see line below)
* `xc(1), yc(1), zc(1), r(1)` : bubble/droplet center and radius
* `cbcvof(0,1 ),cbcvof(1,1 ),cbcvof(0,2 ),cbcvof(1,2 ),cbcvof(0,3 ),cbcvof(1,3 )` : VoF BC type
* `bcvof(0,1  ),bcvof(1,1  ),bcvof(0,2  ),bcvof(1,2  ),bcvof(0,3  ),bcvof(1,3  )` : VoF BC values
* `gacc(1), gacc(2), gacc(3)` : gravity vector
* `sigma` : surface tension coefficient
* `dpdl(1), dpdl(2), dpdl(3)` : imposed pressure gradient
* `late_init, i_late_init` : initialized the VoF field later at timestep `i_late_init`

***bub.in*** enabled when `_USE_VOF` is used, `inivof='bub'` and  `bub>1` in `vof.in`
* `xc(1), yc(1), zc(1), r(1)` : center and radius for the i-th bubble/droplet. This line should be repeated for as many elements as `nbub`

***forcing.in*** enabled when `_TURB_FORCING` is used
* `turbType` : type of turbulence forcing. Now ABC and TGV are available
* `u0_t`     : initial velocity magnitude
* `f0_t`     : intensity of the forcing term
* `k0_t`     : wavenumber for energy injection
* `abc(1), abc(2), abc(3)` : values of A, B and C coefficients for ABC forcing
* `add_noise_abc` : add disturbances to the initial condition to trigger transition

***heat_transfer.in*** enabled if `_HEAT_TRANSFER` is used
* `initmp` : type of temperature field initialization
* `tl0,tg0`: initial liquid and gas temperature
* `cp1,cp2`: cp of phase 1 and 2
* `cv1,cv2`: cv pf phase 1 and 2
* `kappa1,kappa2` : thermal conductivity of phase 1 and 2
* `cbctmp(0,1 ),cbctmp(1,1 ),cbctmp(0,2 ),cbctmp(1,2 ),cbctmp(0,3 ),cbctmp(1,3 )`: temperature BC type
* `bctmp(0,1  ),bctmp(1,1  ),bctmp(0,2  ),bctmp(1,2  ),bctmp(0,3  ),bctmp(1,3  )`: temperature BC value
* `if _BOUSSINESQ`, `tmp0,beta_th`: if `_BOUSSINESQ` is enabled it uses values from these line, which are the arithmetic mean temperature and thermal expansion coeffient 


