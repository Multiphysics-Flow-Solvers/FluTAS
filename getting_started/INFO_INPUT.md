### Enabling input files
Depending on the pre-processor flags considered different input files are required. In the following, descriptions of the input files are provided. NOTE: line position is important in input files, hence use and, if necessary, modify the provided templates available in the [`examples`](./../examples).

The input file ***dns.in*** must be always provide to run a simulation. It contains:
* `itot,jtot,ktot` : grid points in each direction
* `lx,ly,lz`       : domain dimensions in each direction
* `gr`             : stretching parameter (SUPPORTED only in single-phase)
* `cfl,dt_input`   : cfl and employed constant time-step, dt (the constant dt is used only if `constant_dt` is true)
* `constant_dt`    : prescribe a constant dt (`T`) or constant CFL (`F`)
* `time_scheme`, `space_scheme_mom`: scheme for the time and space discretization of the momentum equation
* `rho_sp`,`mu_sp` : fluid density and dynamic viscosity (relevant only single-phase, i.e.,`!defined(_USE_VOF)`. If the flag `defined(_USE_VOF)` is employed during the compilation, the values of `rho_sp` and `mu_sp` are overwritten by the corresponding ones defined in `vof.in` for phase 2)
* `inivel`         : type of initialized velocity field
* `is_wallturb`    : initializes velocity conditions for faster turbulence transition in wall-bounded flows
* `wallturb_type`  : initial condition for triggering turbulent transition in wall-bounded flows
* `bulk_ftype`     : type of forcing to sustain the flow in canonical wall-bounded domains
* `nstep,time_max,tw_max` : stopping criteria, i.e. maximum number of time-steps, maximum simulation time, maximum wall-time
* `stop_type(1),stop_type(2),stop_type(3)` : enables stopping criteria for the simulation
* `restart` : restart simulation from last dump file (restarting files are printed in `data/restart_dir/`)
* `icheck,iout0d,iout1d,iout2d,iout3d,isave` : set the time-step frequency employed to: i) perform numerical stability checks, i.e., time-step restriction and velocity divergence (`icheck`), ii) print zero, one, two, three dimensional outputs (`iout0d,iout1d,iout2d,iout3d`) and iii) save the restarting files (`isave`)
* `cbcvel(0,1,1),cbcvel(1,1,1),cbcvel(0,2,1),cbcvel(1,2,1),cbcvel(0,3,1),cbcvel(1,3,1)`: U velocity BC type
* `cbcvel(0,1,2),cbcvel(1,1,2),cbcvel(0,2,2),cbcvel(1,2,2),cbcvel(0,3,2),cbcvel(1,3,2)`: V velocity BC type
* `cbcvel(0,1,3),cbcvel(1,1,3),cbcvel(0,2,3),cbcvel(1,2,3),cbcvel(0,3,3),cbcvel(1,3,3)`: W velocity BC type
* `cbcpre(0,1  ),cbcpre(1,1  ),cbcpre(0,2  ),cbcpre(1,2  ),cbcpre(0,3  ),cbcpre(1,3  )`: pressure BC type
* ` bcvel(0,1,1), bcvel(1,1,1), bcvel(0,2,1), bcvel(1,2,1), bcvel(0,3,1), bcvel(1,3,1)`: U velocity BC value
* ` bcvel(0,1,2), bcvel(1,1,2), bcvel(0,2,2), bcvel(1,2,2), bcvel(0,3,2), bcvel(1,3,2)`: V velocity BC value
* ` bcvel(0,1,3), bcvel(1,1,3), bcvel(0,2,3), bcvel(1,2,3), bcvel(0,3,3), bcvel(1,3,3)`: W velocity BC value
* ` bcpre(0,1  ), bcpre(1,1  ), bcpre(0,2  ), bcpre(1,2  ), bcpre(0,3  ), bcpre(1,3  )`: pressure BC value
* ` is_forced(1),is_forced(2),is_forced(3)` : choose the direction along which a velocity (`bulk_ftype='cfr'`) or pressure gradient (`bulk_ftype='cpg'`) is imposed to sustain the flow 
* ` gacc_x,gacc_y,gacc_z` : gravity acceleration
* ` bvel_x,bvel_y,bvel_z` : value of the imposed velocity
* ` dpdl_x,dpdl_y,dpdl_z` : value of the imposed pressure gradient
* ` is_outflow(0,1),is_outflow(1,1),is_outflow(0,2),is_outflow(1,2),is_outflow(0,3),is_outflow(1,3)` : set if the physical boundary is an outflow (`T`) or not (`F`)
* `dims_in(1),dims_in(2)`: number of cores/GPU per direction (note that GPU should be parallelized in slabs)
* `nthreadsmax`: maximum number of threads for OpenMP (UNTESTED)       

The input file ***vof.in*** is required when `_USE_VOF` is used. It contains:
* `rho1,rho2,mu1,mu2` : density and dynamic viscosity for phase 1 and 2
* `inivof` : type of vof initialization
* `nbub` : number of bubbles/droplet initialized if `inivof='bub'`. Note that if `nbub>1` the file bub.in is required with the center coordinates and radii of each emulsion, droplet or bubble. Otherwise, if `nbub=1`, this information is provided in the next line.
* `xc(1), yc(1), zc(1), r(1)` : center and radius of the single emulsion, droplet or bubble present in the domain
* `cbcvof(0,1 ),cbcvof(1,1 ),cbcvof(0,2 ),cbcvof(1,2 ),cbcvof(0,3 ),cbcvof(1,3 )` : VoF BC type
* `bcvof(0,1  ),bcvof(1,1  ),bcvof(0,2  ),bcvof(1,2  ),bcvof(0,3  ),bcvof(1,3  )` : VoF BC values
* `sigma` : surface tension coefficient
* `late_init, i_late_init` : initialized the VoF field later at time-step `i_late_init`

Note that by convention phase 1 and 2 represent the dispersed and the continous phase, respectively. 

The input file ***bub.in*** is required when `_USE_VOF` is used, `inivof='bub'` and  `bub>1` in `vof.in`. It contains:
* `xc(i), yc(i), zc(i), r(i)` : center and radius for the i-th bubble/droplet. This line should be repeated for as many elements as `nbub`. The user can either manually provide the center coordinates and radius of each droplet or, alternatively, generate a random distribution of droplets with equal radius using the script `randomDroplet.py` available [`here`](./../utils/preprocessing/).

The input file ***forcing.in*** is required when `_TURB_FORCING` is used. It contains:
* `turb_type`: type of turbulence forcing. Currently available: Arnold-Beltrami-Childress ('abc') and Taylor-Green Vortex ('tgv')
* `u0_t`     : initial velocity magnitude
* `f0_t`     : intensity of the forcing term
* `k0_t`     : wavenumber at which energy is injected
* `abc_x, abc_y, abc_z` : values of A, B and C coefficients for ABC forcing
* `add_noise_abc` : add disturbances to the initial condition to trigger transition

The input file ***heat_transfer.in*** is required if `_HEAT_TRANSFER` is used. It contains:
* `initmp` : type of temperature field initialization
* `tl0,tg0`: initial liquid and gas temperature
* `cp1,cp2`: specific heat capacity at constant pressure of phase 1 and 2
* `cv1,cv2`: specific heat capacity at constant volume of phase 1 and 2
* `kappa1,kappa2` : thermal conductivity of phase 1 and 2
* `cbctmp(0,1 ),cbctmp(1,1 ),cbctmp(0,2 ),cbctmp(1,2 ),cbctmp(0,3 ),cbctmp(1,3 )`: temperature BC type
* `bctmp(0,1  ),bctmp(1,1  ),bctmp(0,2  ),bctmp(1,2  ),bctmp(0,3  ),bctmp(1,3  )`: temperature BC value
* `if _BOUSSINESQ`, `tmp0,beta1_th,beta2_th`: if `_BOUSSINESQ` is enabled it uses values from these line, which are the arithmetic mean temperature and the thermal expansion coefficient of both phases

Note that in case the app `single_phase` is chosen and heat transfer effects are included, i.e.,`HEAT_TRANSFER=1`, the input file ***heat_transfer_sp.in*** becomes:

* `initmp` : type of temperature field initialization
* `tmp0`: initial temperature
* `cp_sp`: specific heat capacity at constant pressure
* `cv_sp`: specific heat capacity at constant volume
* `kappa_sp` : thermal conductivity
* `cbctmp(0,1 ),cbctmp(1,1 ),cbctmp(0,2 ),cbctmp(1,2 ),cbctmp(0,3 ),cbctmp(1,3 )`: temperature BC type
* `bctmp(0,1  ),bctmp(1,1  ),bctmp(0,2  ),bctmp(1,2  ),bctmp(0,3  ),bctmp(1,3  )`: temperature BC value
* `beta_th`: if `_BOUSSINESQ` is enabled, the code requires also to provide the thermal expansion coefficient

