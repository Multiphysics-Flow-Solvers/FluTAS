### Enabling input files
Depending on the pre-processor flags considered, different input files are required. In the following, descriptions of the input files are provided. Important premises:
* **FluTAS** strictly employs dimensional inputs based on the [SI units](https://en.wikipedia.org/wiki/International_System_of_Units). The user should choose the input (physical) quantities to match the desired dimensionless parameters of the study. 
    1. **Example A**: for a channel flow configuration, the user should input the forcing bulk velocity (or the imposed pressure gradient) $U_b$, the height of the channel $l_z$, the dynamic viscosity $\mu$, and the density of the fluid $\rho$ to match the desired Reynolds number, i.e. $Re = \rho U_bl_z/\mu$;
    2. **Example B**: for a heat transfer problem, the user should input the dynamic viscosity $\mu$, the specific heat capacity $c_p$, and the thermal conductivity $k$ to match the desired Prandtl number, i.e. $Pr=\mu c_p/k$.

  The above reasoning applies to any other physical dimensionless parameter of the problem under consideration.
* Line position is essential when constructing the input files. As a reference, use and, if necessary, modify the provided templates available in the [`examples`](./../examples) folder.

The input file ***dns.in*** must be always provide to run a simulation. It contains:
* `itot,jtot,ktot` : grid points in each direction
* `lx,ly,lz`       : domain dimensions in each direction
* `gr`             : stretching parameter (SUPPORTED only in single-phase)
* `cfl,dt_input`   : CFL and employed constant time-step. Note that `dt_input` is used only if `constant_dt` is true
* `constant_dt`    : prescribe a constant dt (`T`) or constant CFL (`F`)
* `time_scheme`, `space_scheme_mom`: scheme for the time and space discretization of the momentum equation
* `rho_sp`,`mu_sp` : density and dynamic viscosity of the fluid (relevant only single-phase, i.e.,`!defined(_USE_VOF)`. If the flag `_USE_VOF` is employed during the compilation, the values of `rho_sp` and `mu_sp` are overwritten by the corresponding ones defined in `vof.in` for phase 2)
* `inivel`, `is_noise_vel`, `noise_vel`: type of initialized velocity field. If the logical variable `is_noise_vel` is true, a random noise of absolute maximum magnitude `noise_vel` is superimposed on the initial velocity field
* `is_wallturb`    : initializes velocity conditions for faster turbulence transition in wall-bounded flows
* `wallturb_type`  : initial condition for triggering turbulent transition in wall-bounded flows
* `bulk_ftype`     : type of forcing to sustain the flow in canonical wall-bounded domains
* `nstep,time_max,tw_max` : stopping criteria, i.e. maximum number of time-steps, maximum simulation time, maximum wall-time
* `stop_type(1),stop_type(2),stop_type(3)` : enables stopping criteria for the simulation
* `restart, num_max_chkpt, input_chkpt, latest` : `restart` is a logical variable, which is true when the user wants to restart the simulation from a dump file (restarting files are printed in `data/restart_dir/restart_subdir_???`). `num_max_chkpt` is the maximum number of saved checkpoint files for restarting. `input_chkpt` is a user-chosen checkpoint file from which the simulation should be restarted (if available). If the logical variable `latest` is true, the simulation restarts from the latest available/saved checkpoint and, therefore, the variable `input_chkpt` is not used
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
 * `dims_in(1),dims_in(2)`: number of CPU cores or GPUs per parallelized direction. Note that **FluTAS** employs a two-dimensional parallelization with pencils aligned along the non-parallelized direction. Three options are available:
   1. Pencils oriented along the **x** direction. The pre-processor flag `_DECOMP_X` controls this option. In this case, a number of `dims_in(1)` and `dims_in(2)` processors are distributed along the **y** and the **z** directions, respectively. This option should be typically preferred both for CPUs and GPUs runs since it minimizes the number of `all-to-all` operations in the Poisson solver. For GPUs runs, the option `_DECOMP_X` allows to employ a slab-decomposition only;
   2. Pencils oriented along the **y** direction. The pre-processor flag `_DECOMP_Y` controls this option. In this case, a number of `dims_in(1)` and `dims_in(2)` processors are distributed along the **x** and the **z** directions, respectively. For GPUs runs, the option `_DECOMP_X` allows to employ a slab-decomposition only;
   3. Pencils oriented along the **z** direction. The pre-processor flag `_DECOMP_Z` controls this option. In this case, a number of `dims_in(1)` and `dims_in(2)` processors are distributed along the **x** and the **y** directions, respectively. For GPUs runs, the option `_DECOMP_Z` allows both a slab and a pencil decomposition.
* `nthreadsmax`: maximum number of threads for OpenMP (UNTESTED)       

The input file ***vof.in*** is required when `_USE_VOF` is used. It contains:
* `rho1,rho2` and `mu1,mu2` : density and dynamic viscosity for phase 1 and 2
* `inivof` : type of vof initialization
* `nbub` : number of bubbles/droplet initialized if `inivof='bub'`. Note that if `nbub>1` the file bub.in is required with the center coordinates and radii of each emulsion, droplet or bubble. Otherwise, if `nbub=1`, this information is provided in the next line.
* `xc(1), yc(1), zc(1), r(1)` : center and radius of the single emulsion, droplet or bubble present in the domain
* `cbcvof(0,1 ),cbcvof(1,1 ),cbcvof(0,2 ),cbcvof(1,2 ),cbcvof(0,3 ),cbcvof(1,3 )` : VoF BC type
* `bcvof(0,1  ),bcvof(1,1  ),bcvof(0,2  ),bcvof(1,2  ),bcvof(0,3  ),bcvof(1,3  )` : VoF BC values
* `sigma` : surface tension coefficient
* `late_init, i_late_init` : if `late_init` is true, the VoF field is initialized at time-step `i_late_init` instead of at `istep=0`. This option is typically useful when the user wants to run a precursor single-phase simulation and add the disperse phase at later stage.

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
* `initmp`, `is_noise_tmp`, `noise_tmp`: type of initialized temperature field. If the logical variable `is_noise_tmp` is true, a random noise of absolute maximum magnitude `noise_tmp` is superimposed on the initial temperature field
* `tl0,tg0`: initial liquid and gas temperature
* `cp1,cp2`: specific heat capacity at constant pressure of phase 1 and 2
* `cv1,cv2`: specific heat capacity at constant volume of phase 1 and 2
* `kappa1,kappa2` : thermal conductivity of phase 1 and 2
* `cbctmp(0,1 ),cbctmp(1,1 ),cbctmp(0,2 ),cbctmp(1,2 ),cbctmp(0,3 ),cbctmp(1,3 )`: temperature BC type
* `bctmp(0,1  ),bctmp(1,1  ),bctmp(0,2  ),bctmp(1,2  ),bctmp(0,3  ),bctmp(1,3  )`: temperature BC value
* `tmp0,beta1_th,beta2_th`: if `_BOUSSINESQ` is enabled, the code requires also to provide a reference temperature (e.g. the arithmetic mean between top and bottom wall), the thermal expansion coefficient for phase 1 and phase 2

Note that in case the app `single_phase` is chosen and heat transfer effects are included, i.e.,`HEAT_TRANSFER=1`, the input file ***heat_transfer_sp.in*** must be provided together with ***dns.in***. In particular, ***heat_transfer_sp.in*** contains:
* `initmp`, `is_noise_tmp`, `noise_tmp`: type of initialized temperature field. If the logical variable `is_noise_tmp` is true, a random noise of absolute maximum magnitude `noise_tmp` is superimposed on the initial temperature field
* `tmp0`: initial temperature
* `cp_sp`: specific heat capacity at constant pressure
* `cv_sp`: specific heat capacity at constant volume
* `kappa_sp` : thermal conductivity of the fluid
* `cbctmp(0,1 ),cbctmp(1,1 ),cbctmp(0,2 ),cbctmp(1,2 ),cbctmp(0,3 ),cbctmp(1,3 )`: temperature BC type
* `bctmp(0,1  ),bctmp(1,1  ),bctmp(0,2  ),bctmp(1,2  ),bctmp(0,3  ),bctmp(1,3  )`: temperature BC value
* `beta_th`: if `_BOUSSINESQ` is enabled, the code requires also to provide the thermal expansion coefficient of the fluid

