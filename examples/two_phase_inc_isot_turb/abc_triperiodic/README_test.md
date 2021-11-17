This case helps to validate the code with respect to single phase incompressible and isothermal code. The setup is the Arnold-Beltrami-Childress (ABC) flow as described in e.g.:
"Podvigina, O & Pouquet, A1994 On the non-linear stability of the 1: 1: 1 abc flow.PhysicaD: Nonlinear Phenomena75(4), 471–508."

Observables:
1. First column: time-steps;
2. Second column: time;
3. Third column: mean kinetic energy;

We are interested in the evolution of the mean kinetic energy in time, which, given the laminar regime, should be kept as close as possible to 1.5. To visualize it:
plot "ke_t_ref.out" using 2:3
