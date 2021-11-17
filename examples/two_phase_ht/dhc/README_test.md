This case helps to validate the code with respect to single phase, buoyancy driven flows, using the Boussinesq approximation. The setup is similar to the case with Rayleigh number 10^6 in:
"Accurate solutions to the square thermally driven cavity at high Rayleigh number"
P. Le Quere, Computers and Fluids Vol.20, No. 1, pp 29-41, 1991.

Observables:
1. First column: time;
2. Second column: Nusselt on the hot wall;
3. Third column: Nusselt on the cold wall.

To visualize using e.g., gnuplot:
plot "nusselt_ref.out" using 1:2, "nusselt_ref.out" using 1:($3*(-1))
