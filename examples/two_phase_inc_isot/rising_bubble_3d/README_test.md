This case helps to validate the code with respect to two-phase incompressible and isothermal code. The setup is the Rising Bubble test case: "Hysing, S.; Turek, S.; Kuzmin, D.; Parolini, N.; Burman, E.; Ganesan, S.; Tobiska, L.: Quantitative benchmark computations of two-dimensional bubble dynamics, International Journal for Numerical Methods in Fluids, Volume 60 Issue 11, Pages 1259-1288, DOI: 10.1002/fld.1934, 2009".

Observables:
1. First column: time;
2. Second column: x-position of bubble center of mass;
3. Third column: y-position of bubble center of mass;
4. Fourth column: z-position of bubble center of mass;
5. Fifth column: x-velocity of the bubble center of mass;
6. Sixth column: y-velocity of the bubble center of mass;
7. Seventh column: z-velocity of the bubble center of mass;

We are interested in z-position and z-velocity of the bubble center of mass. To visualize them:

  * for the position: plot "pos_vt_ref.out" using 1:4
  * for the velocity: plot "pos_vt_ref.out" using 1:7
