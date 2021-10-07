#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

struct sphere_params {double xc; double yc; double zc; double r; double beta; 
double h;};

double
my_f (double x[], size_t dim, void * p) {
   struct sphere_params * fp = (struct sphere_params *)p;
   double psi;
   // double beta = 1.0;
   // double h = 0.005;

   if (dim != 3)
      {
        fprintf (stderr, "error: dim != 3");
        abort ();
      }
   psi = fp->r - sqrt(
            pow((x[0]-fp->xc),2)
            + pow((x[1]-fp->yc),2)
            + pow((x[2]-fp->zc),2));

  return 0.5*(1+tanh(fp->beta/fp->h*psi));
}

// int
// main (void)
extern "C" void init_MTHINC_(double *xc, double *yc, double *zc, double *rs,
  double *xl, double *xu, double *dx, double *beta, double *vof)
{
  double err;

  // double xl[3] = { 0, 0, 0 };
  // double xu[3] = { 1, 1, 1 };

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function F;
  struct sphere_params params = { *xc, *yc, *zc, *rs, *beta, *dx};

  F.f = &my_f;
  F.dim = 3;
  F.params = &params;

  size_t calls = 500000;

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  // {
  //   gsl_monte_plain_state *s = gsl_monte_plain_alloc (3);
  //   gsl_monte_plain_integrate (&F, xl, xu, 3, calls, r, s, 
  //                              &res, &err);
  //   gsl_monte_plain_free (s);
  //   printf(" % .6f\n",res );

  //   // display_results ("plain", res, err);
  // }

  {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
    gsl_monte_miser_integrate (&F, xl, xu, 3, calls, r, s,
                               &*vof, &err);
    gsl_monte_miser_free (s);

    // display_results ("miser", res, err);
    // printf(" % .6f\n",*vof );
  }

  // {
  //   gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

  //   gsl_monte_vegas_integrate (&F, xl, xu, 3, 10000, r, s,
  //                              &res, &err);
  //   display_results ("vegas warm-up", res, err);

  //   printf ("converging...\n");
  //   printf(" % .6f\n",res );

  //   do
  //     {
  //       gsl_monte_vegas_integrate (&F, xl, xu, 3, calls/5, r, s,
  //                                  &res, &err);
  //       // printf ("result = % .6f sigma = % .6f "
  //               "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
  //     }
  // //   while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

  // //   // display_results ("vegas final", res, err);

  // //   gsl_monte_vegas_free (s);
  // // }

  // gsl_rng_free (r);

  // return 0;
}
